#include "simulation_thread.h"
#include <qmessagebox.h>
#include "utility.h"

simulation_thread::simulation_thread(QObject* parent)
	: QThread(parent)
{
}

simulation_thread::~simulation_thread()
{
}


void simulation_thread::run()
{
	stop = false;
	std::ofstream out;
	if (!CUtility::file_exists(filename))
		out.open(filename.c_str());
	else
	{
		std::string errormessage = "File " + filename + " already exists.";
		emit ErrorMessage(errormessage.c_str());
		return;
	}

	//initialize simulation class
	ridesharing_sim sim(number_of_nodes, number_of_buses, bus_type, 0);
	sim.init_network(number_of_buses, number_of_nodes, topology);

	if (number_of_bus_calculations > 1)
	{
		QVector<double> vB, vE;
		QVector<double> vx, vC, vB_equiv_x_inf, vB_equiv_x;
		QVector<QPair<ULL, double>> vBEdata;

		if (number_of_request_rates > 1)
		{
			// read pre-calculated data for efficiency as a function of B
			std::string line;
			std::ifstream myfile("torus_25_Efficiency_x_7.5.dateff.dat");
			bool NumberNext = false;
			if (myfile.is_open())
			{
				while (std::getline(myfile, line))
				{
					if (line != "")
					{
						// file structure : B <tab> E
						size_t findtab = line.find('\t');
						ULL B = std::stoi(line.substr(0, findtab));
						double E = std::stod(line.substr(findtab + 1));
						if (vBEdata.size() > 0)
						{
							// do not save repeated entries
							if (vBEdata.back().first != B)
								vBEdata.push_back(QPair<ULL, double>(B, E));
						}
						else
							vBEdata.push_back(QPair<ULL, double>(B, E));
					}
				}
			}
			else
			{
				std::string errormessage = "Could not open file torus_25_Efficiency_x_7.5.dateff.dat";
				emit ErrorMessage(errormessage.c_str());
				return;
			}
		}

		std::ofstream outplot;
		std::string plotfilename = filename + "eff.dat";
		if (!CUtility::file_exists(plotfilename))
			outplot.open(plotfilename.c_str());
		else
		{
			std::string errormessage = "File " + plotfilename + " already exists.";
			emit ErrorMessage(errormessage.c_str());
			return;
		}

		emit ProcessTextChanged("Starting simulation");
		emit GraphChanged(vB, vE);
		//ui.textEdit_6->setPlainText(ui.textEdit_6->toPlainText() + "\nStarting simulation\n");

		double index_low_B = log(number_of_buses_from);
		double index_high_B = log(number_of_buses_to);
		for (double iter_B = index_low_B; iter_B <= (index_high_B + 0.0001); iter_B += (index_high_B - index_low_B) / (number_of_bus_calculations - 1))
		{
			if (stop)
				break;
			number_of_buses = (int)(exp(iter_B));
			QString newText = "Simulating " + QString::fromStdString(std::to_string(number_of_buses)) + " buses\n";
			emit ProcessTextChanged(newText);
			emit GraphChanged(vB, vE);
			
			if (number_of_request_rates > 1)
			{
				double index_low_x = log(normalized_request_rate_from);
				double index_high_x = log(normalized_request_rate_to);
				for (double iter_x = index_low_x; iter_x <= (index_high_x + 0.0001); iter_x += (index_high_x - index_low_x) / (number_of_request_rates - 1))
				{
					if (stop)
						break;
					normalized_request_rate = exp(iter_x);
					QString newText = "Simulating B = " + QString::fromStdString(std::to_string(number_of_buses)) + ", x = " + QString::fromStdString(std::to_string(normalized_request_rate));
					emit ProcessTextChanged(newText);
					emit GraphChanged(vx, vC);
					//ui.textEdit_6->setPlainText(newText.c_str());

					sim.init_new_sim(number_of_buses, number_of_nodes, topology, normalized_request_rate, bus_type, 500 * number_of_buses);

					//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
					sim.enable_measurements((1.0 * number_of_buses) / sim.request_rate);
					//simulate (and measure) for 100 requests per bus (at least 10000 requests)
					sim.run_sim_requests(std::max((ULL)10000, 500 * number_of_buses));

					//output results
					if (save)
					{
						sim.print_params(out, true);
						sim.print_measurements(out);
						out << std::endl << std::endl;
					}

					vB.push_back((double) number_of_buses);
					double Eff = normalized_request_rate / sim.measurements.get_av_scheduled_customers();
					vE.push_back(Eff);
					vx.push_back((double)normalized_request_rate);
					vC.push_back(sim.measurements.get_av_scheduled_customers());
					// find equivalent number of buses in unlimited-capacity-simulation
					for (int iter_pre = 0; iter_pre < vBEdata.size(); iter_pre++)
					{
						if (iter_pre < vBEdata.size() - 1)
						{
							if (Eff > vBEdata[iter_pre].second && Eff < vBEdata[iter_pre + 1].second)
							{
								// linear interpolation
								double B_equiv = ((double)vBEdata[iter_pre].first) + ((double) (vBEdata[iter_pre + 1].first - vBEdata[iter_pre].first)) * (Eff - vBEdata[iter_pre].second) / (vBEdata[iter_pre + 1].second - vBEdata[iter_pre].second);
								vB_equiv_x_inf.push_back(B_equiv);
								break;
							}
							else if (Eff < vBEdata[iter_pre].second && Eff > vBEdata[iter_pre + 1].second)
							{
								vB_equiv_x_inf.push_back((vBEdata[iter_pre].second + vBEdata[iter_pre + 1].second) / 2);
								break;
							}
						}
						else
						{
							vB_equiv_x_inf.push_back((double) vBEdata[iter_pre].first);
						}
					}

					outplot << vB.back() << '\t' << vE.back() << '\t' << vx.back() << '\t' << vC.back() << '\t' << vB_equiv_x_inf.back() << '\n';
				}
			}
			else
			{
				sim.init_new_sim(number_of_buses, number_of_nodes, topology, normalized_request_rate, bus_type, 100 * number_of_buses);

				//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
				sim.enable_measurements((2.0 * number_of_buses) / sim.request_rate);
				//simulate (and measure) for 100 requests per bus (at least 10000 requests)
				sim.run_sim_requests(std::max((ULL)10000, 500 * number_of_buses));

				//output results
				if (save)
				{
					sim.print_params(out, true);
					sim.print_measurements(out);
					out << std::endl << std::endl;
				}

				vB.push_back((double)number_of_buses);
				double Eff = normalized_request_rate / sim.measurements.get_av_scheduled_customers();
				vE.push_back(Eff);

				outplot << vB.back() << '\t' << vE.back() << '\n';
			}
		}

		//plt::loglog(vB.toStdVector(), vE.toStdVector());
		//plt::show();

		outplot.close();

		emit ProcessTextChanged("Done simulating");

	}
	else if (number_of_request_rates > 1)
	{
		QVector<double> vx, vC;

		emit ProcessTextChanged("Starting simulation");
		emit GraphChanged(vx, vC);
		//ui.textEdit_6->setPlainText(ui.textEdit_6->toPlainText() + "\nStarting simulation\n");

		std::ofstream outplot;
		std::string plotfilename = filename + "plot.dat";
		if (!CUtility::file_exists(plotfilename))
			outplot.open(plotfilename.c_str());
		else
		{
			std::string errormessage = "File " + plotfilename + " already exists.";
			emit ErrorMessage(errormessage.c_str());
			return;
		}

		double index_low = log(normalized_request_rate_from);
		double index_high = log(normalized_request_rate_to);
		for (double iter = index_low; iter <= (index_high + 0.0001); iter += (index_high - index_low) / (number_of_request_rates - 1))
		{
			if (stop)
				break;
			normalized_request_rate = exp(iter);
			QString newText = "Simulating x = " + QString::fromStdString(std::to_string(normalized_request_rate));
			emit ProcessTextChanged(newText);
			emit GraphChanged(vx, vC);
			//ui.textEdit_6->setPlainText(newText.c_str());

			sim.init_new_sim(number_of_buses, number_of_nodes, topology, normalized_request_rate, bus_type, 1000 * number_of_buses);

			//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
			sim.enable_measurements((1.0 * number_of_buses) / sim.request_rate);
			//simulate (and measure) for 100 requests per bus (at least 10000 requests)
			sim.run_sim_requests(std::max((ULL)30000, 1000 * number_of_buses));

			//output results
			if (save)
			{
				sim.print_params(out, true);
				sim.print_measurements(out);
				out << std::endl << std::endl;
			}

			vx.push_back((double)normalized_request_rate);
			vC.push_back(sim.measurements.get_av_scheduled_customers());

			outplot << vx.back() << '\t' 
				<< vC.back() << '\t' 
				<< sim.measurements.get_stddev_scheduled_customers() << '\t' 
				<< CUtility::two_node_av_scheduled_customers(4, normalized_request_rate, 2) << '\t'
				<< CUtility::two_node_stddev_scheduled_customers(4, normalized_request_rate, 2) << '\n';
		}

		//plt::loglog(vB.toStdVector(), vE.toStdVector());
		//plt::show();

		outplot.close();

		emit ProcessTextChanged("Done simulating");
	}
	else
	{
		emit ProcessTextChanged("Starting simulation: ");

		sim.init_new_sim(number_of_buses, number_of_nodes, topology, normalized_request_rate, bus_type, 10 * number_of_buses);

		//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
		sim.enable_measurements((1.0 * number_of_buses) / sim.request_rate);
		//simulate (and measure) for 100 requests per bus (at least 10000 requests)
		sim.run_sim_requests(std::max((ULL)10000, 100 * number_of_buses));

		emit ProcessTextChanged("Done simulating!");

		//output results
		sim.print_params(out, true);
		sim.print_measurements(out);

		out << std::endl << std::endl;
	}

	out.close();
}