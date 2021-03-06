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
	par.stop = false;
	std::ofstream out;
	if (!CUtility::file_exists(par.filename))
		out.open(par.filename.c_str());
	else
	{
		std::string errormessage = "File " + par.filename + " already exists.";
		emit ErrorMessage(errormessage.c_str());
		return;
	}

	//initialize simulation class
	ridesharing_sim sim(&par);
	sim.init_network();

	if (par.number_of_bus_calculations > 1)
	{
		QVector<double> vB, vE;
		QVector<double> vx, vC, vB_equiv_x_inf, vB_equiv_x;
		QVector<QPair<ULL, double>> vBEdata;

		if (par.number_of_request_rates > 1)
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
		std::string plotfilename = par.filename + "eff.dat";
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

		double index_low_B = log(par.number_of_buses_from);
		double index_high_B = log(par.number_of_buses_to);
		for (double iter_B = index_low_B; iter_B <= (index_high_B + 0.0001); iter_B += (index_high_B - index_low_B) / (par.number_of_bus_calculations - 1))
		{
			if (par.stop)
				break;
			ULL number_of_buses = (int)(exp(iter_B) + 0.5);
			QString newText = "Simulating " + QString::fromStdString(std::to_string(number_of_buses)) + " buses\n";
			emit ProcessTextChanged(newText);
			emit GraphChanged(vB, vE);
			
			if (par.number_of_request_rates > 1)
			{
				double index_low_x = log(par.normalized_request_rate_from);
				double index_high_x = log(par.normalized_request_rate_to);
				for (double iter_x = index_low_x; iter_x <= (index_high_x + 0.0001); iter_x += (index_high_x - index_low_x) / (par.number_of_request_rates - 1))
				{
					if (par.stop)
						break;
					double normalized_request_rate = exp(iter_x);
					QString newText = "Simulating B = " + QString::fromStdString(std::to_string(number_of_buses)) + ", x = " + QString::fromStdString(std::to_string(normalized_request_rate));
					emit ProcessTextChanged(newText);
					emit GraphChanged(vx, vC);
					//ui.textEdit_6->setPlainText(newText.c_str());

					sim.init_new_sim(number_of_buses, normalized_request_rate, 500 * number_of_buses);

					//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
					sim.enable_measurements((1.0 * number_of_buses) / sim.request_rate);
					//simulate (and measure) for 100 requests per bus (at least 10000 requests)
					sim.run_sim_requests(std::max((ULL)10000, 500 * number_of_buses));

					//output results
					if (par.save)
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
				sim.init_new_sim(number_of_buses, par.normalized_request_rate, std::max((ULL)50000, 100 * number_of_buses));
				
				//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
				sim.enable_measurements((2.0 * number_of_buses) / sim.request_rate);
				
				if (par.simulate_until_exact)
				{
					bool exact_result = false;
					while (!exact_result)
					{
						if (par.stop)
							break;

						//simulate (and measure) for 100 requests per bus (at least 10000 requests)
						sim.run_sim_requests(std::max((ULL)50000, 500 * number_of_buses));

						double total_C_av = sim.measurements.get_av_scheduled_customers();
						if (abs(sim.measurements.C_averages.back() - total_C_av) < 0.1 * total_C_av)
							exact_result = true;
						else
						{
							QString newText = QString::fromStdString(std::to_string(total_C_av)) + ", " + QString::fromStdString(std::to_string(sim.measurements.C_averages.back()));
							emit ProcessTextChanged(newText);
						}
					}
				}
				else
				{
					//simulate (and measure) for 100 requests per bus (at least 10000 requests)
					sim.run_sim_requests(std::max((ULL)10000, 500 * number_of_buses));
				}

				//output results
				if (par.save)
				{
					sim.print_params(out, true);
					sim.print_measurements(out);
					out << std::endl << std::endl;
				}

				vB.push_back((double)number_of_buses);
				double total_C_av = sim.measurements.get_av_scheduled_customers();
				double Eff = par.normalized_request_rate / total_C_av;
				vE.push_back(Eff);

				outplot << vB.back() << '\t' << vE.back();
				if (abs(sim.measurements.C_averages.back() - total_C_av) > 0.01 * total_C_av)
				{
					outplot << '\t' << "x" << '\t' << par.normalized_request_rate / sim.measurements.C_averages.back();
				}
				outplot << '\n';
			}
		}

		//plt::loglog(vB.toStdVector(), vE.toStdVector());
		//plt::show();

		outplot.close();

		emit ProcessTextChanged("Done simulating");

	}
	else if (par.number_of_request_rates > 1)
	{
		QVector<double> vx, vC;

		emit ProcessTextChanged("Starting simulation");
		emit GraphChanged(vx, vC);
		//ui.textEdit_6->setPlainText(ui.textEdit_6->toPlainText() + "\nStarting simulation\n");

		std::ofstream outplot;
		std::string plotfilename = par.filename + "plot.dat";
		if (!CUtility::file_exists(plotfilename))
			outplot.open(plotfilename.c_str());
		else
		{
			std::string errormessage = "File " + plotfilename + " already exists.";
			emit ErrorMessage(errormessage.c_str());
			return;
		}

		double index_low = log(par.normalized_request_rate_from);
		double index_high = log(par.normalized_request_rate_to);
		for (double iter = index_low; iter <= (index_high + 0.0001); iter += (index_high - index_low) / (par.number_of_request_rates - 1))
		{
			if (par.stop)
				break;
			double normalized_request_rate = exp(iter);
			QString newText = "Simulating x = " + QString::fromStdString(std::to_string(normalized_request_rate));
			emit ProcessTextChanged(newText);
			emit GraphChanged(vx, vC);
			//ui.textEdit_6->setPlainText(newText.c_str());

			sim.init_new_sim(par.number_of_buses, normalized_request_rate, std::max((ULL)30000, 1000 * par.number_of_buses));

			//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
			sim.enable_measurements((1.0 * par.number_of_buses) / sim.request_rate);

			if (par.simulate_until_exact)
			{
				bool exact_result = false;
				while (!exact_result)
				{
					//simulate (and measure) for 100 requests per bus (at least 10000 requests)
					sim.run_sim_requests(std::max((ULL)50000, 1000 * par.number_of_buses));

					double total_C_av = sim.measurements.get_av_scheduled_customers();
					if (abs(sim.measurements.C_averages.back() - total_C_av) < 0.1 * total_C_av)
						exact_result = true;
				}
			}
			else
			{
				//simulate (and measure) for 100 requests per bus (at least 10000 requests)
				sim.run_sim_requests(std::max((ULL)50000, 1000 * par.number_of_buses));
			}

			//output results
			if (par.save)
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
				<< CUtility::two_node_av_scheduled_customers(5, normalized_request_rate, 2) << '\t'
				<< CUtility::two_node_stddev_scheduled_customers(5, normalized_request_rate, 2);

			if (abs(sim.measurements.C_averages.back() - vC.back()) > 0.01 * vC.back())
			{
				outplot << '\t' << "x" << '\t' << sim.measurements.C_averages.back();
			}

			outplot << '\n';
		}

		//plt::loglog(vB.toStdVector(), vE.toStdVector());
		//plt::show();

		outplot.close();

		emit ProcessTextChanged("Done simulating");
	}
	else
	{
		emit ProcessTextChanged("Starting simulation: ");

		sim.init_new_sim(par.number_of_buses, par.normalized_request_rate, 10 * par.number_of_buses);

		//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
		sim.enable_measurements((1.0 * par.number_of_buses) / sim.request_rate);
		//simulate (and measure) for 100 requests per bus (at least 10000 requests)
		sim.run_sim_requests(std::max((ULL)10000, 100 * par.number_of_buses));

		emit ProcessTextChanged("Done simulating!");

		//output results
		sim.print_params(out, true);
		sim.print_measurements(out);

		out << std::endl << std::endl;
	}

	out.close();
}