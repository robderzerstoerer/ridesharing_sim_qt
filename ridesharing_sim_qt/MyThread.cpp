#include "simulation_thread.h"

simulation_thread::simulation_thread(QObject *parent)
	: QThread(parent)
{
}

simulation_thread::~simulation_thread()
{
}


void simulation_thread::run()
{
	stop = false;
	std::ofstream out(filename.c_str());
	//initialize simulation class
	ridesharing_sim sim(number_of_nodes, number_of_buses, 0);
	sim.init_network(number_of_buses, number_of_nodes, topology);

	if (number_of_bus_calculations > 1)
	{
		QVector<double> vB, vE;

		emit ProcessTextChanged("Starting simulation");
		emit GraphChanged(vB, vE);
		//ui.textEdit_6->setPlainText(ui.textEdit_6->toPlainText() + "\nStarting simulation\n");
		
		double index_low = log(number_of_buses_from);
		double index_high = log(number_of_buses_to);
		for (double iter = std::max(1, index_low); iter <= (index_high + 0.0001); iter += (index_high-index_low)/(number_of_bus_calculations-1))
		{
			if (stop)
				break;
			number_of_buses = (int)(10 * exp((double)iter / 5));
			QString newText = "Simulating " + QString::fromStdString(std::to_string(number_of_buses)) + " buses\n";
			emit ProcessTextChanged(newText);
			emit GraphChanged(vB, vE);
			//ui.textEdit_6->setPlainText(newText.c_str());
			
			//Returning focus to the edit
			//ui.textEdit_6->setFocus();
			//!!! Here I want to move the cursor 4 characters left to place it before the </i> tag.
			//ui.textEdit_6->textCursor().movePosition(QTextCursor::EndOfBlock, QTextCursor::MoveAnchor);

			sim.init_new_sim(number_of_buses, number_of_nodes, topology, normalized_request_rate, 10 * number_of_buses);

			//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
			sim.enable_measurements((1.0 * number_of_buses) / sim.request_rate);
			//simulate (and measure) for 100 requests per bus (at least 10000 requests)
			sim.run_sim_requests(std::max((ULL)10000, 100 * number_of_buses));

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
		}

		//plt::loglog(vB.toStdVector(), vE.toStdVector());
		//plt::show();

		emit ProcessTextChanged("Done simulating");

	}
	else
	{
		emit ProcessTextChanged("Starting simulation: ");

		sim.init_new_sim(number_of_buses, number_of_nodes, topology, normalized_request_rate, 10 * number_of_buses);

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