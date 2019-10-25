#pragma once

#include <QThread>
#include "ridesharing_sim.h"

class simulation_thread : public QThread
{
	Q_OBJECT

public:
	explicit simulation_thread(QObject *parent = 0);
	~simulation_thread();
	void run();

	void simulate_B_list(simulation_parameters& sim_par,
		std::ofstream& out,
		std::ofstream& outplot,
		std::vector<std::pair<ULL, double>>& vBEdata = std::vector<std::pair<ULL, double>>());

	// returns Efficiency and  if job cancelled or failed
	double single_simulation(simulation_parameters &sim_par, 
		std::ofstream& out, 
		std::ofstream& outplot, 
		std::vector<std::pair<ULL, double>>& vBEdata = std::vector<std::pair<ULL, double>>());

	program_parameters par;

signals:
	void ProcessTextChanged(QString);
	void GraphChanged(QVector<double>, QVector<double>);
	void ErrorMessage(QString);
};
