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

	std::string topology;
	ULL number_of_buses = 0;
	ULL number_of_buses_from;
	ULL number_of_buses_to;
	ULL number_of_bus_calculations;
	ULL number_of_nodes;
	double normalized_request_rate;
	ULL number_of_request_rates;
	double normalized_request_rate_from;
	double normalized_request_rate_to;
	ULL bus_type;
	bool save;
	std::string filename;

	bool stop;

signals:
	void ProcessTextChanged(QString);
	void GraphChanged(QVector<double>, QVector<double>);
};
