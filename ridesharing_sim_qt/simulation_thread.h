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

	program_parameters par;

signals:
	void ProcessTextChanged(QString);
	void GraphChanged(QVector<double>, QVector<double>);
	void ErrorMessage(QString);
};
