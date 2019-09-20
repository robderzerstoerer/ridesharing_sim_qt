#include <cstddef>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <queue>
#include <map>
#include <list>
#include <set>
#include <tuple>
#include <vector>
#include <cmath>
#include <random>
#include <cassert>
#include <qmessagebox.h>

//#include "matplotlib.h"
//namespace plt = matplotlibcpp;

#include "ridesharing_sim.h"


#ifndef _INTEGER_TYPES
#define ULL uint64_t
#define LL int64_t
#define _INTEGER_TYPES
#endif

#ifndef _EPSILON
#define MACRO_EPSILON 0.000001
#define _EPSILON
#endif

constexpr double pi = 3.14159265358979323846;

#include "ridesharing_sim_qt.h"
#include <qmessagebox.h>

ridesharing_sim_qt::ridesharing_sim_qt(QWidget *parent)
	: QMainWindow(parent)
{
	mThread = new simulation_thread;
	ui.setupUi(this);
	connect(ui.pushButton, SIGNAL(clicked()), this, SLOT(ifkl()));
	connect(ui.pushButton_2, SIGNAL(clicked()), this, SLOT(ifklstopped()));
	connect(ui.checkBox, SIGNAL(stateChanged(int)), this, SLOT(B_to_checkBoxChanged(int)));
	connect(ui.checkBox_2, SIGNAL(stateChanged(int)), this, SLOT(SaveCheckBoxChanged(int)));

	qRegisterMetaType<QVector<double>>("QVector<double>");
	connect(mThread, SIGNAL(ProcessTextChanged(QString)), this, SLOT(onProcessTextChanged(QString)));
	connect(mThread, SIGNAL(GraphChanged(QVector<double>, QVector<double>)), this, SLOT(onGraphChanged(QVector<double>, QVector<double>)));
	
	ui.listWidget->addItem("torus");
	ui.listWidget->addItem("ring");
	ui.listWidget->addItem("two_nodes");

	ui.listWidget->setItemSelected(ui.listWidget->item(0),true);

	ui.checkBox->setChecked(true);
	ui.checkBox_2->setChecked(true);

	ui.textEdit->setText("25");
	ui.textEdit_2->setText("30");
	ui.textEdit_3->setText("1500");
	ui.textEdit_4->setText("7.5");
	ui.textEdit_5->setText("test.dat");
	ui.textEdit_6->setDisabled(true);
	ui.textEdit_7->setText("12");

	ui.customPlot->xAxis->setLabel("x");
	ui.customPlot->yAxis->setLabel("y");
	// set axes ranges, so we see all data:
	ui.customPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
	ui.customPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
	QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
	ui.customPlot->xAxis->setTicker(logTicker);
	ui.customPlot->yAxis->setTicker(logTicker);
	ui.customPlot->xAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
	ui.customPlot->xAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
	ui.customPlot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
	ui.customPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
}

void ridesharing_sim_qt::ifkl()
{
	ui.customPlot->addGraph();
	ui.customPlot->graph(0)->setData(QVector<double>(), QVector<double>());
	ui.customPlot->graph(0)->setPen(QPen(Qt::blue));
	ui.customPlot->replot();
	
	// Choose topology
	mThread->topology = ui.listWidget->selectedItems().front()->text().toStdString();

	if (ui.checkBox->isChecked())
	{
		mThread->number_of_buses = 0;
		mThread->number_of_buses_from = std::stoi(ui.textEdit_2->toPlainText().toStdString());
		mThread->number_of_buses_to = std::stoi(ui.textEdit_3->toPlainText().toStdString());
		mThread->number_of_bus_calculations = std::stoi(ui.textEdit_7->toPlainText().toStdString());
	}
	else
	{
		mThread->number_of_buses = std::stoi(ui.textEdit_2->toPlainText().toStdString());
		mThread->number_of_buses_from = 0;
		mThread->number_of_buses_to = 0;
		mThread->number_of_bus_calculations = 1;
	}
	mThread->number_of_nodes = std::stoi(ui.textEdit->toPlainText().toStdString());
	mThread->normalized_request_rate = std::stod(ui.textEdit_4->toPlainText().toStdString());
	
	mThread->save = ui.checkBox_2->isChecked();
	if (mThread->save)
	{
		mThread->filename = ui.textEdit_5->toPlainText().toStdString();
	}

	// Start Thread
	mThread->start();

	
}

void ridesharing_sim_qt::B_to_checkBoxChanged(int arg1)
{
	if (ui.checkBox->isChecked())
	{
		ui.textEdit_3->setDisabled(false);
		ui.textEdit_7->setDisabled(false);
	}
	else
	{
		ui.textEdit_3->setDisabled(true);
		ui.textEdit_7->setDisabled(true);
	}
}

void ridesharing_sim_qt::SaveCheckBoxChanged(int arg1)
{
	if (ui.checkBox_2->isChecked())
		ui.textEdit_5->setDisabled(false);
	else
		ui.textEdit_5->setDisabled(true);
}

void ridesharing_sim_qt::ifklstopped()
{
	mThread->stop = true;
}

void ridesharing_sim_qt::onProcessTextChanged(QString newText)
{
	ui.textEdit_6->setPlainText(newText);
}

void ridesharing_sim_qt::onGraphChanged(QVector<double> vB, QVector<double> vE)
{
	ui.customPlot->graph(0)->setData(vB, vE);
	if (!vB.isEmpty() && !vE.isEmpty())
	{
		ui.customPlot->xAxis->setRange(vB.first() / 2, vB.last() + 0.1);
		ui.customPlot->yAxis->setRange(vE.first() / 2, vE.last() + 0.1);
	}
	ui.customPlot->replot();
}