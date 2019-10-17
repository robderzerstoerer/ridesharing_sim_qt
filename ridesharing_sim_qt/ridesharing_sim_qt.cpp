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
#include <qfiledialog.h>
#include <qdir.h>

//#include "matplotlib.h"
//namespace plt = matplotlibcpp;

#include "ridesharing_sim.h"
#include "utility.h"


#ifndef _INTEGER_TYPES
#define ULL uint64_t
#define LL int64_t
#define _INTEGER_TYPES
#endif

#ifndef _EPSILON
#define MACRO_EPSILON 0.000001
#define _EPSILON
#endif


#include "ridesharing_sim_qt.h"
#include <qmessagebox.h>

ridesharing_sim_qt::ridesharing_sim_qt(QWidget *parent)
	: QMainWindow(parent)
{
	mThread = new simulation_thread;
	ui.setupUi(this);
	connect(ui.pushButton, SIGNAL(clicked()), this, SLOT(onButtonSimulateClicked()));
	connect(ui.pushButton_2, SIGNAL(clicked()), this, SLOT(onButtonStopClicked()));
	connect(ui.checkBox, SIGNAL(stateChanged(int)), this, SLOT(B_to_checkBoxChanged(int)));
	connect(ui.pushButton_3, SIGNAL(clicked()), this, SLOT(onButtonChooseFile1Clicked()));
	connect(ui.pushButton_4, SIGNAL(clicked()), this, SLOT(onButtonChooseFile2Clicked()));
	connect(ui.pushButton_5, SIGNAL(clicked()), this, SLOT(onButtonSwapXYClicked()));
	connect(ui.pushButton_6, SIGNAL(clicked()), this, SLOT(onButtonCompareClicked()));
	connect(ui.pushButton_7, SIGNAL(clicked()), this, SLOT(onButtonSavePlotClicked()));

	qRegisterMetaType<QVector<double>>("QVector<double>");
	connect(mThread, SIGNAL(ProcessTextChanged(QString)), this, SLOT(onProcessTextChanged(QString)));
	connect(mThread, SIGNAL(GraphChanged(QVector<double>, QVector<double>)), this, SLOT(onGraphChanged(QVector<double>, QVector<double>)));
	connect(mThread, SIGNAL(ErrorMessage(QString)), this, SLOT(onErrorMessage(QString)));
	
	ui.listWidget->addItem("torus");
	ui.listWidget->addItem("ring");
	ui.listWidget->addItem("directed ring");
	ui.listWidget->addItem("two_nodes");

	ui.listWidget->setItemSelected(ui.listWidget->item(0),true);

	ui.checkBox->setChecked(true);
	ui.checkBox_3->setChecked(false);
	ui.checkBox_4->setChecked(true);

	ui.textEdit->setText("100");
	ui.textEdit_2->setText("1");
	ui.textEdit_3->setText("3000");
	ui.textEdit_4->setText("2.0");
	ui.textEdit_5->setText("test.dat");
	ui.textEdit_6->setDisabled(true);
	ui.textEdit_7->setText("12");
	ui.textEdit_8->setDisabled(true);
	ui.textEdit_9->setDisabled(true);
	ui.textEdit_12->setText("5");

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
	ui.customPlot->addGraph();

	
	//QMessageBox Msgbox(QMessageBox::Icon::Information, "Output", message.c_str());
	//Msgbox.exec();
}

void ridesharing_sim_qt::onButtonSimulateClicked()
{
	ui.customPlot->graph(0)->setData(QVector<double>(), QVector<double>());
	ui.customPlot->graph(0)->setPen(QPen(Qt::blue));
	ui.customPlot->replot();
	
	// Choose topology
	mThread->par.topology = ui.listWidget->selectedItems().front()->text().toStdString();

	if (ui.checkBox->isChecked())
	{
		mThread->par.number_of_buses = 0;
		mThread->par.number_of_buses_from = std::stoi(ui.textEdit_2->toPlainText().toStdString());
		mThread->par.number_of_buses_to = std::stoi(ui.textEdit_3->toPlainText().toStdString());
		mThread->par.number_of_bus_calculations = std::stoi(ui.textEdit_7->toPlainText().toStdString());
	}
	else
	{
		mThread->par.number_of_buses = std::stoi(ui.textEdit_2->toPlainText().toStdString());
		mThread->par.number_of_buses_from = 0;
		mThread->par.number_of_buses_to = 0;
		mThread->par.number_of_bus_calculations = 1;
	}
	mThread->par.number_of_nodes = std::stoi(ui.textEdit->toPlainText().toStdString());
	

	if (ui.checkBox_3->isChecked())
	{
		mThread->par.number_of_request_rates = std::stoi(ui.textEdit_10->toPlainText().toStdString());
		mThread->par.normalized_request_rate = std::stod(ui.textEdit_4->toPlainText().toStdString());
		mThread->par.normalized_request_rate_from = std::stod(ui.textEdit_4->toPlainText().toStdString());
		mThread->par.normalized_request_rate_to = std::stod(ui.textEdit_11->toPlainText().toStdString());
	}
	else
	{
		mThread->par.number_of_request_rates = 1;
		mThread->par.normalized_request_rate = std::stod(ui.textEdit_4->toPlainText().toStdString());
		mThread->par.normalized_request_rate_from = mThread->par.normalized_request_rate;
		mThread->par.normalized_request_rate_to = mThread->par.normalized_request_rate;
	}

	int cap = std::stoi(ui.textEdit_12->toPlainText().toStdString());
	if (cap == -1)
	{
		mThread->par.bus_type = 0;  // unlimited capacity
	}
	else
	{
		mThread->par.bus_type = 1;  // limited capacity
	}

	mThread->par.simulate_until_exact = true;

	mThread->par.calc_p_full = ui.checkBox_4->isChecked();
	mThread->par.calc_cap_delay = false;
	
	mThread->par.save = true;

	if (mThread->par.save)
	{
		mThread->par.filename = ui.textEdit_5->toPlainText().toStdString();
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


void ridesharing_sim_qt::onButtonStopClicked()
{
	mThread->par.stop_thread = true;
}

void ridesharing_sim_qt::onProcessTextChanged(QString newText)
{
	ui.textEdit_6->setPlainText(newText);
}

void ridesharing_sim_qt::onErrorMessage(QString errorText)
{
	QMessageBox Msgbox(QMessageBox::Icon::Critical, "Error!", errorText);
	Msgbox.exec();
}

void ridesharing_sim_qt::onGraphChanged(QVector<double> vB, QVector<double> vE)
{
	ui.customPlot->graph(0)->setData(vB, vE);
	if (!vB.isEmpty() && !vE.isEmpty())
	{
		double minx = 1e10;
		double maxx = -1e10;
		double miny = 1e10;
		double maxy = -1e10;
		for (int i = 0; i < vB.size(); i++)
		{
			if (vB[i] < minx)
				minx = vB[i];
			if (vB[i] > maxx)
				maxx = vB[i];
		}
		for (int i = 0; i < vE.size(); i++)
		{
			if (vE[i] < miny)
				miny = vE[i];
			if (vE[i] > maxy)
				maxy = vE[i];
		}
		ui.customPlot->xAxis->setRange(minx - 0.1 * abs(minx), maxx + 0.1 * abs(maxx));
		ui.customPlot->yAxis->setRange(miny - 0.1 * abs(miny), maxy + 0.1 * abs(maxy));
	}
	ui.customPlot->replot();
}

void ridesharing_sim_qt::read_file(std::string filename, int iFile)
{
	for (int i = 0; i < 3; i++)
	{
		averages[iFile].clear();
		stddevs[iFile].clear();
		counts[iFile].clear();
	}

	std::string line;
	std::ifstream myfile(filename);
	std::string prevLine = "";
	bool NumberNext = false;
	if (myfile.is_open())
	{
		while (std::getline(myfile, line))
		{
			if (line == "PARAMS" || line == "MEASUREMENTS" || line == "")
			{
				NumberNext = false;
			}
			else if (NumberNext == true)
			{
				// extract n, av, stddev
				int find1 = line.find('\t');
				if (find1 != std::string::npos)
				{
					int find2 = line.find('\t', find1 + 1);
					counts[iFile][prevLine].push_back(std::stod(line.substr(0, find1)));
					averages[iFile][prevLine].push_back(std::stod(line.substr(find1 + 1, find2 - find1 - 1)));
					stddevs[iFile][prevLine].push_back(std::stod(line.substr(find2 + 1, line.length() - find2 - 1)));
				}
				else	// no standard deviation or multiple measurements stored
				{
					averages[iFile][prevLine].push_back(std::stod(line));
					stddevs[iFile][prevLine].push_back(-1.0);
					counts[iFile][prevLine].push_back(-1.0);
				}

				NumberNext = false;
			}
			else // new variable incoming
			{
				prevLine = line;
				NumberNext = true;
			}
		}
		myfile.close();

		QListWidget* list = NULL;
		if (iFile == 0)
			list = ui.listWidget_2;
		else if (iFile == 1)
			list = ui.listWidget_3;

		list->clear();
		for (std::map<std::string, std::vector<double>>::iterator it = averages[iFile].begin(); it != averages[iFile].end(); it++)
		{
			list->addItem(QString::fromStdString(it->first));
		}
	}
	else
	{
		QMessageBox Msgbox(QMessageBox::Icon::Critical, "Error!", "Could not open File.");
		Msgbox.exec();
	}
}


void ridesharing_sim_qt::onButtonChooseFile1Clicked()
{
	QString file1 = QFileDialog::getOpenFileName(this, "Open File 1 to compare", QDir::homePath());
	if (file1 != "")
		ui.textEdit_8->setPlainText(file1);

	read_file(file1.toStdString(), 0);
}

void ridesharing_sim_qt::onButtonChooseFile2Clicked()
{
	QString file2 = QFileDialog::getOpenFileName(this, "Open File 2 to compare", QDir::homePath());
	if (file2 != "")
		ui.textEdit_9->setPlainText(file2);

	read_file(file2.toStdString(), 1);
}

void ridesharing_sim_qt::onButtonSwapXYClicked()
{
	QString left = ui.label_8->text();
	ui.label_8->setText(ui.label_9->text());
	ui.label_9->setText(left);
}

void ridesharing_sim_qt::onButtonCompareClicked()
{
	if (ui.listWidget_2->selectedItems().size() != 1 || ui.listWidget_3->selectedItems().size() != 1)
	{
		QMessageBox Msgbox(QMessageBox::Icon::Critical, "Error!", "Please select items to be compared.");
		Msgbox.exec();
		return;
	}
	QVector<double> vec1 = QVector<double>::fromStdVector(averages[0][ui.listWidget_2->selectedItems().front()->text().toStdString()]);
	QVector<double> vec2 = QVector<double>::fromStdVector(averages[1][ui.listWidget_3->selectedItems().front()->text().toStdString()]);
	
	if (ui.label_8->text() == "x")
		onGraphChanged(vec1, vec2);
	else
		onGraphChanged(vec2, vec1);
}

void ridesharing_sim_qt::onButtonSavePlotClicked()
{
	// std::ofstream out(mThread->filename + "plot.dat");
}