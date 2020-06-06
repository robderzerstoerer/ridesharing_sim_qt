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
	mThread[0] = new simulation_thread;
	mThread[1] = new simulation_thread;
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
	connect(mThread[0], SIGNAL(ProcessTextChanged(QString)), this, SLOT(onProcessTextChanged(QString)));
	connect(mThread[0], SIGNAL(GraphChanged(QVector<double>, QVector<double>)), this, SLOT(onGraphChanged(QVector<double>, QVector<double>)));
	connect(mThread[0], SIGNAL(ErrorMessage(QString)), this, SLOT(onErrorMessage(QString)));
	connect(mThread[1], SIGNAL(ProcessTextChanged(QString)), this, SLOT(onProcessTextChanged(QString)));
	connect(mThread[1], SIGNAL(GraphChanged(QVector<double>, QVector<double>)), this, SLOT(onGraphChanged(QVector<double>, QVector<double>)));
	connect(mThread[1], SIGNAL(ErrorMessage(QString)), this, SLOT(onErrorMessage(QString)));
	
	ui.listWidget->addItem("torus");
	ui.listWidget->addItem("ring");
	ui.listWidget->addItem("directed ring");
	ui.listWidget->addItem("two_nodes");
	ui.listWidget->addItem("double ring");
	ui.listWidget->addItem("ladder");
	ui.listWidget->addItem("3ladder");
	ui.listWidget->addItem("4ladder");
	ui.listWidget->addItem("star");
	ui.listWidget->addItem("grid");
	ui.listWidget->addItem("simplified_city");
	ui.listWidget->addItem("complete_graph");
	ui.listWidget->addItem("3cayley_tree");
	ui.listWidget->addItem("poisson_random");
	ui.listWidget->addItem("delaunay_random_torus");
	ui.listWidget->addItem("gabriel_random_torus");
	ui.listWidget->addItem("berlin_mitte");
	ui.listWidget->addItem("bornholm");
	ui.listWidget->addItem("goettingen");
	ui.listWidget->addItem("isle_of_man");
	ui.listWidget->addItem("jersey");
	ui.listWidget->addItem("korfu");
	ui.listWidget->addItem("mallorca");
	ui.listWidget->addItem("manhattan");
	ui.listWidget->addItem("modena");
	ui.listWidget->addItem("ostfriesland");
	ui.listWidget->addItem("rural_poland");
	
	ui.listWidget_4->addItem("SIMPLE_SIM");
	ui.listWidget_4->addItem("EVERYTHING");
	ui.listWidget_4->addItem("B_AIM");


	ui.listWidget->setItemSelected(ui.listWidget->item(0),true);
	ui.listWidget_4->setItemSelected(ui.listWidget_4->item(0), true);

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
	ui.textEdit_13->setText("50");
	ui.textEdit_14->setText("0.5");

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

	std::string message;

	traffic_network n;
	std::string testnetwork = "complete_graph";
	int test_node_number = 5;
	n.init(testnetwork.c_str(), test_node_number);

	std::ofstream edgefile;
	edgefile.open(testnetwork + "_N_" + std::to_string(test_node_number) + "_edgelist.dat");

	std::map< ULL, std::set< std::pair<ULL, double> > > edgelist = n.get_edgelist();
	std::vector<std::tuple<ULL, ULL, double>> edgevector;
	
	for (std::map< ULL, std::set< std::pair<ULL, double> > >::iterator it = edgelist.begin(); it != edgelist.end(); it++)
	{
		for (std::set< std::pair<ULL, double> >::iterator itt = it->second.begin(); itt != it->second.end(); itt++)
		{
			if (itt->first > it->first)
				edgevector.push_back(std::make_tuple(it->first, itt->first, itt->second));
		}
	}

	for (int i = 0; i < edgevector.size(); i++)
	{
		ULL from;
		ULL to;
		double dist;
		std::tie(from, to, dist) = edgevector[i];
		edgefile << "$" <<from + 1 << "$" << " & " << "$" << to + 1 << "$" << " & " << "$" << dist << "$" << "\\\\ \\hline \n";
	}

	edgefile.close();

	message += "p_l1 = " + std::to_string(n.get_p_l1()) + '\n' + 
		"p_2n = " + std::to_string(n.get_p_2n()) + '\n'+
		"dav = " + std::to_string(n.get_dav()) + '\n';

	message += std::to_string(120) + ": " + std::to_string(CUtility::calc_p_full(44, 140, 8)) + '\n';
	/*
	for (int a_c = 0; a_c <= 15; a_c++)
	{
		message += std::to_string(a_c*100) + ": " + std::to_string(CUtility::calc_p_full(300, a_c*100, 5)) + '\n';
	} */
	
	QMessageBox Msgbox(QMessageBox::Icon::Information, "Output", message.c_str());
	Msgbox.exec();
}

void ridesharing_sim_qt::onButtonSimulateClicked()
{
	// make sure we do not use old values of B, x, capacities and the topologies
	mThread[0]->par.B_list.clear();
	mThread[0]->par.x_list.clear();
	mThread[0]->par.topology_list.clear();
	mThread[0]->par.capacity_list.clear();

	ui.customPlot->graph(0)->setData(QVector<double>(), QVector<double>());
	ui.customPlot->graph(0)->setPen(QPen(Qt::blue));
	ui.customPlot->replot();

	// global simulation variables
	mThread[0]->par.simulate_until_exact = true;
	mThread[0]->par.foldername = "allsim_pfull_correct\\";
	mThread[0]->par.num_requests_per_bus_init = 2000;
	mThread[0]->par.num_requests_per_bus_sim = 2000;
	mThread[0]->par.calc_cap_delay = false;

	if (ui.listWidget_4->selectedItems().front()->text().toStdString() == "SIMPLE_SIM")
		mThread[0]->par.sim_mode = SIMPLE_SIM;
	else if (ui.listWidget_4->selectedItems().front()->text().toStdString() == "EVERYTHING")
		mThread[0]->par.sim_mode = EVERYTHING;
	else if (ui.listWidget_4->selectedItems().front()->text().toStdString() == "B_AIM")
		mThread[0]->par.sim_mode = B_AIM;

	if (mThread[0]->par.sim_mode == EVERYTHING)
	{
		mThread[0]->par.capacity_list = {1, 2, 3, 4, 5, 6, 7, 8};

		mThread[0]->par.scan_point_sim_time = 60; // 1min
		mThread[0]->par.max_point_sim_time = 6*60; // 6min
	}
	else if (mThread[0]->par.sim_mode == SIMPLE_SIM)
	{
		mThread[0]->par.capacity_list = { std::stoi(ui.textEdit_12->toPlainText().toStdString()) };

		// use topology input from ListWidget and TextEdit
		mThread[0]->par.topology_list = { std::make_pair(ui.listWidget->selectedItems().front()->text().toStdString(), std::stoi(ui.textEdit->toPlainText().toStdString())) };

		mThread[0]->par.scan_point_sim_time = 60 * 60; // 1h
		mThread[0]->par.max_point_sim_time = 4 * 60 * 60; // 4h
	}
	else if (mThread[0]->par.sim_mode == B_AIM)
	{
		int cap_from_str = std::stoi(ui.textEdit_12->toPlainText().toStdString());
		if (cap_from_str != -1)
		{
			int max_cap = 12;
			int min_cap = 1;
			for (int i = max_cap; i >= min_cap; i--)
				mThread[0]->par.capacity_list.push_back(i);
		}
		else
		{
			mThread[0]->par.capacity_list = { -1 };
		}

		int N = std::stoi(ui.textEdit->toPlainText().toStdString());

		mThread[0]->par.simulate_until_exact = false;
		mThread[0]->par.scan_point_sim_time = 1e10;
		mThread[0]->par.max_point_sim_time = 1e10;

		mThread[0]->par.lambda = std::stod(ui.textEdit_13->toPlainText().toStdString());
		mThread[0]->par.E_aim = std::stod(ui.textEdit_14->toPlainText().toStdString());

		//traffic_network net();
		//net.init(topo, N);

	}
	

	if (ui.checkBox->isChecked())
	{
		int B_from = std::stoi(ui.textEdit_2->toPlainText().toStdString());
		int B_to = std::stoi(ui.textEdit_3->toPlainText().toStdString());
		int num_B = std::stoi(ui.textEdit_7->toPlainText().toStdString());
		
		// logarithmically spaced Bs
		double index_low_B = log(B_from);
		double index_high_B = log(B_to);
		for (double iter_B = index_low_B; iter_B <= (index_high_B + 0.0001); iter_B += (index_high_B - index_low_B) / (num_B - 1))
		{
			int B = (int)(exp(iter_B) + 0.5);
			mThread[0]->par.B_list.push_back(B);
		}
	}
	else
	{
		// just use this single value of B
		mThread[0]->par.B_list.push_back(std::stoi(ui.textEdit_2->toPlainText().toStdString()));
	}	
	

	if (ui.checkBox_3->isChecked())
	{
		int num_x = std::stoi(ui.textEdit_10->toPlainText().toStdString());
		double x_from = std::stod(ui.textEdit_4->toPlainText().toStdString());
		double x_to = std::stod(ui.textEdit_11->toPlainText().toStdString());

		// use linearly spaced values for x
		for (double x = x_from; x <= (x_to + 0.0001); x += (x_to - x_from) / (num_x -1))
		{
			mThread[0]->par.x_list.push_back(x);
		}
	}
	else
	{
		// just use this single value of x
		mThread[0]->par.x_list.push_back(std::stod(ui.textEdit_4->toPlainText().toStdString()));
	}

	mThread[0]->par.calc_p_full = ui.checkBox_4->isChecked();
	
	mThread[0]->par.filename = ui.textEdit_5->toPlainText().toStdString();

	mThread[1]->par = mThread[0]->par;

	if (mThread[0]->par.sim_mode == EVERYTHING)
	{
		mThread[0]->par.topology_list = {
			std::make_pair("two_nodes", 2),
			std::make_pair("simplified_city", 16),
			std::make_pair("ring", 25),
			std::make_pair("ring", 100),
			std::make_pair("torus", 25),
			std::make_pair("torus", 100),
			// std::make_pair("grid", 100)
			// std::make_pair("complete_graph", 5),
			// std::make_pair("3cayley_tree", 46),
			// std::make_pair("delaunay_random_torus", 25),
			// std::make_pair("delaunay_random_torus", 100),
			// std::make_pair("bornholm", 100),
		};
		mThread[1]->par.topology_list = {
			// std::make_pair("two_nodes", 2),
			// std::make_pair("simplified_city", 16),
			// std::make_pair("ring", 25),
			// std::make_pair("ring", 100),
			// std::make_pair("torus", 25),
			// std::make_pair("torus", 100),
			std::make_pair("grid", 100),
			std::make_pair("complete_graph", 5),
			std::make_pair("3cayley_tree", 46),
			std::make_pair("delaunay_random_torus", 25),
			std::make_pair("delaunay_random_torus", 100)
			// std::make_pair("bornholm", 100),
		};
	}
	else if (mThread[0]->par.sim_mode == B_AIM)
	{
		mThread[0]->par.topology_list =
		{
			std::make_pair("torus", 100),
			std::make_pair("grid", 100),
		};

		mThread[1]->par.topology_list =
		{
			std::make_pair("3cayley_tree", 46),
			std::make_pair("delaunay_random_torus", 25),
			std::make_pair("delaunay_random_torus", 100)
		};
	}
	else if (mThread[0]->par.sim_mode == SIMPLE_SIM)
	{
		/*
		mThread[0]->par.topology_list = {
			// std::make_pair("two_nodes", 2),
			std::make_pair("simplified_city", 16),
			// std::make_pair("ring", 25),
			// std::make_pair("ring", 100),
			std::make_pair("torus", 25),
			std::make_pair("torus", 100),
			std::make_pair("grid", 100)
			// std::make_pair("complete_graph", 5),
			// std::make_pair("3cayley_tree", 46),
			// std::make_pair("delaunay_random_torus", 25),
			// std::make_pair("delaunay_random_torus", 100),
			// std::make_pair("bornholm", 100),
		};
		mThread[1]->par.topology_list = {
			// std::make_pair("two_nodes", 2),
			// std::make_pair("simplified_city", 16),
			// std::make_pair("ring", 25),
			// std::make_pair("ring", 100),
			// std::make_pair("torus", 25),
			// std::make_pair("torus", 100),
			// std::make_pair("grid", 100),
			// std::make_pair("complete_graph", 5),
			std::make_pair("3cayley_tree", 46),
			std::make_pair("delaunay_random_torus", 25),
			std::make_pair("delaunay_random_torus", 100)
			// std::make_pair("bornholm", 100),
		};
		*/
	}

	// Start Thread
	mThread[0]->start();
	// multithreading only if parallel computing makes sense
	if (mThread[0]->par.sim_mode == EVERYTHING || mThread[0]->par.sim_mode == B_AIM || mThread[0]->par.sim_mode == SIMPLE_SIM)
		mThread[1]->start();
	
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
	mThread[0]->par.stop_thread = true;
	mThread[1]->par.stop_thread = true;
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
	if (CUtility::read_file(filename, averages[iFile], stddevs[iFile], counts[iFile]))
	{
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