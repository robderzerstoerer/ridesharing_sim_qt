#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_ridesharing_sim_qt.h"
#include "simulation_thread.h"


class ridesharing_sim_qt : public QMainWindow
{
	Q_OBJECT

public:
	ridesharing_sim_qt(QWidget *parent = Q_NULLPTR);
	simulation_thread* mThread[2];
	

private:
	Ui::ridesharing_sim_qtClass ui;

	void read_file(std::string filename, int iFile);

	std::map<std::string, std::vector<double>> averages[2];
	std::map<std::string, std::vector<double>> stddevs[2];
	std::map<std::string, std::vector<double>> counts[2];


private slots:
	void onButtonSimulateClicked();
	void onButtonStopClicked();
	void onButtonSwapXYClicked();
	void onButtonChooseFile1Clicked();
	void onButtonChooseFile2Clicked();
	void onButtonCompareClicked();
	void onButtonSavePlotClicked();

	void B_to_checkBoxChanged(int arg1);

public slots:
	void onProcessTextChanged(QString newText);
	void onGraphChanged(QVector<double> vB, QVector<double> vE);
	void onErrorMessage(QString messageText);
};
