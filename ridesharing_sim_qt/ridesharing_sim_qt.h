#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_ridesharing_sim_qt.h"
#include "simulation_thread.h"


class ridesharing_sim_qt : public QMainWindow
{
	Q_OBJECT

public:
	ridesharing_sim_qt(QWidget *parent = Q_NULLPTR);
	simulation_thread* mThread;

private:
	Ui::ridesharing_sim_qtClass ui;

private slots:
	void ifkl();
	void ifklstopped();
	void B_to_checkBoxChanged(int arg1);
	void SaveCheckBoxChanged(int arg1);

public slots:
	void onProcessTextChanged(QString newText);
	void onGraphChanged(QVector<double> vB, QVector<double> vE);
};
