/********************************************************************************
** Form generated from reading UI file 'ridesharing_sim_qt.ui'
**
** Created by: Qt User Interface Compiler version 5.13.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_RIDESHARING_SIM_QT_H
#define UI_RIDESHARING_SIM_QT_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_ridesharing_sim_qtClass
{
public:
    QWidget *centralWidget;
    QPushButton *pushButton;
    QListWidget *listWidget;
    QTextEdit *textEdit;
    QTextEdit *textEdit_2;
    QTextEdit *textEdit_3;
    QTextEdit *textEdit_4;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QLabel *label_4;
    QCheckBox *checkBox;
    QLabel *label_5;
    QTextEdit *textEdit_5;
    QLabel *label_6;
    QCheckBox *checkBox_2;
    QCustomPlot *customPlot;
    QTextEdit *textEdit_6;
    QPushButton *pushButton_2;
    QTextEdit *textEdit_7;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *ridesharing_sim_qtClass)
    {
        if (ridesharing_sim_qtClass->objectName().isEmpty())
            ridesharing_sim_qtClass->setObjectName(QString::fromUtf8("ridesharing_sim_qtClass"));
        ridesharing_sim_qtClass->resize(1186, 758);
        centralWidget = new QWidget(ridesharing_sim_qtClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        pushButton = new QPushButton(centralWidget);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        pushButton->setGeometry(QRect(880, 570, 121, 31));
        listWidget = new QListWidget(centralWidget);
        listWidget->setObjectName(QString::fromUtf8("listWidget"));
        listWidget->setGeometry(QRect(880, 20, 256, 192));
        textEdit = new QTextEdit(centralWidget);
        textEdit->setObjectName(QString::fromUtf8("textEdit"));
        textEdit->setGeometry(QRect(880, 230, 251, 31));
        textEdit_2 = new QTextEdit(centralWidget);
        textEdit_2->setObjectName(QString::fromUtf8("textEdit_2"));
        textEdit_2->setGeometry(QRect(880, 280, 251, 31));
        textEdit_3 = new QTextEdit(centralWidget);
        textEdit_3->setObjectName(QString::fromUtf8("textEdit_3"));
        textEdit_3->setGeometry(QRect(880, 380, 251, 31));
        textEdit_4 = new QTextEdit(centralWidget);
        textEdit_4->setObjectName(QString::fromUtf8("textEdit_4"));
        textEdit_4->setGeometry(QRect(880, 430, 251, 31));
        label = new QLabel(centralWidget);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(760, 20, 81, 31));
        label_2 = new QLabel(centralWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(760, 230, 81, 31));
        label_3 = new QLabel(centralWidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(760, 270, 81, 31));
        label_4 = new QLabel(centralWidget);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(760, 380, 81, 31));
        checkBox = new QCheckBox(centralWidget);
        checkBox->setObjectName(QString::fromUtf8("checkBox"));
        checkBox->setGeometry(QRect(890, 330, 191, 31));
        label_5 = new QLabel(centralWidget);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(720, 430, 161, 31));
        textEdit_5 = new QTextEdit(centralWidget);
        textEdit_5->setObjectName(QString::fromUtf8("textEdit_5"));
        textEdit_5->setGeometry(QRect(880, 510, 251, 31));
        label_6 = new QLabel(centralWidget);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(760, 510, 101, 31));
        checkBox_2 = new QCheckBox(centralWidget);
        checkBox_2->setObjectName(QString::fromUtf8("checkBox_2"));
        checkBox_2->setGeometry(QRect(890, 470, 191, 31));
        customPlot = new QCustomPlot(centralWidget);
        customPlot->setObjectName(QString::fromUtf8("customPlot"));
        customPlot->setGeometry(QRect(19, 19, 691, 671));
        textEdit_6 = new QTextEdit(centralWidget);
        textEdit_6->setObjectName(QString::fromUtf8("textEdit_6"));
        textEdit_6->setGeometry(QRect(880, 620, 251, 91));
        pushButton_2 = new QPushButton(centralWidget);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));
        pushButton_2->setGeometry(QRect(1010, 570, 121, 31));
        textEdit_7 = new QTextEdit(centralWidget);
        textEdit_7->setObjectName(QString::fromUtf8("textEdit_7"));
        textEdit_7->setGeometry(QRect(1020, 330, 111, 31));
        ridesharing_sim_qtClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(ridesharing_sim_qtClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1186, 18));
        ridesharing_sim_qtClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(ridesharing_sim_qtClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        ridesharing_sim_qtClass->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(ridesharing_sim_qtClass);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        ridesharing_sim_qtClass->setStatusBar(statusBar);

        retranslateUi(ridesharing_sim_qtClass);

        QMetaObject::connectSlotsByName(ridesharing_sim_qtClass);
    } // setupUi

    void retranslateUi(QMainWindow *ridesharing_sim_qtClass)
    {
        ridesharing_sim_qtClass->setWindowTitle(QCoreApplication::translate("ridesharing_sim_qtClass", "ridesharing_sim_qt", nullptr));
        pushButton->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "Simulate", nullptr));
        label->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "Topology", nullptr));
        label_2->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "N", nullptr));
        label_3->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "B", nullptr));
        label_4->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "until B", nullptr));
        checkBox->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "Num_B:", nullptr));
        label_5->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "Normalized request rate", nullptr));
        label_6->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "Filename", nullptr));
        checkBox_2->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "Save", nullptr));
        pushButton_2->setText(QCoreApplication::translate("ridesharing_sim_qtClass", "Stop Simulation", nullptr));
    } // retranslateUi

};

namespace Ui {
    class ridesharing_sim_qtClass: public Ui_ridesharing_sim_qtClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RIDESHARING_SIM_QT_H
