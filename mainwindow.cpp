#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "ibex.h"
#include "sivia.h"

#include <QElapsedTimer>

double xmin=-10,xmax=10,ymin=-10,ymax=10;
double epsilon;
double err[5]={0.5,0.5,0.5,0.5,0.5};
double Qinter;
int bfind=0;

MainWindow::MainWindow(QWidget *parent) :  QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    //MainWindow::on_ButtonFindSol_clicked();
}

void MainWindow::Init() {
    epsilon=ui->EpsilonSpinBox->value();
    Qinter=ui->InterSpinBox->value();
}

void MainWindow::repaint()
{
    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);
    Sivia sivia(*R,Qinter,bfind,err,epsilon);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_ButtonFindSol_clicked()
{
    QElapsedTimer timer;
    timer.start();
    Init();

    // Build the frame

    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);
    Sivia sivia(*R,Qinter,bfind,err,epsilon);
    // run SIVIA
//    while(bfind==0){
//        Qinter=Qinter-1;
//        ui->InterSpinBox->setValue(Qinter);
//        Sivia sivia(*R,Qinter,bfind,err,epsilon);
//    }
    //QString valueAsString = "Robot localized with ";
    //valueAsString.append(QString::number(Qinter));valueAsString.append("/5 inliers");
    //QMessageBox::information(this,"Info",valueAsString);

    while(bfind==1){
        for (uint i=0;i<(sizeof(err)/sizeof(*err));i++){
            err[i]=err[i]-0.05;
            Sivia sivia(*R,Qinter,bfind,err,epsilon);
        }

    }
    //err=err+0.05;
//    ui->ErrSpinBox_2->setValue(err[0]);
//    Sivia sivia2(*R,Qinter,bfind,err,epsilon);
    //QString mess2 = "Robot localized with ";
    //mess2.append(QString::number(err));mess2.append(" error");
    //QMessageBox::information(this,"Info",mess2);
    qDebug() << "This operation took" << timer.elapsed() << "milliseconds";
}


void MainWindow::on_ButtonStart_clicked()
{
    QElapsedTimer timer;
    timer.start();
    Init();
    for (uint i=0;i<(sizeof(err)/sizeof(*err));i++){
       err[i] = ui->ErrSpinBox_2->value();
    }
    // Build the frame

    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);
    // run SIVIA
    Sivia sivia(*R,Qinter,bfind,err,epsilon);
    qDebug() << "This operation took" << timer.elapsed() << "milliseconds";
}


void MainWindow::on_ErrSpinBox_2_valueChanged(double arg1)
{
    err[0] = arg1;
}

void MainWindow::on_InterSpinBox_valueChanged(int arg1){}

void MainWindow::on_Zoomplus_clicked()
{
    xmin /= 2;
    xmax /= 2;
    ymin /= 2;
    ymax /= 2;
    repaint();
}

void MainWindow::on_Zoomminus_clicked()
{
    xmin *= 2;
    xmax *= 2;
    ymin *= 2;
    ymax *= 2;
    repaint();
}

void MainWindow::on_ZoomZone_clicked()
{
    xmin = 2.5;
    xmax = 3.5;
    ymin = -6;
    ymax = -4;
    repaint();
}

void MainWindow::on_ZoomReset_clicked()
{
    xmin = -10;
    xmax = 10;
    ymin = -10;
    ymax = 10;
    repaint();
}

