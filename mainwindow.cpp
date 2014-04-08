#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "ibex.h"
#include "sivia.h"

#include <QElapsedTimer>

int timeinfo=1;
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
    for (uint i=0;i<(sizeof(err)/sizeof(*err));i++){
       err[i] = 0.00;
    }

    // Build the frame

    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);
    Sivia sivia(*R,Qinter,bfind,err,epsilon);

    uint i=0;
    //double startstep=0.05+floor(10*epsilon)/10-floor(10*epsilon)/20;
    double startstep=1;
    //qDebug()<<"start: "<<startstep<<endl;
    double step = 0.5;
    int nstep = 2;
    int stepctr=0;

    while(step>0.05){
        int forw=0;
        int back=0;
        // "Forward"
        qDebug()<<"Start Step: "<<startstep<<endl;
        qDebug()<<"Step: "<<step<<endl;
        // "Backward"
        // The idea here is to developp a 'forward/backward' like method.
        // Maybe we iterate with a high step in forward and lower the step to find a solution in backward.
        while(bfind==1){

            for (uint j=0;j<(sizeof(err)/sizeof(*err));j++){
                err[j]=startstep;
                if(i==j) err[j]=startstep-((stepctr+1))*step;
            }

            Sivia sivia(*R,Qinter,bfind,err,epsilon);
            stepctr=(stepctr+1)%nstep;
            if (stepctr==0){
                i++;
                i = i % (sizeof(err)/sizeof(*err));
                if(i==0)    startstep=startstep-step;
            }
            back++;
        }
        qDebug()<<"err back: "<<"is "<<err[0]<<";"<<err[1]<<";"<<err[2]<<";"<<err[3]<<";"<<err[4]<<endl;

        while(bfind==0){

            for (uint j=0;j<(sizeof(err)/sizeof(*err));j++){
                err[j]=startstep;
                if(i==j) err[j]=startstep+((stepctr+1))*step;
            }

            Sivia sivia(*R,Qinter,bfind,err,epsilon);
            stepctr=(stepctr+1)%nstep;
            if (stepctr==0){
                i++;
                i = i % (sizeof(err)/sizeof(*err));
                if(i==0)    startstep=startstep+step;
            }
            forw++;
        }
        qDebug()<<"err for: "<<"is "<<err[0]<<";"<<err[1]<<";"<<err[2]<<";"<<err[3]<<";"<<err[4]<<endl;
        if(back>forw)
            startstep/=0.5;
        else
            startstep*=0.5;
        step/=2;

    }
    ui->ErrSpinBox_1->setValue(err[0]);
    ui->ErrSpinBox_2->setValue(err[1]);
    ui->ErrSpinBox_3->setValue(err[2]);
    ui->ErrSpinBox_4->setValue(err[3]);
    ui->ErrSpinBox_5->setValue(err[4]);

    if (timeinfo){
        QString mess = "Execution time : ";
        mess.append(QString::number(timer.elapsed()));mess.append(" ms");
        QMessageBox::information(this,"Info",mess);
    }
}


void MainWindow::on_ButtonStart_clicked()
{
    QElapsedTimer timer;
    timer.start();
    Init();

    // Build the frame

    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);
    // run SIVIA
    Sivia sivia(*R,Qinter,bfind,err,epsilon);
    if (timeinfo){
        QString mess = "Execution time : ";
        mess.append(QString::number(timer.elapsed()));mess.append(" ms");
        QMessageBox::information(this,"Info",mess);
    }
}

void MainWindow::on_ErrSpinBox_1_valueChanged(double arg1)
{
    err[0] = arg1;
}

void MainWindow::on_ErrSpinBox_2_valueChanged(double arg1)
{
    err[1] = arg1;
}
void MainWindow::on_ErrSpinBox_3_valueChanged(double arg1)
{
    err[2] = arg1;
}
void MainWindow::on_ErrSpinBox_4_valueChanged(double arg1)
{
    err[3] = arg1;
}
void MainWindow::on_ErrSpinBox_5_valueChanged(double arg1)
{
    err[4] = arg1;
}
void MainWindow::on_InterSpinBox_valueChanged(int arg1)
{
    Qinter = arg1;
}

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
    xmin = 3.0;
    xmax = 3.75;
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


void MainWindow::on_EpsilonSpinBox_valueChanged(double arg1)
{
    epsilon = arg1;
}


void MainWindow::on_checkBox_toggled(bool checked)
{
    if(checked) timeinfo=1;
    else timeinfo=0;
}
