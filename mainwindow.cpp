#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "ibex.h"
#include "sivia.h"

#include <QElapsedTimer>
vector<double> u,xv,yv,thetav;
double t;
int timeinfo=1;
double xmin=-25,xmax=25,ymin=-25,ymax=25;
double epsilon,erroutlier;
double err[5]={0.5,0.5,0.5,0.5,0.5};
double Qinter;
int bfind=0;
double rpos[3]={3,-5,0};

MainWindow::MainWindow(QWidget *parent) :  QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    GenTraj();
    //MainWindow::on_ButtonFindSol_clicked();
}

void MainWindow::Init() {
    epsilon=ui->EpsilonSpinBox->value();
    Qinter=ui->InterSpinBox->value();
}

double MainWindow::sign(double a){
    if (a >= 0)
        return 1;
    else
        return -1;
}
void MainWindow::GenTraj(){
    int nb_step=6500;
    xv.resize(nb_step); yv.resize(nb_step); thetav.resize(nb_step); u.resize(nb_step);
    xv[0] = 0; yv[0] = 0; thetav[0]=0;
    double dt=0.02;
    double w=0.05;
    double tt=0;
    for(uint i= 1; i < nb_step; i++){
        tt = tt+dt;
        u[i-1]= 0.1*sign(sin(w*tt));
        yv[i] = yv[i-1] + dt*sin(thetav[i-1]);
        xv[i] = xv[i-1] + dt*cos(thetav[i-1]);
        thetav[i] = thetav[i-1] + dt*(u[i-1]);
    }
}

void MainWindow::RobotTraj(){
    //double thetav= (atan(sqrt(7)-4*cos(t))+atan(sqrt(7)+4*cos(t)))-(atan(4-sqrt(7))+atan(4+sqrt(7)));
    // 8 curve
    rpos[0]=xv[t];
    rpos[1]=yv[t];
    rpos[2]=thetav[t];

}
void MainWindow::Simu(){
    ui->checkBox->setChecked(false);

    for(double i=0;i<6500;i=i+200){
        t=i;

        on_ButtonFindSol_clicked();
        ui->TSlider->setValue(t);
        Zoom();
        delay();
    }
    QString messend = "End of simulation";
    QMessageBox::information(this,"Info",messend);

}

void MainWindow::repaint()
{
    RobotTraj();
    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);
    Sivia sivia(*R,rpos,Qinter,bfind,err,epsilon,erroutlier);
    for(double i=0;i<6500-111;i=i+10){
        R->DrawLine(xv[i],yv[i],xv[i+10],yv[i+10],QPen(Qt::red));
    }

    R->DrawRobot(rpos[0],rpos[1],rpos[2]);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_ButtonFindSol_clicked()
{
    QElapsedTimer timer;
    timer.start();
    RobotTraj();
    Init();
    // Build the frame
    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);


    if(Qinter==5){
        for (uint i=0;i<(sizeof(err)/sizeof(*err));i++){
           err[i] = 0.00;
        }
        Sivia sivia(*R,rpos,Qinter,bfind,err,epsilon,erroutlier);
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
    //        qDebug()<<"Start Step: "<<startstep<<endl;
    //        qDebug()<<"Step: "<<step<<endl;
            // "Backward"
            // The idea here is to developp a 'forward/backward' like method.
            // Maybe we iterate with a high step in forward and lower the step to find a solution in backward.
            while(bfind==1){

                for (uint j=0;j<(sizeof(err)/sizeof(*err));j++){
                    err[j]=startstep;
                    if(i==j) err[j]=startstep-((stepctr+1))*step;
                }

                Sivia sivia(*R,rpos,Qinter,bfind,err,epsilon,erroutlier);
                stepctr=(stepctr+1)%nstep;
                if (stepctr==0){
                    i++;
                    i = i % (sizeof(err)/sizeof(*err));
                    if(i==0)    startstep=startstep-step;
                }
                back++;
            }
           // qDebug()<<"err back: "<<"is "<<err[0]<<";"<<err[1]<<";"<<err[2]<<";"<<err[3]<<";"<<err[4]<<endl;

            while(bfind==0){

                for (uint j=0;j<(sizeof(err)/sizeof(*err));j++){
                    err[j]=startstep;
                    if(i==j) err[j]=startstep+((stepctr+1))*step;
                }

                Sivia sivia(*R,rpos,Qinter,bfind,err,epsilon,erroutlier);
                stepctr=(stepctr+1)%nstep;
                if (stepctr==0){
                    i++;
                    i = i % (sizeof(err)/sizeof(*err));
                    if(i==0)    startstep=startstep+step;
                }
                forw++;
            }
            //qDebug()<<"err for: "<<"is "<<err[0]<<";"<<err[1]<<";"<<err[2]<<";"<<err[3]<<";"<<err[4]<<endl;
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
        repaint();

        if (timeinfo){
            QString mess = "Execution time : ";
            mess.append(QString::number(timer.elapsed()));mess.append(" ms");
            QMessageBox::information(this,"Info",mess);
        }
    }
    else{
        on_ButtonStartParam_clicked();
    }
}


void MainWindow::on_ButtonStartParam_clicked()
{
    QElapsedTimer timer;
    timer.start();
    RobotTraj();
    Init();
    if(Qinter!=5){
        ui->ErrSpinBox_1->setValue(0.3);
        ui->ErrSpinBox_2->setValue(0.3);
        ui->ErrSpinBox_3->setValue(0.3);
        ui->ErrSpinBox_4->setValue(0.3);
        ui->ErrSpinBox_5->setValue(0.3);
    }

    // Build the frame

    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);
    // run SIVIA
    Sivia sivia(*R,rpos,Qinter,bfind,err,epsilon,erroutlier);
    R->DrawRobot(rpos[0],rpos[1],rpos[2]);
    repaint();
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
void MainWindow::Zoom()
{
    xmin = xv[t+220]-5;
    xmax = xv[t+220]+5;
    ymin = yv[t+220]-5;
    ymax = yv[t+220]+5;
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
    xmin = rpos[0]-5;
    xmax = rpos[0]+5;
    ymin = rpos[1]-5;
    ymax = rpos[1]+5;
    repaint();
}

void MainWindow::on_ZoomReset_clicked()
{
    xmin = -25;
    xmax = 25;
    ymin = -25;
    ymax = 25;
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

void MainWindow::on_TSlider_valueChanged(int value)
{
    t=double(value);
}

void MainWindow::on_OutlierSpinBox_valueChanged(double arg1)
{
    erroutlier = arg1;
}

void MainWindow::on_Tplot_clicked()
{
    Simu();
}
void MainWindow::delay()
{
    QTime dieTime= QTime::currentTime().addMSecs(10);
    while( QTime::currentTime() < dieTime )
    QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
}



