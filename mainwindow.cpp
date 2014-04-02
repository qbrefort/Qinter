#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "ibex.h"
#include "sivia.h"

double epsilon;
double Qinter;
int bfind=0;

MainWindow::MainWindow(QWidget *parent) :  QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    Init();

    // Build the frame
    double xmin=-10;
    double xmax=10;
    double ymin=-10;
    double ymax=10;

    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);
    Sivia sivia(*R,Qinter,bfind,epsilon);
    QString valueAsString = "Qinter = ";
    valueAsString.append(QString::number(Qinter));
    QMessageBox::information(this,"Info",valueAsString);
    // run SIVIA
    while(bfind==0){
        Qinter=Qinter-1;
        QString valueAsString = "Qinter = ";
        valueAsString.append(QString::number(Qinter));
        QMessageBox::information(this,"Info",valueAsString);
        Sivia sivia(*R,Qinter,bfind,epsilon);
    }

}

void MainWindow::Init() {
    epsilon=ui->EpsilonSpinBox->value();
    Qinter=ui->InterSpinBox->value();
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_ButtonStart_clicked()
{

    Init();

    // Build the frame
    double xmin=-10;
    double xmax=10;
    double ymin=-10;
    double ymax=10;

    repere* R = new repere(this,ui->graphicsView,xmin,xmax,ymin,ymax);
    // run SIVIA
    Sivia sivia(*R,Qinter,bfind,epsilon);
}


void MainWindow::on_InterSpinBox_valueChanged(int arg1)
{

}
