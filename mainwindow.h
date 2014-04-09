#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QDebug>
#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);

    ~MainWindow();

    void Init();
    void GenTraj();
    void RobotTraj(int traj);
    void Simu();
    double sign(double a);

private:
    Ui::MainWindow *ui;

signals:


private slots:
    void on_ButtonStart_clicked();
    void on_InterSpinBox_valueChanged(int arg1);
    void on_Zoomplus_clicked();
    void on_Zoomminus_clicked();
    void on_ZoomZone_clicked();
    void on_ZoomReset_clicked();
    void repaint();
    void on_ButtonFindSol_clicked();
    void on_ErrSpinBox_1_valueChanged(double arg1);
    void on_ErrSpinBox_2_valueChanged(double arg1);
    void on_ErrSpinBox_3_valueChanged(double arg1);
    void on_ErrSpinBox_4_valueChanged(double arg1);
    void on_ErrSpinBox_5_valueChanged(double arg1);
    void on_EpsilonSpinBox_valueChanged(double arg1);
    void on_checkBox_toggled(bool checked);

    void on_TSlider_valueChanged(int value);
    void on_OutlierSpinBox_valueChanged(double arg1);
    void on_Tplot_clicked();
    void delay();
};



#endif // MAINWINDOW_H
