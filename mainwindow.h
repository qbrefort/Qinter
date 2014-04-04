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

private:
    Ui::MainWindow *ui;

signals:


private slots:
    void on_ButtonStart_clicked();
    void on_InterSpinBox_valueChanged(int arg1);
    void on_ErrSpinBox_2_valueChanged(double arg1);
    void on_Zoomplus_clicked();
    void on_Zoomminus_clicked();
    void on_ZoomZone_clicked();
    void on_ZoomReset_clicked();
    void repaint();
};



#endif // MAINWINDOW_H
