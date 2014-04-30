#ifndef SIVIA_H
#define SIVIA_H

#include "repere.h"
#include "ibex.h"

using namespace ibex;
using namespace std;

struct sivia_struct{
    double *x;
    double *y;
    double *z;
    double robot_position_found[3];
    double robot_position[4];
    double *theta;
    double *speedx;
    double *speedy;
    int isinside;
    int q;
    int nb_beacon;
    int in_perhaps;
    double *err;
    double epsilon_sivia;
    int *outliers;
    double erroutlier;
    int iteration;
    vector<IntervalVector> box;
    vector<IntervalVector> vin;
    vector<IntervalVector> vin_prev;
};

class Sivia {
public:

    /*
     * Run the SIVIA algorithm.
     *
     * Parameters:
     * R:   where to draw the boxes.
     * epsilon: precision downto which boxes are bisected.
     */
    Sivia(repere& R,struct sivia_struct *my_struct);

    /*
     * Contract "box" with "c" and draw the trace (i.e., the difference between box and c(box))
     * with the colors "pencolor" and "brushcolor".
     */
    void contract_and_draw(Ctc& c, IntervalVector& box,IntervalVector& iinside,int isctcinsside,struct sivia_struct *my_struct,int& nbox,  const QColor & pencolor, const QColor & brushcolor);

private:
    repere& R;
};

#endif // SIVIA_H
