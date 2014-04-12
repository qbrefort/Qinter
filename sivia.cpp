#include "sivia.h"
#include <math.h>

void Sivia::contract_and_draw(Ctc& c, IntervalVector& X,int isctcinside,int& isinside, const QColor & pencolor, const QColor & brushcolor) {
    IntervalVector X0=X;       // get a copy
    try {
        c.contract(X);
        if (X==X0) return;     // nothing contracted.
        IntervalVector* rest;
        int n=X0.diff(X,rest); // calculate the set difference
        for (int i=0; i<n; i++) {     // display the boxes
            R.DrawBox(rest[i][0].lb(),rest[i][0].ub(), rest[i][1].lb(),rest[i][1].ub(),QPen(pencolor),QBrush(brushcolor));
            if (isctcinside==1) isinside=1;
        }
        delete[] rest;
    } catch(EmptyBoxException&) {
        R.DrawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),QPen(pencolor),QBrush(brushcolor));
    }
}

Sivia::Sivia(repere& R,double *rpos,int Qinter,int &isinside,int &Sperhaps,double *err, double epsilon, double erroutlier) : R(R) {

    //min g(x)=sum(err[i])
    //all my constraints

    isinside=0;
    // Create the function we want to apply SIVIA on.
    Variable x,y;
    double x1,y1,x2,y2,x3,y3,x4,y4,x5,y5;
    double r1,r2,r3,r4,r5;
    //r1=9.0;r2=2.0;r3=6.0;r4=6.0;r5=10.0;
    y1=x1=14;x2=4;y2=-7;x3=7;y3=10;x4=y4=-10;x5=-4;y5=12;
    double xr=rpos[0],yr=rpos[1];
    r1=sqrt(pow((xr-x1),2)+pow((yr-y1),2));
    r2=sqrt(pow((xr-x2),2)+pow((yr-y2),2))*(1+erroutlier/100);
    r3=sqrt(pow((xr-x3),2)+pow((yr-y3),2));
    r4=sqrt(pow((xr-x4),2)+pow((yr-y4),2));
    r5=sqrt(pow((xr-x5),2)+pow((yr-y5),2));
    Function f(x,y,sqr(x-x1)+sqr(y-y1));
    Function f2(x,y,sqr(x-x2)+sqr(y-y2));
    Function f3(x,y,sqr(x-x3)+sqr(y-y3));
    Function f4(x,y,sqr(x-x4)+sqr(y-y4));
    Function f5(x,y,sqr(x-x5)+sqr(y-y5));

    NumConstraint c11(x,y,f(x,y)<=(r1+err[0])*(r1+err[0]));
    NumConstraint c12(x,y,f(x,y)>=(r1-err[0])*(r1-err[0]));
    NumConstraint c13(x,y,f(x,y)>(r1+err[0])*(r1+err[0]));
    NumConstraint c14(x,y,f(x,y)<(r1-err[0])*(r1-err[0]));


    NumConstraint c21(x,y,f2(x,y)<=(r2+err[1])*(r2+err[1]));
    NumConstraint c22(x,y,f2(x,y)>=(r2-err[1])*(r2-err[1]));
    NumConstraint c23(x,y,f2(x,y)>(r2+err[1])*(r2+err[1]));
    NumConstraint c24(x,y,f2(x,y)<(r2-err[1])*(r2-err[1]));

    NumConstraint c31(x,y,f3(x,y)<=(r3+err[2])*(r3+err[2]));
    NumConstraint c32(x,y,f3(x,y)>=(r3-err[2])*(r3-err[2]));
    NumConstraint c33(x,y,f3(x,y)>(r3+err[2])*(r3+err[2]));
    NumConstraint c34(x,y,f3(x,y)<(r3-err[2])*(r3-err[2]));

    NumConstraint c41(x,y,f4(x,y)<=(r4+err[3])*(r4+err[3]));
    NumConstraint c42(x,y,f4(x,y)>=(r4-err[3])*(r4-err[3]));
    NumConstraint c43(x,y,f4(x,y)>(r4+err[3])*(r4+err[3]));
    NumConstraint c44(x,y,f4(x,y)<(r4-err[3])*(r4-err[3]));

    NumConstraint c51(x,y,f5(x,y)<=(r5+err[4])*(r5+err[4]));
    NumConstraint c52(x,y,f5(x,y)>=(r5-err[4])*(r5-err[4]));
    NumConstraint c53(x,y,f5(x,y)>(r5+err[4])*(r5+err[4]));
    NumConstraint c54(x,y,f5(x,y)<(r5-err[4])*(r5-err[4]));


    // Create contractors with respect to each
    // of the previous constraints.
    CtcFwdBwd out11(c11);CtcFwdBwd out12(c12);
    CtcFwdBwd in11(c13); CtcFwdBwd in12(c14);

    CtcFwdBwd out21(c21);CtcFwdBwd out22(c22);
    CtcFwdBwd in21(c23);CtcFwdBwd in22(c24);


    CtcFwdBwd out31(c31);CtcFwdBwd out32(c32);
    CtcFwdBwd in31(c33);CtcFwdBwd in32(c34);


    CtcFwdBwd out41(c41);CtcFwdBwd out42(c42);
    CtcFwdBwd in41(c43);CtcFwdBwd in42(c44);


    CtcFwdBwd out51(c51);CtcFwdBwd out52(c52);
    CtcFwdBwd in51(c53);CtcFwdBwd in52(c54);


    // Create a contractor that removes all the points
    // that do not satisfy either f(x,y)<=(r+epsilon)^2 or f(x,y)>=(r-epsilon)^2.
    // These points are "outside" of the solution set.
    CtcCompo outside1(out11,out12);
    CtcCompo outside2(out21,out22);
    CtcCompo outside3(out31,out32);
    CtcCompo outside4(out41,out42);
    CtcCompo outside5(out51,out52);

    // Create a contractor that removes all the points
    // that do not satisfy both f(x,y)>(r+epsilon)^2 or f(x,y)<(r-epsilon)^2.
    // These points are "inside" the solution set.
    CtcUnion inside1(in11,in12);
    CtcUnion inside2(in21,in22);
    CtcUnion inside3(in31,in32);
    CtcUnion inside4(in41,in42);
    CtcUnion inside5(in51,in52);

    Array<Ctc> Ain(inside1,inside2,inside3,inside4,inside5);

    Array<Ctc> Aout(outside1,outside2,outside3,outside4,outside5);
    int maxq = 5; //nb of contractors
    int ctcq = maxq - Qinter + 1; //nb for q-relaxed function of Ibex
    CtcQInter inside(Ain,ctcq);
    CtcQInter outside(Aout,Qinter);


    // Build the initial box.
    IntervalVector box(2);
    box[0]=Interval(-25,25);
    box[1]=Interval(-25,25);

    // Build the way boxes will be bisected.
    // "LargestFirst" means that the dimension bisected
    // is always the largest one.
    LargestFirst lf;

    stack<IntervalVector> s;
    s.push(box);
    int q = maxq - Qinter;
    int K_ok[5]={0,0,0,0,0};
    int K_ind[5]={0,0,0,0,0};
    Sperhaps=0;
    while (!s.empty()) {
        IntervalVector box=s.top();

        s.pop();
        contract_and_draw(inside,box,1,isinside,Qt::magenta,Qt::red);
        if (box.is_empty()) { continue; }

        contract_and_draw(outside,box,0,isinside,Qt::darkBlue,Qt::cyan);
        if (box.is_empty()) { continue; }

        if (box.max_diam()<epsilon) {
            R.DrawBox(box[0].lb(),box[0].ub(),box[1].lb(),box[1].ub(),QPen(Qt::yellow),QBrush(Qt::white));
            Sperhaps=1;
        } else {
            pair<IntervalVector,IntervalVector> boxes=lf.bisect(box);
            s.push(boxes.first);
            s.push(boxes.second);
        }
    }

    double re=0.5;
    // Draw the beacons
    R.DrawEllipse(x1,y1,re,QPen(Qt::black),QBrush(Qt::NoBrush));
    R.DrawEllipse(x2,y2,re,QPen(Qt::black),QBrush(Qt::NoBrush));
    R.DrawEllipse(x3,y3,re,QPen(Qt::black),QBrush(Qt::NoBrush));
    R.DrawEllipse(x4,y4,re,QPen(Qt::black),QBrush(Qt::NoBrush));
    R.DrawEllipse(x5,y5,re,QPen(Qt::black),QBrush(Qt::NoBrush));

}
