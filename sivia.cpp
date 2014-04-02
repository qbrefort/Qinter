#include "sivia.h"


void Sivia::contract_and_draw(Ctc& c, IntervalVector& X,int& bfind, const QColor & pencolor, const QColor & brushcolor) {
    IntervalVector X0=X;       // get a copy
    try {
        c.contract(X);
        if (X==X0) return;     // nothing contracted.
        IntervalVector* rest;
        int n=X0.diff(X,rest); // calculate the set difference
        for (int i=0; i<n; i++) {     // display the boxes
            R.DrawBox(rest[i][0].lb(),rest[i][0].ub(), rest[i][1].lb(),rest[i][1].ub(),QPen(pencolor),QBrush(brushcolor));
            bfind=1;
        }
        delete[] rest;
    } catch(EmptyBoxException&) {
        R.DrawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),QPen(pencolor),QBrush(brushcolor));
    }
}


Sivia::Sivia(repere& R,int Qinter,int &bfind, double epsilon) : R(R) {
    qDebug("Run SIVIA");
    bfind=0;
    // Create the function we want to apply SIVIA on.
    Variable x,y;
    int x1,y1,x2,y2,x3,y3,x4,y4,x5,y5;
    y1=x1=4;x2=4;y2=-7;x3=7;y3=0;x4=y4=-2;x5=y5=6;
    Function f(x,y,sqr(x-x1)+sqr(y-y1));
    Function f2(x,y,sqr(x-x2)+sqr(y-y2));
    Function f3(x,y,sqr(x-x3)+sqr(y-y3));
    Function f4(x,y,sqr(x-x4)+sqr(y-y4));
    Function f5(x,y,sqr(x-x5)+sqr(y-y5));

    NumConstraint c11(x,y,f(x,y)<=9.5*9.5);
    NumConstraint c12(x,y,f(x,y)>=8.5*8.5);
    NumConstraint c13(x,y,f(x,y)>9.5*9.5);
    NumConstraint c14(x,y,f(x,y)<8.5*8.5);


    NumConstraint c21(x,y,f2(x,y)<=2.5*2.5);
    NumConstraint c22(x,y,f2(x,y)>=1.5*1.5);
    NumConstraint c23(x,y,f2(x,y)>2.5*2.5);
    NumConstraint c24(x,y,f2(x,y)<1.5*1.5);

    NumConstraint c31(x,y,f3(x,y)<=6.5*6.5);
    NumConstraint c32(x,y,f3(x,y)>=5.5*5.5);
    NumConstraint c33(x,y,f3(x,y)>6.5*6.5);
    NumConstraint c34(x,y,f3(x,y)<5.5*5.5);

    NumConstraint c41(x,y,f4(x,y)<=6.5*6.5);
    NumConstraint c42(x,y,f4(x,y)>=5.5*5.5);
    NumConstraint c43(x,y,f4(x,y)>6.5*6.5);
    NumConstraint c44(x,y,f4(x,y)<5.5*5.5);

    NumConstraint c51(x,y,f5(x,y)<=9.5*9.5);
    NumConstraint c52(x,y,f5(x,y)>=8.5*8.5);
    NumConstraint c53(x,y,f5(x,y)>9.5*9.5);
    NumConstraint c54(x,y,f5(x,y)<8.5*8.5);


    // Create contractors with respect to each
    // of the previous constraints.
    CtcFwdBwd out11(c11);
    CtcFwdBwd out12(c12);
    CtcFwdBwd in11(c13);
    CtcFwdBwd in12(c14);

    CtcFwdBwd out21(c21);
    CtcFwdBwd out22(c22);
    CtcFwdBwd in21(c23);
    CtcFwdBwd in22(c24);

    CtcFwdBwd out31(c31);
    CtcFwdBwd out32(c32);
    CtcFwdBwd in31(c33);
    CtcFwdBwd in32(c34);

    CtcFwdBwd out41(c41);
    CtcFwdBwd out42(c42);
    CtcFwdBwd in41(c43);
    CtcFwdBwd in42(c44);

    CtcFwdBwd out51(c51);
    CtcFwdBwd out52(c52);
    CtcFwdBwd in51(c53);
    CtcFwdBwd in52(c54);

    // Create a contractor that removes all the points
    // that do not satisfy either f(x,y)<=2 or f(x,y)>=0.
    // These points are "outside" of the solution set.
    CtcCompo outside1(out11,out12);
    CtcCompo outside2(out21,out22);
    CtcCompo outside3(out31,out32);
    CtcCompo outside4(out41,out42);
    CtcCompo outside5(out51,out52);

    // Create a contractor that removes all the points
    // that do not satisfy both f(x,y)>2 or f(x,y)<0.
    // These points are "inside" the solution set.
    CtcUnion inside1(in11,in12);
    CtcUnion inside2(in21,in22);
    CtcUnion inside3(in31,in32);
    CtcUnion inside4(in41,in42);
    CtcUnion inside5(in51,in52);

    Array<Ctc> Ain(inside1,inside2,inside3,inside4,inside5);

    Array<Ctc> Aout(outside1,outside2,outside3,outside4,outside5);
    //int n = 1; //nb of intersected contractor
    int q = 5 - Qinter + 1;
    CtcQInter inside(Ain,q);
    CtcQInter outside(Aout,Qinter);

//    CtcQInter outside(outside1,outside3);
//    CtcQInter inside(inside1,inside3);


    // Build the initial box.
    IntervalVector box(2);
    box[0]=Interval(-10,10);
    box[1]=Interval(-10,10);

    // Build the way boxes will be bisected.
    // "LargestFirst" means that the dimension bisected
    // is always the largest one.
    LargestFirst lf;

    stack<IntervalVector> s;
    s.push(box);

    while (!s.empty()) {
        IntervalVector box=s.top();
        s.pop();
        contract_and_draw(inside,box,bfind,Qt::magenta,Qt::red);
        if (box.is_empty()) { continue; }

        contract_and_draw(outside,box,bfind,Qt::darkBlue,Qt::cyan);
        if (box.is_empty()) { continue; }

        if (box.max_diam()<epsilon) {
            R.DrawBox(box[0].lb(),box[0].ub(),box[1].lb(),box[1].ub(),QPen(Qt::yellow),QBrush(Qt::NoBrush));
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


    R.Save("paving");
}
