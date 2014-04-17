#include "sivia.h"
#include <math.h>

#include <stdlib.h>

void Sivia::contract_and_draw(Ctc& c, IntervalVector& X,IntervalVector& iinside,int isctcinside,int& isinside, const QColor & pencolor, const QColor & brushcolor) {
    IntervalVector X0= X;       // get a copy
    try {
        c.contract(X);
        if (X==X0) return;     // nothing contracted.
        IntervalVector* rest;
        int n=X0.diff(X,rest); // calculate the set difference
        for (int i=0; i<n; i++) {     // display the boxes
            R.DrawBox(rest[i][0].lb(),rest[i][0].ub(), rest[i][1].lb(),rest[i][1].ub(),QPen(pencolor),QBrush(brushcolor));
            if (isctcinside==1) {
                isinside=1;
                iinside=rest[i].mid();
            }
        }
        delete[] rest;
    } catch(EmptyBoxException&) {
        R.DrawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),QPen(pencolor),QBrush(brushcolor));
    }
}

Sivia::Sivia(repere& R,IntervalVector& iinside,double *rpos,int Qinter,int nbeacon,int &isinside,int &Sperhaps,double *err, double epsilon,int *outlier, double erroutlier) : R(R) {

    //min g(x)=sum(err[i])
    //all my constraints
    isinside=0;
    Variable xvar,yvar;
    int n = nbeacon;
    double *x=new double[n]; // vecteur des abcisses des donnees
    double *y=new double[n]; // vecteur des ordonnees des donnees
    double *r=new double[n]; // vecteur des rayons
    //r1=9.0;r2=2.0;r3=6.0;r4=6.0;r5=10.0;
    //random config (original config)
    //y1=x1=14;x2=4;y2=-7;x3=7;y3=10;x4=y4=-10;x5=-4;y5=12;x6=0;y6=0;
    //square config
//    x1=-24;y1=24;x2=-24;y2=-24;x3=24;y3=24;x4=24;y4=-24;x5=y5=0;x6=y6=0;
    //square + 1 close
    if (n>=2){
        x[0]=-24;y[0]=24;x[1]=-24;y[1]=22;
    }
    if (n>=3){
        x[2]=-24;y[2]=-24;
    }
    if (n>=4){
        x[3]=24;y[3]=24;
    }
    if (n>=5){
        x[4]=24;y[4]=-24;
    }
    if (n>=6){
        x[5]=y[5]=0;
    }
    if (n>=7){
        for(int i=6;i<n;i++){
            x[i]= 25 - rand() % 50;
            y[i]= 25 - rand() % 50;
        }
    }
    double xr=rpos[0],yr=rpos[1];

    for (int i=0;i<n;i++) {
        r[i]= sqrt(pow(xr-x[i],2)+pow(yr-y[i],2));
        if (outlier[i]==1)
            r[i] *= (1+erroutlier/100);
    }

    vector<Function*> f;
    for(int i=0;i<n;i++) {
        f.push_back(new Function(xvar,yvar,sqrt(sqr(xvar-x[i])+sqr(yvar-y[i]))));
    }

    vector<Ctc*> vec_out;
    vector<Ctc*> vec_in;
    for(int i=0;i<n;i++) {
        vec_out.push_back(new CtcIn(*(f[i]),(r[i]+Interval(-1,1)*err[i])));
        vec_in.push_back(new CtcNotIn(*(f[i]),(r[i]+Interval(-1,1)*err[i])));
    }

    double re=0.5;
    int maxq = nbeacon; //nb of contractors
    int ctcq = maxq - Qinter + 1; //nb for q-relaxed function of Ibex

    CtcQInter inside(vec_in,ctcq);
    CtcQInter outside(vec_out,Qinter);
    IntervalVector box(2, Interval(-25,25));

    LargestFirst lf;
    stack<IntervalVector> s;
    s.push(box);
    Sperhaps=0;
    while (!s.empty()) {
        IntervalVector box=s.top();

        s.pop();
        contract_and_draw(inside,box,iinside,1,isinside,Qt::magenta,Qt::red);
        if (box.is_empty()) { continue; }

        contract_and_draw(outside,box,iinside,0,isinside,Qt::darkBlue,Qt::cyan);
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
    for(int i=0;i<n;i++)
        R.DrawEllipse(x[i],y[i],re,QPen(Qt::black),QBrush(Qt::NoBrush));


    vec_out.clear();
    vec_in.clear();
    f.clear();
}

