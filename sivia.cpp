#include "sivia.h"
#include <math.h>
#include <algorithm>

#include <stdlib.h>



unsigned Sivia::nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

void Sivia::contract_and_draw(Ctc& c, IntervalVector& X,IntervalVector& viinside,int isctcinside,struct sivia_struct *my_struct,int& nbox, const QColor & pencolor, const QColor & brushcolor) {
    IntervalVector X0= X;       // get a copy
    try {
        c.contract(X);
        if (X==X0) return;     // nothing contracted.
        IntervalVector* rest;
        int n=X0.diff(X,rest); // calculate the set difference
        for (int i=0; i<n; i++) {     // display the boxes
            R.DrawBox(rest[i][0].lb(),rest[i][0].ub(), rest[i][1].lb(),rest[i][1].ub(),QPen(pencolor),QBrush(brushcolor));
            if (isctcinside==1) {
                viinside = rest[i];
                my_struct->vin.push_back(viinside);
                nbox++;
                my_struct->isinside=1;
                my_struct->areain += viinside[0].diam()*viinside[1].diam();
                //cout<<viinside<<endl;
            }
            if (isctcinside==0) {
                my_struct->vper.push_back(rest[i]);
            }
        }
        delete[] rest;
    } catch(EmptyBoxException&) {
        R.DrawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),QPen(pencolor),QBrush(brushcolor));
    }
}

void Sivia::Sivia_Pair(struct sivia_struct *my_struct){
    my_struct->areain = 0;
    my_struct->areap = 0;
    my_struct->isinside=0;
    Variable xvar,yvar;
    int n = my_struct->nb_beacon;
    double *x=my_struct->x; // vecteur des abcisses des donnees
    double *y=my_struct->y; // vecteur des ordonnees des donnees
    double *r=new double[n]; // vecteur des rayons
    double xr=0,yr=0;
    xr=my_struct->robot_position[0];
    yr=my_struct->robot_position[1];
    for (int i=0;i<n;i++) {
//        r[i]= sqrt(pow(xr-x[i],2)+pow(yr-y[i],2)+pow(zr-z[i],2));
        r[i]= sqrt(pow(xr-x[i],2)+pow(yr-y[i],2));
        if (my_struct->outliers[i]!=0)
            r[i] *= (1+my_struct->outliers[i]*my_struct->erroutlier/100);
    }

    double re=0.5;

    double rxmin,rxmax,rymin,rymax;


        str_tab *my_tabx = new str_tab();
        str_tab *my_taby = new str_tab();

        int n_p, p;
        n_p = my_struct->nb_beacon;
        p=2;

        std::vector<bool> v(n_p);
        std::fill(v.begin() + p, v.end(), true);
        unsigned nc = nChoosek(n_p,p);
        int comb[nc*2];
        int cpt=0;
        do {
            for (int i = 0; i < n_p; ++i) {
                if (!v[i]) {
                    comb[cpt]=i+1;
                    cpt++;
                }
            }
        } while (std::next_permutation(v.begin(), v.end()));
        int comb2[nc][2];
        int j=0;
        for(int i=0;i<nc*2;i=i+2){
            comb2[j][0] = comb[i];
            comb2[j][1] = comb[i+1];
            j++;
        }


        my_tabx->lb = new int[2*nc];
        my_tabx->up = new int[2*nc];
        my_tabx->pos1 = new int[2*nc];
        my_tabx->pos2 = new int[2*nc];
        my_tabx->value = new double[2*nc];


        my_taby->lb = new int[2*nc];
        my_taby->up = new int[2*nc];
        my_taby->pos1 = new int[2*nc];
        my_taby->pos2 = new int[2*nc];
        my_taby->value = new double[2*nc];

        str_tab *my_temp_tabx = new str_tab();
        str_tab *my_temp_taby = new str_tab();

        my_temp_tabx->lb = new int[2*nc];
        my_temp_tabx->up = new int[2*nc];
        my_temp_tabx->pos1 = new int[2*nc];
        my_temp_tabx->pos2 = new int[2*nc];
        my_temp_tabx->value = new double[2*nc];


        my_temp_taby->lb = new int[2*nc];
        my_temp_taby->up = new int[2*nc];
        my_temp_taby->pos1 = new int[2*nc];
        my_temp_taby->pos2 = new int[2*nc];
        my_temp_taby->value = new double[2*nc];


        for(int i=0;i<nc;i++){
            //cout<<comb2[i][0]<<" "<<comb2[i][1]<<endl;
            my_struct->comb1 = int(comb2[i][0]);
            my_struct->comb2 = int(comb2[i][1]);

            //cout << my_struct->intervalIn[0]<<";"<<my_struct->intervalIn[1]  <<endl;

            int ci = my_struct->comb1 -1;
            int cj = my_struct->comb2 -1;
            double d_tilde_i = r[ci];
            double d_tilde_j = r[cj];

            double Delta_i = my_struct->beacon_interval*r[ci]/100;
            double Delta_j = my_struct->beacon_interval*r[cj]/100;
            double x0,y0;
            if(my_struct->iteration == 0){
                x0=y0=0;
                my_struct->r_pos_found_prev[0]=0;
                my_struct->r_pos_found_prev[1]=0;
            }
            else{
//                x0=my_struct->robot_position_found[0];
//                y0=my_struct->robot_position_found[1];
                x0=my_struct->r_pos_found_prev[0];
                y0=my_struct->r_pos_found_prev[1];
            }

            double sx_i = x[ci];
            double sy_i = y[ci];

            double sx_j = x[cj];
            double sy_j = y[cj];

            double speedi = my_struct->maxspeed;

            double eps_vdt = speedi;
            double eps0 = my_struct->epsilon_sivia;

            double v_tilde_i = pow(d_tilde_i,2) + pow(Delta_i,2) - pow(x0-sx_i,2) - pow(y0-sy_i,2) - 0.5*pow(eps_vdt+eps0,2);
            double v_tilde_j = pow(d_tilde_j,2) + pow(Delta_j,2) - pow(x0-sx_j,2) - pow(y0-sy_j,2) - 0.5*pow(eps_vdt+eps0,2);

            double delta_i = 2*d_tilde_i*Delta_i + 0.5*pow(eps_vdt+eps0,2);
            double delta_j = 2*d_tilde_i*Delta_j + 0.5*pow(eps_vdt+eps0,2);

            double ai1 = 2*(x0-sx_i);
            double ai2 = 2*(y0-sy_i);
            double aj1 = 2*(x0-sx_j);
            double aj2 = 2*(y0-sy_j);

            double det = 1/(ai1*aj2-ai2*aj1);

            double b1i = det*aj2;
            double b2i = -det*aj1;
            double b1j = -det*ai2;
            double b2j = det*ai1;

            double x1_tilde = b1i*v_tilde_i +b1j*v_tilde_j;
            double x2_tilde = b2i*v_tilde_i +b2j*v_tilde_j;

            double r_1 = fabs(b1i)*delta_i + fabs(b1j)*delta_j;
            double r_2 = fabs(b2i)*delta_i + fabs(b2j)*delta_j;

            double x1_lb = x1_tilde - r_1;
            double x1_ub = x1_tilde + r_1;

            double x2_lb = x2_tilde - r_2;
            double x2_ub = x2_tilde + r_2;

            double lbx,ubx,lby,uby;
            lbx = x1_lb;
            ubx = x1_ub;
            lby = x2_lb;
            uby = x2_ub;
            //cout << lbx<<";"<<ubx<<";"<<lby<<";"<<uby <<endl;

            my_tabx->lb[i] = 1;
            my_tabx->lb[nc+i-1] = 0;
            my_tabx->up[i] = 0;
            my_tabx->up[nc+i-1] = 1;
            my_tabx->pos1[i] = my_struct->comb1;
            my_tabx->pos2[i] = my_struct->comb2;
            my_tabx->pos1[nc+i-1] = my_struct->comb1;
            my_tabx->pos2[nc+i-1] = my_struct->comb2;
            my_tabx->value[nc+i] = ubx;
            my_tabx->value[i] = lbx;

            my_taby->lb[i] = 1;
            my_taby->lb[i+nc-1] = 0;
            my_taby->up[i] = 0;
            my_taby->up[nc+i-1] = 1;
            my_taby->pos1[i] = my_struct->comb1;
            my_taby->pos2[i] = my_struct->comb2;
            my_taby->pos1[nc+i-1] = my_struct->comb1;
            my_taby->pos2[nc+i-1] = my_struct->comb2;
            my_taby->value[nc+i-1] = uby;
            my_taby->value[i] = lby;

            my_temp_tabx->lb[i] = 1;
            my_temp_tabx->lb[nc+i] = 0;
            my_temp_tabx->up[i] = 0;
            my_temp_tabx->up[nc+i] = 1;
            my_temp_tabx->pos1[i] = my_struct->comb1;
            my_temp_tabx->pos2[i] = my_struct->comb2;
            my_temp_tabx->pos1[nc+i] = my_struct->comb1;
            my_temp_tabx->pos2[nc+i] = my_struct->comb2;
            my_temp_tabx->value[nc+i] = ubx;
            my_temp_tabx->value[i] = lbx;

            my_temp_taby->lb[i] = 1;
            my_temp_taby->lb[i+nc] = 0;
            my_temp_taby->up[i] = 0;
            my_temp_taby->up[nc+i] = 1;
            my_temp_taby->pos1[i] = my_struct->comb1;
            my_temp_taby->pos2[i] = my_struct->comb2;
            my_temp_taby->pos1[nc+i] = my_struct->comb1;
            my_temp_taby->pos2[nc+i] = my_struct->comb2;
            my_temp_taby->value[nc+i] = uby;
            my_temp_taby->value[i] = lby;
        }

        std::vector<double> sortedvector (my_tabx->value,my_tabx->value+2*nc);

        // using default comparison (operator <):
        std::sort (sortedvector.begin(), sortedvector.begin()+2*nc);


        for(int i=0;i<2*nc;i++){
            int j=0;
            for (std::vector<double>::iterator it=sortedvector.begin(); it!=sortedvector.end(); ++it){
    //            std::cout << ' ' << *it;
                if(my_temp_tabx->value[i]==*it){
                    my_tabx->lb[j] = my_temp_tabx->lb[i];
                    my_tabx->up[j] = my_temp_tabx->up[i];
                    my_tabx->pos1[j] = my_temp_tabx->pos1[i];
                    my_tabx->pos2[j] = my_temp_tabx->pos2[i];
                    my_tabx->value[j] = my_temp_tabx->value[i];
                }
                j++;
            }
        }
        std::vector<double> sortedvectory (my_taby->value,my_taby->value+2*nc);

        // using default comparison (operator <):
        std::sort (sortedvectory.begin(), sortedvectory.begin()+2*nc);

        for(int i=0;i<2*nc;i++){
            int j=0;
            for (std::vector<double>::iterator it=sortedvectory.begin(); it!=sortedvectory.end(); ++it){
    //            std::cout << ' ' << *it;
                if(my_temp_taby->value[i]==*it){
                    my_taby->lb[j] = my_temp_taby->lb[i];
                    my_taby->up[j] = my_temp_taby->up[i];
                    my_taby->pos1[j] = my_temp_taby->pos1[i];
                    my_taby->pos2[j] = my_temp_taby->pos2[i];
                    my_taby->value[j] = my_temp_taby->value[i];
                }
                j++;
            }
        }
//        std::cout << "myvectorx contains:";
//        for(int i=0;i<2*nc;i++){
//            cout<<my_temp_tabx->value[i]<<";"<<my_temp_tabx->pos1[i]<<my_temp_tabx->pos2[i]<<";"<<my_temp_tabx->lb[i]<<"|";
//        }
//        cout<<"\n";
//        std::cout << "myvectorx2 contains:";
//        for(int i=0;i<2*nc;i++){
//            cout<<my_tabx->value[i]<<";"<<my_tabx->pos1[i]<<my_tabx->pos2[i]<<";"<<my_tabx->lb[i]<<"|";
//        }
//        cout<<"\n";
////        std::cout << "myvectory contains:";
////        for(int i=0;i<2*nc;i++){
////            cout<<my_temp_taby->value[i]<<";"<<my_temp_taby->pos1[i]<<my_temp_taby->pos2[i]<<";"<<my_temp_taby->lb[i]<<"|";
////        }
////        cout<<"\n";
//        std::cout << "myvectory2 contains:";
//        for(int i=0;i<2*nc;i++){
//            cout<<my_taby->value[i]<<";"<<my_taby->pos1[i]<<my_taby->pos2[i]<<";"<<my_taby->lb[i]<<"|";
//        }
//        std::cout << endl<<"Imax:";
        int Imax=0;
        int Imaxy=0;
        int tab_imax[2*nc];
        int tab_imaxy[2*nc];
        for(int i=0;i<2*nc;i++){
            if((my_tabx->up[i])==0){   Imax++;}
            if((my_tabx->up[i])==1){   Imax--;}
            tab_imax[i] = Imax;
            if((my_taby->up[i])==0){   Imaxy++;}
            if((my_taby->up[i])==1){   Imaxy--;}
            tab_imaxy[i] = Imaxy;
            //cout<<Imax<<"-";
        }
        Imax = -20;
        int indmax;
        for(int i=0;i<2*nc;i++){
            if(tab_imax[i]>Imax){
                indmax= i;
                Imax = tab_imax[i];
            }
        }
        //cout<<indmax<<";"<<Imax;
        int indmaxy;
        Imaxy = -20;
        for(int i=0;i<2*nc;i++){
            if(tab_imaxy[i]>Imaxy){
                indmaxy= i;
                Imaxy = tab_imaxy[i];
            }
        }

        double ubx_found=0,uby_found=0;
        for(int i=0;i<2*nc;i++){
            if(my_tabx->pos1[indmax] == my_tabx->pos1[i] && my_tabx->pos2[indmax] == my_tabx->pos2[i] && my_tabx->up[i] == 1)
            {ubx_found = my_tabx->value[i];}
            if(my_taby->pos1[indmaxy] == my_taby->pos1[i] && my_taby->pos2[indmaxy] == my_taby->pos2[i] && my_taby->up[i] == 1)
            {uby_found = my_taby->value[i];}
        }



        rxmin = my_struct->r_pos_found_prev[0] + my_tabx->value[indmax];
        rxmax = my_struct->r_pos_found_prev[0] + ubx_found;
        rymin = my_struct->r_pos_found_prev[1] + my_taby->value[indmaxy];
        rymax = my_struct->r_pos_found_prev[1] + uby_found;

        my_struct->lbx = rxmin;
        my_struct->ubx = rxmax;
        my_struct->lby = rymin;
        my_struct->uby = rymax;

        my_struct->r_pos_found_prev[0] += ( my_tabx->value[indmax] + ubx_found)/2;
        my_struct->r_pos_found_prev[1] += ( my_taby->value[indmaxy] + uby_found )/2;

//        rxmin = my_struct->r_pos_found_prev[0] - (ubx_found + my_tabx->value[indmax])/2;
//        rxmax = my_struct->r_pos_found_prev[0] + (ubx_found + my_tabx->value[indmax])/2;
//        rymin = my_struct->r_pos_found_prev[1] - (uby_found + my_taby->value[indmaxy])/2;
//        rymax = my_struct->r_pos_found_prev[1] + (uby_found + my_taby->value[indmaxy])/2;

//        cout<<"\nMoved of:["<<my_tabx->value[indmax]<<";"<<ubx_found<<"] || ["<<my_taby->value[indmaxy]<<";"<<uby_found<<"]"<<endl<<endl;



    for(int i=0;i<n;i++)
        R.DrawEllipse(x[i],y[i],re,QPen(Qt::black),QBrush(Qt::NoBrush));
    R.DrawEllipse(my_struct->r_pos_found_prev[0],my_struct->r_pos_found_prev[1],0.1,QPen(Qt::red),QBrush(Qt::Dense4Pattern));
    R.DrawBox(rxmin,rxmax,rymin,rymax,QPen(Qt::black),QBrush(Qt::NoBrush));
    sortedvector.clear();
    sortedvectory.clear();
    v.clear();
}

Sivia::Sivia(repere& R,struct sivia_struct *my_struct) : R(R) {
    if(my_struct->pairs==1){
        Sivia_Pair(my_struct);
    }
    else{
        my_struct->areain = 0;
        my_struct->areap = 0;
        my_struct->isinside=0;
        Variable xvar,yvar,zvar,tvar;
        int n = my_struct->nb_beacon;
        double *x=my_struct->x; // vecteur des abcisses des donnees
        double *y=my_struct->y; // vecteur des ordonnees des donnees
        double *z=my_struct->z;
        double *r=new double[n]; // vecteur des rayons
        //r1=9.0;r2=2.0;r3=6.0;r4=6.0;r5=10.0;
        //random config (original config)
        //y1=x1=14;x2=4;y2=-7;x3=7;y3=10;x4=y4=-10;x5=-4;y5=12;x6=0;y6=0;
        //square config
    //    x1=-24;y1=24;x2=-24;y2=-24;x3=24;y3=24;x4=24;y4=-24;x5=y5=0;x6=y6=0;

        //Create a square pattern of beacon with the 5th in the center, after 6 beacons, their position is generated at random
    //    if (n>=2){
    //        x[0]=-24;y[0]=24;x[1]=-24;y[1]=22;
    //    }
    //    if (n>=3){
    //        x[2]=-24;y[2]=-24;
    //    }
    //    if (n>=4){
    //        x[3]=24;y[3]=24;
    //    }
    //    if (n>=5){
    //        x[4]=24;y[4]=-24;
    //    }
    //    if (n>=6){
    //        x[5]=y[5]=0;
    //    }

        double xr=my_struct->robot_position[0],yr=my_struct->robot_position[1],zr=my_struct->robot_position[2];

        for (int i=0;i<n;i++) {
    //        r[i]= sqrt(pow(xr-x[i],2)+pow(yr-y[i],2)+pow(zr-z[i],2));
            r[i]= sqrt(pow(xr-x[i],2)+pow(yr-y[i],2));
            if (my_struct->outliers[i]!=0)
                r[i] *= (1+my_struct->outliers[i]*my_struct->erroutlier/100);
        }

        vector<Function*> f;

        for(int i=0;i<n;i++) {
    //        f.push_back(new Function(xvar,yvar,zvar,sqrt(sqr(xvar-x[i])+sqr(yvar-y[i])+sqr(zvar-z[i]))));
            f.push_back(new Function(xvar,yvar,sqrt(sqr(xvar-Interval(x[i]-my_struct->beacon_interval*r[i]/100,x[i]+my_struct->beacon_interval*r[i]/100))
                                                    +sqr(yvar-Interval(y[i]-my_struct->beacon_interval*r[i]/100,y[i]+my_struct->beacon_interval*r[i]/100)))));
        }

        vector<Ctc*> vec_out;
        vector<Ctc*> vec_in;
        CtcNotIn* intemp1,*intemp2,*intemp3;
        CtcIn* outtemp1,*outtemp2,*outtemp3;
    //    for(int i=0;i<n;i++) {

    //        if (cos(th1[i])>0){
    //            outtemp1 = (new CtcIn(*(theta_1[i]),Interval(0,100)));
    //            intemp1 = (new CtcNotIn(*(theta_1[i]),Interval(0,100)));
    //        }
    //        else{
    //            outtemp1 = (new CtcIn(*(theta_1[i]),Interval(-100,0)));
    //            intemp1 = (new CtcNotIn(*(theta_1[i]),Interval(-100,0)));
    //        }
    //        if (cos(th2[i])>0){
    //            outtemp2 = (new CtcIn(*(theta_2[i]),Interval(-100,0)));
    //            intemp2 =(new CtcNotIn(*(theta_2[i]),Interval(-100,0)));
    //        }
    //        else{
    //            outtemp2 = (new CtcIn(*(theta_2[i]),Interval(0,100)));
    //            intemp2 =(new CtcNotIn(*(theta_2[i]),Interval(0,100)));
    //        }
    //        outtemp3 = (new CtcIn(*(f[i]),(r[i]+Interval(-1,1)*my_struct->err[i])));
    //        intemp3 = (new CtcNotIn(*(f[i]),(r[i]+Interval(-1,1)*my_struct->err[i])));

    //        vec_out.push_back(new CtcUnion(*outtemp1,*outtemp2,*outtemp3));
    //        vec_in.push_back(new CtcCompo(*intemp1,*intemp2,*intemp3));
    //    }
    //    free(intemp1);free(intemp2);free(intemp3);free(outtemp1);free(outtemp2);free(outtemp3);
    //    cout<<"coucou"<<endl;
        for(int i=0;i<n;i++) {
            vec_out.push_back(new CtcIn(*(f[i]),(r[i]+Interval(-1,1)*my_struct->err[i])));
            vec_in.push_back(new CtcNotIn(*(f[i]),(r[i]+Interval(-1,1)*my_struct->err[i])));
        }

        double re=0.5;
        int maxq = my_struct->nb_beacon; //nb of contractors
        int ctcq = maxq - my_struct->q + 1; //nb for q-relaxed function of Ibex
        int ninbox = 0;

        CtcQInter insidetmp(vec_in,ctcq);
        CtcQInter outsidetmp(vec_out,my_struct->q);
        CtcFixPoint inside(insidetmp);
        CtcFixPoint outside(outsidetmp);
        IntervalVector box1 = my_struct->box.back();
    //    IntervalVector box(3);box[0]=box[1]=Interval(-25,25);box[2]=Interval(0,2);
        IntervalVector viinside(3);
        //vin.resize(4);

        LargestFirst lf;
        stack<IntervalVector> s;
        s.push(box1);
        my_struct->in_perhaps=0;
        while (!s.empty()) {
            IntervalVector box=s.top();
            s.pop();
            contract_and_draw(inside,box,viinside,1,my_struct,ninbox,Qt::magenta,Qt::red);
            if (box.is_empty()) { continue; }

            contract_and_draw(outside,box,viinside,0,my_struct,ninbox,Qt::darkBlue,Qt::cyan);
            if (box.is_empty()) { continue; }

            if (box.max_diam()<my_struct->epsilon_sivia) {
                R.DrawBox(box[0].lb(),box[0].ub(),box[1].lb(),box[1].ub(),QPen(Qt::yellow),QBrush(Qt::white));
                my_struct->areap += box[0].diam()*box[1].diam();
                my_struct->in_perhaps=1;
            } else {
                pair<IntervalVector,IntervalVector> boxes=lf.bisect(box);
                s.push(boxes.first);
                s.push(boxes.second);
            }
        }

        double tx[ninbox],ty[ninbox],tz[ninbox];

        //cout<<"next"<<ninbox<<endl;
        for(int i=0;i<ninbox;i++){
            IntervalVector cur = (my_struct->vin.back());
            my_struct->vin_prev.push_back(cur);
            //cout<<cur<<endl;
            Interval xcur=cur[0];
            Interval ycur=cur[1];
            tx[i]=xcur.mid();
            ty[i]=ycur.mid();
            my_struct->vin.pop_back();
        }

        for(int i=0;i<ninbox;i++){
            IntervalVector cur = (my_struct->vper.back());
            my_struct->vin_prev.push_back(cur);
            my_struct->vper.pop_back();
        }

        double xin=0,yin=0,zin=0;
        for(int i=0;i<ninbox;i++){
            xin += tx[i];
            yin += ty[i];
            zin += tz[i];
        }
        xin/=double(ninbox);
        yin/=double(ninbox);
        zin/=double(ninbox);

        my_struct->robot_position_found[0] = xin;
        my_struct->robot_position_found[1] = yin;
        my_struct->robot_position_found[2] = zin;
        for(int i=0;i<n;i++)
            R.DrawEllipse(x[i],y[i],re,QPen(Qt::black),QBrush(Qt::NoBrush));

        my_struct->vin.clear();
        my_struct->vper.clear();
        vec_out.clear();
        vec_in.clear();
        f.clear();
    }

}

