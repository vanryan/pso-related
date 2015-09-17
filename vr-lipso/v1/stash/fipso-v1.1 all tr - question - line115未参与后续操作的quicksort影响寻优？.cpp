#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

const int dim = 30;
const int N = 30;
const int testnum = 5;
const int gennum = 1000;
const int coornum = (dim + 1) * N;

double pbest[coornum],update[coornum];

float coor[coornum];
FILE* pso;

void Fun1(int a)                              //Sphere
{
    coor[(dim+1)*a+dim] = 0;
    int i;
    for(i = 0; i < dim; i++)
        coor[(dim+1)*a+dim] += pow(coor[(dim+1)*a+i], 2);
}

void Fun2(int a)                             //Rosenbrock
{
    coor[(dim+1)*a+dim] = 0;
    int i;
    for(i = 0; i < dim; i++)
    {
        if(i < dim - 1)
            coor[(dim+1)*a+dim] += 100 * pow(( coor[(dim+1)*a+i+1]  -  pow(coor[(dim+1)*a+i] , 2) ) , 2)  + pow(coor[(dim+1)*a+i] - 1 , 2);
        else
            coor[(dim+1)*a+dim] += 0;//100*pow(( coor[(dim+1)*a]  -  pow(coor[(dim+1)*a+i] ,2) ) ,2)  + pow(coor[(dim+1)*a+i] - 1 ,2);
    }

}

void Fun3(int a)                             //Rastrigrin
{
    coor[(dim+1)*a+dim] = 0;
    int i;
    for(i = 0; i < dim; i++)
        coor[(dim+1)*a+dim] += pow(coor[(dim+1)*a+i], 2) - 10 * cos(2 * 3.14 * coor[(dim+1)*a+i]) + 10;
}

void Fun4(int a)                           //Griewank
{
    coor[(dim+1)*a+dim] = 0;
    int i;
    float b = 0, c = 1;
    for(i = 0; i < dim; i++)
    {
        b +=  pow(coor[(dim+1)*a+i]  , 2) / 4000;
        c  *=  cos(coor[(dim+1)*a+i]  / sqrt(i + 1));
    }
    coor[(dim+1)*a+dim] = 1 + b - c;
}

int partion(double a[2][N],int p,int r){
    //rand
    srand((unsigned)time(NULL));
    int e=rand()%(r-p+1)+p;
    double tem;
    int tem2;
    tem=a[0][e];
    tem2=a[1][e];
    a[0][e]=a[0][r];
    a[1][e]=a[1][r];
    a[0][r]=tem;
    a[1][r]=tem2;
    double x=a[0][r];
    int i=p-1;
    for (int j=p;j<r;j++){
        if (a[0][j]<=x){
            tem=a[0][i+1];
            a[0][i+1]=a[0][j];
            a[0][j]=tem;
            tem2=a[1][i+1];
            a[1][i+1]=a[1][j];
            a[1][j]=tem2;
            i++;
        }
    }
    tem=a[0][r];
    a[0][r]=a[0][i+1];
    a[0][i+1]=tem;
    tem2=a[1][r];
    a[1][r]=a[1][i+1];
    a[1][i+1]=tem2;
    return i+1;
}

void QuickSort(double a[2][N],int p,int r){
    if (p<r){
        int q = partion(a,p,r);
        QuickSort(a,p,q-1);
        QuickSort(a,q+1,r);
    }
}

void get_vel_update(int lneighs,int gbp)
{
    int i,j,k;
    double coef[N];
    double phi = 4.1;

    double PBNINDEX[2][N]; // Pbest & Index
    for(k=0; k<N; k++)
    {
        PBNINDEX[0][k] = pbest[(dim+1)*k+dim];  // pbest
        PBNINDEX[1][k] = k;  // The index of particle
    }
    QuickSort(PBNINDEX,0,N-1);  // Done sorting the fitness while keeping the index

    // The quicksort gives a small to large sequence! (And fitness is better when smaller
    //for(k=0; k<N; k++)
    //printf("%f,[%f],",PBNINDEX[0][k],PBNINDEX[1][k]);

    PBNINDEX[1][0]=gbp;

    for(j = 0; j < N; j++)
    {
        for(i = 0; i < dim; i++)
        {
            double sum_coef = 0;
            double pos_ave = 0;

            for( k = 0; k < lneighs; k++) // One particle can only learn from a certain number of neighbors
            {
                coef[(int)PBNINDEX[1][k]] = (double)rand() / RAND_MAX;
                sum_coef += coef[(int)PBNINDEX[1][k]];
                pos_ave += coef[(int)PBNINDEX[1][k]] * pbest[(dim+1)*((int)PBNINDEX[1][k])+i];
            }
            // The particle should also learn from itself
            double temp = (double)rand() / RAND_MAX;
            sum_coef += temp;
            pos_ave += temp * pbest[(dim+1)*j+i];

            // After learning
            pos_ave /= sum_coef;
            double var_phi = sum_coef * phi / (lneighs+1);

            double social_central_tendency = pos_ave - coor[(dim+1)*j+i];
            update[(dim+1)*j+i] = var_phi * social_central_tendency ;
        }
    }
}



int main()
{
    //printf("a");
    srand((unsigned)time(NULL));
    int i, j, k, l, vfun, vtest, vgen;
    float x, edge, Vmax, goal, gb, res = 0,avggennew=0;
    float vel[coornum], gbest[coornum], results[2*testnum];


    pso = fopen("pso1-4.txt", "a+");
    fprintf(pso, "num\tw\tc1\tc2\tfuncnum\tavggen\tratio\tfinalresult\n");

    x = 0.7298;

    for(vfun = 1; vfun < 5; vfun++)
    {
        if(vfun == 1)
        {
            edge = 100;
            Vmax = edge;
            goal = 0.01;
        }
        if(vfun == 2)
        {
            edge = 30;
            Vmax = edge;
            goal = 100;
        }
        if(vfun == 3)
        {
            edge = 5.12;
            Vmax = edge;
            goal = 100;
        }
        if(vfun == 4)
        {
            edge = 600;
            Vmax = edge;
            goal = 0.1;
        }
        int ngoal = 0, avggen = 0;

        for(i = 0; i < testnum; i++)
        {
            results[2*i] = gennum;
            results[2*i+1] = 0;
        }

        for(vtest = 0; vtest < testnum; vtest++)
        {
            int sum = 1;
            for(i = 0; i < coornum; i++)
            {
                pbest[i] = INT_MAX;
                gbest[i] = INT_MAX;
            }

            for(i = 0; i < dim + 1; i++)               //粒子初始化；
            {
                for(j = 0; j < N; j++)
                {
                    if(i < dim)
                    {
                        coor[(dim+1)*j+i] = (double)rand() / RAND_MAX *  2 * edge - edge;
                        vel[(dim+1)*j+i] = (double)rand() / RAND_MAX * 2 * Vmax - Vmax;
                        pbest[(dim+1)*j+i] = coor[(dim+1)*j+i];
                        //printf("%f\t", coor[(dim+1)*j+i]);
                    }
                    else
                    {
                        // i=dim
                        if(vfun == 1)
                            Fun1(j);
                        if(vfun == 2)
                            Fun2(j);
                        if(vfun == 3)
                            Fun3(j);
                        if(vfun == 4)
                            Fun4(j);
                        //printf("\n%d:%f,", j, coor[(dim+1)*j+i]);
                        pbest[(dim+1)*j+i] = coor[(dim+1)*j+i];
                        if(coor[(dim+1)*j+i] < gbest[(dim+1)*j+i])
                        {
                            for(l = 0; l < dim + 1; l++)
                                gbest[(dim+1)*j+l] = coor[(dim+1)*j+l];
                            //printf("gbest:%f\t", gbest[(dim+1)*j+i]);
                        }
                        for(k = 0; k < N; k++)
                            if( gbest[(dim+1)*j+i] < gbest[(dim+1)*k+i])
                            {
                                //printf("t[%d][%d],",j,k);
                                for(l = 0; l < dim + 1; l++)
                                    gbest[(dim+1)*k+l] = gbest[(dim+1)*j+l];
                            }
                    }
                }
            }

            int gbp;
            gb = pbest[dim];
            for(i = 1; i < N; i++)
                if(pbest[(dim+1)*i+dim] < gb)
                {
                    gb = pbest[(dim+1)*i+dim];
                    gbp = i;
                }

            for(vgen = 0; vgen < gennum; vgen++)            //gennum代进化；
            {
                /*
                //px
                int iteqs;
                double PBNINDEX[2][N]; // Pbest & Index
                for(iteqs=0; iteqs<N; iteqs++)
                {
                    PBNINDEX[0][iteqs]=pbest[(dim+1)*iteqs+dim];  // pbest
                    PBNINDEX[1][iteqs]=iteqs;  // The index of particle
                }
                //QuickSort(PBNINDEX,0,N-1);
                //px
*/
                //gbp = (int)PBNINDEX[1][0];
                //printf("gbp:<%d,%f>\n",gbp,pbest[(dim+1)*gbp+dim]);

                get_vel_update(1,gbp);

                double phi = 4.1;
                for(i = 0; i < dim + 1; i++)
                {
                    for(j = 0; j < N; j++)
                    {
                        if(i < dim)
                        {
                            //vel[(dim+1)*j+i] = x *(vel[(dim+1)*j+i] + update[(dim+1)*j+i]);
                            float r1 = (double)rand() / RAND_MAX * 1;
                            float r2 = (double)rand() / RAND_MAX * 1;
                            //vel[(dim+1)*j+i] = x * vel[(dim+1)*j+i] + c1 * r1 * (pbest[(dim+1)*j+i] - coor[(dim+1)*j+i])\
                            + c2 * r2 * (gbest[(dim+1)*j+i] - coor[(dim+1)*j+i]);
                            //vel[(dim+1)*j+i] = x * ( vel[(dim+1)*j+i] + phi*( r1*( pbest[(dim+1)*gbp+i]-coor[(dim+1)*j+i])\
                                                     + r2*(pbest[(dim+1)*j+i]-coor[(dim+1)*j+i]) )/2 );
                              vel[(dim+1)*j+i] = x * (vel[(dim+1)*j+i]+update[(dim+1)*j+i]);
                            /*if(vel[(dim+1)*j+i] > Vmax)
                                vel[(dim+1)*j+i] = Vmax;
                            if(vel[(dim+1)*j+i] < -Vmax)
                                vel[(dim+1)*j+i] = -Vmax;*/

                            coor[(dim+1)*j+i] += vel[(dim+1)*j+i];
                            /*if(coor[(dim+1)*j+i] > edge)
                            {
                                coor[(dim+1)*j+i] = edge;
                                vel[(dim+1)*j+i] = 0;
                            }

                            if(coor[(dim+1)*j+i] < -edge)
                            {
                                coor[(dim+1)*j+i] = -edge;
                                vel[(dim+1)*j+i] = 0;
                            }*/
                        }
                        else
                        {
                            if(vfun == 1)
                                Fun1(j);
                            if(vfun == 2)
                                Fun2(j);
                            if(vfun == 3)
                                Fun3(j);
                            if(vfun == 4)
                                Fun4(j);
                            if(coor[(dim+1)*j+i] < pbest[(dim+1)*j+i])
                                for(l = 0; l < dim + 1; l++)
                                    pbest[(dim+1)*j+l] = coor[(dim+1)*j+l];
                            if(coor[(dim+1)*j+i] < gbest[(dim+1)*j+i])
                            {
                                for(l = 0; l < dim + 1; l++)
                                    gbest[(dim+1)*j+l] = coor[(dim+1)*j+l];
                            }
                            for(k = 0; k < N; k++)
                                if( gbest[(dim+1)*j+i] < gbest[(dim+1)*k+i])
                                    for(l = 0; l < dim + 1; l++)
                                        gbest[(dim+1)*k+l] = gbest[(dim+1)*j+l];
                        }
                    }
                }
                gb = pbest[dim];
                for(i = 1; i < N; i++)
                    if(pbest[(dim+1)*i+dim] < gb)
                    {
                        gb = pbest[(dim+1)*i+dim];
                        gbp = i;
                    }

                if(gb < goal && sum != 0)
                {
                    results[2*vtest] = (double)vgen;
                    results[2*vtest+1] = gb;
                    sum = 0;
                    //ngoal++;
                    break;
                }
                if(vgen == gennum - 1)
                {
                    results[2*vtest] = gennum - 1;
                    results[2*vtest+1] = gb;
                }

            } // End of vgen
        } // End of vtest

        for(i = 0; i < testnum; i++)
        {
            if(results[2*i] != gennum - 1)
            {
                ngoal++;
                avggen += results[2*i];
                res += results[2*i+1];
            }
        }

        if(ngoal == 0){
            avggennew = (double)gennum - 1;
            printf("\tFail.");
        }
        else{
            avggennew = (double)avggen / ngoal;
            printf("\tAve Ites:%d,Suc Tests%d\t", avggen/ngoal, ngoal);
        }

        fprintf(pso, "\n平均代数为：%f\n达优率为：%f%%", avggennew, (double)ngoal / testnum * 100);
    } // end of vfun
    fclose(pso);

    return 0;
}
