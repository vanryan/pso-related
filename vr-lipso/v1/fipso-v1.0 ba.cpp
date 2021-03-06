/*
 FIPSO v1
 Van Ryan

 BA network
 2013.12.10

*/
#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

const int dim = 30;
const int N = 50;
const int testnum = 100;
const int gennum = 5000;
const int coornum = (dim + 1) * N;
int t[N][N]; //Adjacent Matrix

float coor[coornum];
FILE* pso;
FILE* bafp;

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

void Fun5(int a)
{
    coor[(dim+1)*a+dim] = 0;

}

void get_vel_update(float update[coornum], float pbest[coornum], int t[N][N])
{
    int i,j,k;
    double coef[N];
    double phi = 4.01;

    for(j = 0; j < N; j++)
    {
        for(i = 0; i < dim; i++)
        {
            double sum_coef = 0;
            double pos_ave = 0;

            //denominator = N;

            for( k =0; k < N; k++)
            {
                if(k!=j && t[j][k]!=INT_MAX){ // Considering a particle's neighbors
                    coef[k] = (double)rand() / RAND_MAX;
                    sum_coef += coef[k];
                    pos_ave += coef[k] * pbest[(dim+1)*k+i];
                }
            }
            pos_ave /= sum_coef;
            double var_phi = sum_coef * phi / N;

            double social_central_tendency = pos_ave - coor[(dim+1)*j+i];
            update[(dim+1)*j+i] = var_phi * social_central_tendency ;
        }
    }
}

void BAnet_produce()
{
    if((bafp = fopen("fipso-v1-ba-dis.txt", "rt")) == NULL)
    {
        printf("ERROR!");
        return ;
    }
    int i, j;
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
        {
            fscanf(bafp, "%d,", &t[i][j]);
            if(i == 0 && j == 0)
                t[i][j] = 0;
            if(t[i][j] == 2)
                t[i][j] = INT_MAX;

        }
    fclose(bafp);

}


int main()
{
    //printf("a");
    srand((unsigned)time(NULL));
    int i, j, l, vfun, vtest, vgen;
    float x, edge, Vmax, goal, gb, res = 0,avggennew=0;
    float vel[coornum], pbest[coornum], update[coornum], results[2*testnum];
    float AveEndFit[4];
    int sucCounter;
    for(i = 0; i < 4; i++)
        AveEndFit[i]=0;

    BAnet_produce(); // Produce the particle's BA distribution

    pso = fopen("fipso1.0-ba-f1-f4.txt", "a+");
    fprintf(pso, "num\tw\tc1\tc2\tfuncnum\tavggen\tratio\tfinalresult\n");

    x = 0.7298;

    for(vfun = 1; vfun < 5; vfun++)
    {
        for(i = 0; i < testnum; i++)
        {
            results[2*i] = gennum;
            results[2*i+1] = 0;
        }

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
        sucCounter = 0;

        for(vtest = 0; vtest < testnum; vtest++)
        {
            int sum = 1;
            for(i = 0; i < coornum; i++)
            {
                pbest[i] = INT_MAX;
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
                    }
                }
            }


            for(vgen = 0; vgen < gennum; vgen++)            //gennum代进化；
            {
                get_vel_update(update, pbest, t);

                for(i = 0; i < dim + 1; i++)
                {
                    for(j = 0; j < N; j++)
                    {
                        if(i < dim)
                        {
                            vel[(dim+1)*j+i] = x * ( vel[(dim+1)*j+i] + update[(dim+1)*j+i] );

                            coor[(dim+1)*j+i] += vel[(dim+1)*j+i];

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
                        }
                    }
                }
                gb = pbest[dim];
                for(i = 1; i < N; i++)
                    if(pbest[(dim+1)*i+dim] < gb)
                        gb = pbest[(dim+1)*i+dim];
                if(gb < goal && sum != 0)
                {
                    results[2*vtest] = (double)vgen;
                    results[2*vtest+1] = gb;
                    sum = 0;
                    //ngoal++;
                    printf("\ngb:%f,%d\t",gb,vgen);
                    sucCounter+=1;
                }
                if(vgen == gennum - 1)
                {
                    printf("lastgb:%f\n",gb);

                    if(sum == 0) // Succeeded
                        AveEndFit[vfun-1] += gb;
                    else
                        results[2*vtest] = vgen;
                }

            } // End of for VGEN!
        } // End of for VTEST   !

        if(sucCounter!=0)
            fprintf(pso,"\n\n\nLast succeeded average gb fitness:%f\n",AveEndFit[vfun-1]/sucCounter);
        else
            fprintf(pso,"\n\n\nLast succeeded average gb fitness: U gotta be kidding, man!\n");

        for(i = 0; i < testnum; i++)
        {
            if(results[2*i] != gennum - 1)
            {
                ngoal++;
                avggen += results[2*i];
                res += results[2*i+1];
            }
        }
        printf("%d,%d", avggen/vtest, ngoal);
        if(ngoal == 0)
            avggennew = (double)gennum - 1;
        else
            avggennew = (double)avggen / ngoal;
        fprintf(pso, "\n平均代数为：%f\n达优率为：%f%%", avggennew, (double)sucCounter / testnum * 100);
    } // End of for VFUN!

    fclose(pso);
}
