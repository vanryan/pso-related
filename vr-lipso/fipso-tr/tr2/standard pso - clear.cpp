#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

const int dim = 30;
const int N = 30;
const int testnum = 10;
const int gennum = 5000;
const int coornum = (dim + 1) * N;

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



int main()
{
    //printf("a");
    srand((unsigned)time(NULL));
    int i, j, k, l, vfun, vtest, vgen;
    float w, c1, c2, edge, Vmax, goal, gb, res = 0,avggennew=0;
    float vel[coornum], pbest[coornum], gbest[coornum], results[2*testnum];


    pso = fopen("pso1-4.txt", "a+");
    fprintf(pso, "num\tw\tc1\tc2\tfuncnum\tavggen\tratio\tfinalresult\n");

    w = 0.6;
    c1 = 1.7;
    c2 = 1.7;

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

            /*for(i = 0; i < coornum; i++)
            {
                printf("\nabc:%f,%f,%f", coor[i], pbest[i], gbest[i]);
                if(i % (dim + 1) == dim)
                    printf("\n");
            }*/
            for(vgen = 0; vgen < gennum; vgen++)            //gennum代进化；
            {
                for(i = 0; i < dim + 1; i++)
                {
                    for(j = 0; j < N; j++)
                    {
                        if(i < dim)
                        {
                            float r1 = (double)rand() / RAND_MAX * 1;
                            float r2 = (double)rand() / RAND_MAX * 1;
                            vel[(dim+1)*j+i] = w * vel[(dim+1)*j+i] + c1 * r1 * (pbest[(dim+1)*j+i] - coor[(dim+1)*j+i])\
                                               + c2 * r2 * (gbest[(dim+1)*j+i] - coor[(dim+1)*j+i]);
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
                gb = gbest[dim];
                for(i = 1; i < N; i++)
                    if(gbest[(dim+1)*i+dim] < gb)
                        gb = gbest[(dim+1)*i+dim];
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

        printf("\tAve Ites:%d,Suc Tests%d\t", avggen/ngoal, ngoal);

        if(ngoal == 0)
            avggennew = (double)gennum - 1;
        else
            avggennew = (double)avggen / ngoal;
        fprintf(pso, "\n平均代数为：%f\n达优率为：%f%%", avggennew, (double)ngoal / testnum * 100);
    } // end of vfun
    fclose(pso);

    return 0;
}
