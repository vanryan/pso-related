/*
 FIPSO with perceptive radius
 Van Ryan

 For A Specific Function

*/
#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

#include"../hybsettings.h"

#include<fstream>
using namespace std;

#include"ortho_rot_matrix.h"
#include"basicfuncs.h"
#include"hybridfuncs.h"

void get_vel_update(double update[coornum], double pbest[coornum], double c1, double c2, int adjmat[N][N])
{
    int i,j,k;
    double coef[N];
    double phi = 4.01;

    for(j = 0; j < N; j++)
    {
        for(i = 0; i < dim; i++)
        {
            int neighbors = N;
            double sum_coef = 0;
            double pos_ave = 0;

            //denominator = N;

            for( k =0; k < N; k++)
            {
                if(adjmat[j][k]==1)
                {
                    coef[k] = (double)rand() / RAND_MAX;
                    sum_coef += coef[k];
                    pos_ave += coef[k] * pbest[(dim+1)*k+i];
                }
                else neighbors--;

            }
            pos_ave /= sum_coef;
            double var_phi = sum_coef * phi / neighbors;

            double social_central_tendency = pos_ave - coor[(dim+1)*j+i];
            update[(dim+1)*j+i] = var_phi * social_central_tendency ;
        }
    }
}


int main()
{
    srand((unsigned)time(NULL));
    /*
    Variables for the composition functions
    */
    double fmax[basic_func_num], or_matrix[dim][dim], opt_arr[basic_func_num][dim], opt_arr_old[basic_func_num][dim];

    /*
        Generic & Procedural Vars
    */

    int i, j, k, l, vfun, vtest, vgen;
    double edge, Vmax, goal, gb;
    double vel[coornum], pbest[coornum], update[coornum];

    /*
        Vars for Radius
    */
    double absRadius; // Absolute perceptive radius
    double relRadius; // A iteration var for radius


    /*
        Vars for Real-time Recording
    */
    int adjmat[N][N]; // Adjacent Matrix
    double results[2*testnum],res = 0,avggennew=0; // Results[] can be used to record gb
    int sucCounter, radRec;

    /*
       Vars for Storage Section
    */
    double AveOptGen[relRange];  // FuncNumber*RadiusNumber: Every 19[radius] with a func
    double AveOptRate[relRange];  // FuncNumber*RadiusNumber: Every 19[radius] with a func
    //So the truth is, we're getting the average(a matrix) of average(a test)
    //double AveOptFit[relRange][gennum/aveDegInterval+1]; // The average of the gbest fitness of the interval generations
    double AveEndFit[relRange]; // End fit
    double AveSucFit[relRange];// SUCCESSFUL EndFit The average of the gbest fitness in the end of each SUCCESSFUL test (last generation)
    double FIndex[relRange], recordInd[testnum];
    /* End of Storage Section */

    for(i=0; i<relRange; i++)
    {
        AveOptGen[i]=0;
        AveOptRate[i]=0;
        AveSucFit[i]=0;
        AveEndFit[i]=0;
        FIndex[i]=99;
    }

    pso = fopen("../../vicfipso.txt", "a+");

    vfun = 4;

    /*
    Initiate the hybrid composition function
    */
    hyb_func4_init(fmax,or_matrix,opt_arr,opt_arr_old);

    edge = hyb_func_goalnedgenvmax[vfun-1][1];
    Vmax = hyb_func_goalnedgenvmax[vfun-1][2];
    goal = hyb_func_goalnedgenvmax[vfun-1][0];

    radRec = 0; // Recording radius

    for(relRadius = relStart; relRadius < relEnd; relRadius += relStep)
    {

        for(i = 0; i < testnum; i++)
        {
            results[2*i] = gennum;
            results[2*i+1] = 0;
        }


        int ngoal = 0, avggen = 0;
        sucCounter = 0;

        for(vtest = 0; vtest < testnum; vtest++)
        {
            double neoDiag=0;
            double    dist=0;
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
                            hyb_func1(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                        if(vfun == 2)
                            hyb_func2(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                        if(vfun == 3)
                            hyb_func3(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                        if(vfun == 4)
                            hyb_func4(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                        if(vfun == 5)
                            hyb_func5(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                        if(vfun == 6)
                            hyb_func6(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                        //printf("\n%d:%f,", j, coor[(dim+1)*j+i]);
                        pbest[(dim+1)*j+i] = coor[(dim+1)*j+i];
                    }
                }
            }

            absRadius = relRadius * sqrt(dim * pow(2 * edge , 2));

            for(j=0; j<N; j++)
                for(k=0; k<N; k++) // Avoid Relapse
                {
                    if( getDistance(j,k) <= absRadius && j!=k ) // k Can be perceived by j
                    {
                        adjmat[j][k]=1;
                    }
                    else
                    {
                        adjmat[j][k]=0;
                    }
                }

            //End of initialization


            for(vgen = 0; vgen < gennum; vgen++)            //gennum代进化；
            {
                get_vel_update(update, pbest, c1, c2, adjmat);

                for(i = 0; i < dim + 1; i++)
                {
                    for(j = 0; j < N; j++)
                    {
                        if(i < dim)
                        {
                            vel[(dim+1)*j+i] = x * ( vel[(dim+1)*j+i] + update[(dim+1)*j+i] );
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
                                hyb_func1(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                            if(vfun == 2)
                                hyb_func2(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                            if(vfun == 3)
                                hyb_func3(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                            if(vfun == 4)
                                hyb_func4(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                            if(vfun == 5)
                                hyb_func5(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                            if(vfun == 6)
                                hyb_func6(coor,j,or_matrix,fmax,opt_arr,opt_arr_old);
                            if(coor[(dim+1)*j+i] < pbest[(dim+1)*j+i])
                                for(l = 0; l < dim + 1; l++)
                                    pbest[(dim+1)*j+l] = coor[(dim+1)*j+l];
                        }
                    }
                }

                if( (vgen+1) % (changeInte) == 0 && vgen!=0)
                {
                    //printf("%d,%f\n",vgen,neoDiag);
                    absRadius=relRadius*neoDiag;
                }

                neoDiag=0;// reset

                for(j=0; j<N; j++)
                    for(k=0; k<N; k++) // Avoid Relapse
                    {
                        dist=getDistance(j,k);
                        if( getDistance(j,k) <= absRadius && j!=k ) // k Can be perceived by j
                        {
                            adjmat[j][k]=1;
                        }
                        else
                        {
                            adjmat[j][k]=0;
                        }
                        if( (vgen+1) % changeInte == (changeInte-1) )
                            if(dist>neoDiag)
                            {
                                neoDiag=dist;
                                //printf("%d,%d,dist:%f\t",j,k,dist);
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
                    //printf("\ngb:%f,%d\t",gb,vgen);
                    sucCounter+=1;
                }

                if(vgen == 1000)
                    recordInd[vtest]=gb;

                if(vgen == gennum - 1)
                {
                    printf("lastgb:%1.30f\n",gb);

                    AveEndFit[radRec] += gb;

                    if(sum == 0) // Succeeded
                        AveSucFit[radRec] += gb;
                    else
                        results[2*vtest] = vgen;
                }

            } // End of for VGEN!

            //printf("%1.30f",gb);
        } // End of for VTEST   !

        FIndex[radRec]=countFIndex(recordInd);

        AveEndFit[radRec] = AveEndFit[radRec]/testnum;

        if(sucCounter!=0)
            AveSucFit[radRec] = AveSucFit[radRec]/sucCounter;

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

        AveOptGen[radRec] = avggennew;
        AveOptRate[radRec] = (double)sucCounter / testnum * 100;
        //fprintf(pso, "\n平均达优代数为：%f\n达优率为：%f%%", avggennew, (double)sucCounter / testnum * 100);


        if(AveOptRate[radRec]==0)  //Failed in this radius
        {
            printf("%1.30f,-1,%f,%f,%f\n",AveEndFit[radRec],AveOptGen[radRec],AveOptRate[radRec],FIndex[radRec]); // Record a radius in a row, literally
            fprintf(pso,"%1.30f,-1,%f,%f,%f\n",AveEndFit[radRec],AveOptGen[radRec],AveOptRate[radRec],FIndex[radRec]); // Record a radius in a row, literally
        }
        else
        {
            printf("%1.30f,%1.30f,%f,%f,%f\n",AveEndFit[radRec],AveSucFit[radRec],AveOptGen[radRec],AveOptRate[radRec],FIndex[radRec]); // Record a radius in a row, literally
            fprintf(pso,"%1.30f,%1.30f,%f,%f,%f\n",AveEndFit[radRec],AveSucFit[radRec],AveOptGen[radRec],AveOptRate[radRec],FIndex[radRec]); // Record a radius in a row, literally
        }

        radRec++;

    } // End of relRadius

    fclose(pso);
}
