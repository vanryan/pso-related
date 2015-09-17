/*
 Fitness-Distance-Ratio Particle Swarm Intelligence

 By Van Ryan

 @Constant Inertia Weight

 @No restraints for velocities or positions

 @Initialization: Vmax=edge(Xrange)

 -For A Specific Function
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

    /*
        Vars for Radius
    */
    float absRadius; // Absolute perceptive radius
    float relRadius; // A iteration var for radius

    /*
        Vars for Real-time Recording
    */
    int adjmat[N][N]; // Adjacent Matrix
    int sucCounter;
    int radRec; // A var to count the index of the radius

    /*
        Vars for Storage Section
     */
    double AveOptGen[relRange]; // Average number of generations of success optimization
    double AveOptRate[relRange]; // Ratio of success tests/ total tests
    double AveEndFit[relRange]; // End fit
    double AveSucFit[relRange];// SUCCESSFUL EndFit The average of the gbest fitness in the end of each SUCCESSFUL test (last generation)
    double FIndex[relRange], recordInd[testnum];
    /* End of Storage Section */

    /*
        Initialization of some variables
    */
    for(i=0; i<relRange; i++)
    {
        AveOptGen[i]=0;
        AveOptRate[i]=0;
        AveSucFit[i]=0;
        AveEndFit[i]=0;
        FIndex[i]=99; // Just a dummy and impossible value
    }


    pso = fopen("../../vicfdrpso.txt", "a+");  // We cannot use 'w' on the server in the lab in lack of privileges.

    vfun=4;

    /*
     Initiate the hybrid composition function
    */
    hyb_func4_init(fmax,or_matrix,opt_arr,opt_arr_old);

    edge = hyb_func_goalnedgenvmax[vfun-1][1];
    Vmax = hyb_func_goalnedgenvmax[vfun-1][2];
    goal = hyb_func_goalnedgenvmax[vfun-1][0];

    radRec=0;

    for(relRadius=relStart; relRadius<relEnd; relRadius+=relStep) // Loop layer 1 (relRadius)
    {
        sucCounter = 0;

        for(vtest = 0; vtest < testnum; vtest++)
        {
            int sum = 1;
            for(i = 0; i < coornum; i++)
            {
                pbest[i] = INT_MAX;
                gbest[i] = INT_MAX;
            }

            for(i = 0; i < dim + 1; i++)         //Initialization
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

            // Updating each bird's gbest:
            for(j=0; j<N; j++)
            {
                for(k=0; k<N; k++)
                    if( adjmat[j][k] != 0 ) // k Can be perceived by j
                    {
                        // Updating the gbest of j
                        //printf("%f\t",getDistance(j,k));
                        if( coor[(dim+1)*k+dim] < gbest[(dim+1)*j+dim] )
                        {
                            for(l=0; l<=dim; l++)
                                gbest[(dim+1)*j+l] = coor[(dim+1)*k+l];
                        }
                    }
                //printf("gbest of %d:%f\nfitness of %d:%f\n",j,gbest[(dim+1)*j+dim],j,coor[(dim+1)*j+dim]);
            }

            gb = gbest[dim]; // Updating the global gbest
            for(i = 1; i < N; i++)
                if(gbest[(dim+1)*i+dim] < gb)
                    gb = gbest[(dim+1)*i+dim];

            //End of initialization

            for(vgen = 0; vgen < gennum; vgen++)            //gennum´ú½ø»¯£»
            {
                float inert_w;
                inert_w=(wmax-wmin)*(gennum-vgen-1)/(gennum+wmin);

                for(i = 0; i < dim + 1; i++)
                {
                    for(j = 0; j < N; j++)
                    {
                        if(i < dim)
                        {
                            int lneigh; // The neighbor for the particle to learn from
                            lneigh=findNeigh_fdr(j,i);

                            double r1 = (double)rand() / RAND_MAX * 1;
                            double r2 = (double)rand() / RAND_MAX * 1;
                            double r3 = (double)rand() / RAND_MAX * 1;
                            vel[(dim+1)*j+i] = inert_w * vel[(dim+1)*j+i] + c1 * r1 * (pbest[(dim+1)*j+i] - coor[(dim+1)*j+i])\
                                               + c2 * r2 * (gbest[(dim+1)*j+i] - coor[(dim+1)*j+i])+ c3 * r3 * (pbest[(dim+1)*lneigh+i] - coor[(dim+1)*j+i]);
                            //vel[(dim+1)*j+i] = getMin( getMax(vel[(dim+1)*j+i],(-Vmax) ),Vmax);

                            coor[(dim+1)*j+i] += vel[(dim+1)*j+i]; // Moving
                            //coor[(dim+1)*j+i]  = getMin( getMax(coor[(dim+1)*j+i], (-edge) ) ,edge);
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

                            // When i = dim, updating pbest:
                            if(coor[(dim+1)*j+i] < pbest[(dim+1)*j+i])
                                for(l = 0; l < dim + 1; l++)
                                    pbest[(dim+1)*j+l] = coor[(dim+1)*j+l];

                            // Comparing itself with gbest:
                            if(coor[(dim+1)*j+i] < gbest[(dim+1)*j+i])
                            {
                                for(l = 0; l < dim + 1; l++)
                                    gbest[(dim+1)*j+l] = coor[(dim+1)*j+l];
                            }

                        } // End of else
                    } // for birds
                } // for dims

                // Getting absolute radius
                 absRadius = relRadius * sqrt(dim * pow(2 * edge , 2));

                // Producing the initiating Adjacent Matrix
                for( j=0; j<N ; j++ )
                    for( k=0; k<N; k++) // Avoid Relapse
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

                // Updating each bird's gbest:
                for( j=0; j<N ; j++ )
                {
                    for( k=0; k<N ; k++ )
                        if( adjmat[j][k] != 0 ) // k Can be perceived by j
                            // Updating the gbest of j
                            if( coor[(dim+1)*k+dim] < gbest[(dim+1)*j+dim] )
                                for(l=0; l<=dim; l++)
                                    gbest[(dim+1)*j+l] = coor[(dim+1)*k+l];
                    //printf("gbest of %d:%f\nfitness of %d:%f\n",j,gbest[(dim+1)*j+dim],j,coor[(dim+1)*j+dim]);
                }
                // End of Updating each bird's gbest

                gb = gbest[dim]; // Updating the global gbest
                for(i = 1; i < N; i++)
                    if(gbest[(dim+1)*i+dim] < gb)
                        gb = gbest[(dim+1)*i+dim];

                if(gb < goal && sum != 0)
                {
                    sum = 0; // Change the flag
                    //ngoal++;
                    sucCounter += 1;
                    AveOptGen[radRec] += vgen;
                    printf("suc gb:%1.15f\n",gb);
                }

                if(vgen+1 == 1000)
                    recordInd[vtest]=gb;

                if(vgen == gennum - 1)
                {
                    printf("last gb:%1.15f\n",gb);

                    AveEndFit[radRec] += gb;

                    if(sum == 0) // Succeeded
                        AveSucFit[radRec] += gb;
                }

            } // End of for VGEN!
        } // End of vtest

        FIndex[radRec]=countFIndex(recordInd);

        AveEndFit[radRec] = AveEndFit[radRec]/testnum;

        if(sucCounter!=0) // If succeeded
        {
            AveSucFit[radRec] = AveSucFit[radRec]/sucCounter;
            AveOptGen[radRec] = AveOptGen[radRec]/sucCounter;
        }
        else
        {
            AveSucFit[radRec] = -1; // A flag to show the failure
            AveOptGen[radRec] = -1; // A flag to show the failure
        }


        AveOptRate[radRec]=(double) sucCounter / testnum * 100;

        printf("%1.15f,%1.15f,%f,%f,%f\n",AveEndFit[radRec],AveSucFit[radRec],AveOptGen[radRec],AveOptRate[radRec],FIndex[radRec]);
        fprintf(pso,"%1.15f,%1.15f,%f,%f,%f\n",AveEndFit[radRec],AveSucFit[radRec],AveOptGen[radRec],AveOptRate[radRec],FIndex[radRec]);
        // AveSucFit is the Average Endfit for only successful runs

        radRec++;

        printf("\n\nRadius %d is tested\n\n",radRec);

    } // End of Loop layer 1 (relRadius)

    printf("\n\n Done! You got it, dude!");

    fclose(pso);

    return 0;
}
