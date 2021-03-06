/*
 Comprehensive Learning Particle Swarm Optimization

 By Van Ryan

 @Constant Inertia Weight

 @No restraints for velocities or positions

 -For A Specific Function
*/

#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

#include"../settings.h"

int main()
{
    srand((unsigned)time(NULL));

    /*
        Generic & Procedural Vars
    */
    int i, j, k, l, vfun, vtest, vgen;
    double edge, Vmax, goal, gb;

    /*
        Vars for Real-time Recording
    */
    int sucCounter;

    /*
        Vars for Storage Section
     */
    double AveOptGen; // Average number of generations of success optimization
    double AveOptRate; // Ratio of success tests/ total tests
    double AveEndFit; // End fit
    double AveSucFit;// SUCCESSFUL EndFit The average of the gbest fitness in the end of each SUCCESSFUL test (last generation)
    double FIndex, recordInd[testnum];
    /* End of Storage Section */

    /*
        Initialization of some variables
    */

    AveOptGen=0;
    AveOptRate=0;
    AveSucFit=0;
    AveEndFit=0;
    FIndex=99; // Just a dummy and impossible value

    pso = fopen("clpso.txt", "a+");  // We cannot use 'w' on the server in the lab in lack of privileges.

    vfun=4;

    switch(vfun)
    {
    case 1:
        edge = 100;
        goal = 0.01;
        break;
    case 2:
        edge = 30;
        goal = 100;
        break;
    case 3:
        edge = 5.12;
        goal = 100;
        break;
    case 4:
        edge = 600;
        goal = 0.1;
        break;
    default:
        edge = 100;
        goal = 0.01;
    }

    Vmax = edge;
    //Vmax = edge/2;

    //int ngoal = 0, avggen = 0;
    sucCounter = 0;

    for(vtest = 0; vtest < testnum; vtest++)
    {
        int sum = 1;
        for(i = 0; i < coornum; i++)
        {
            pbest[i] = INT_MAX;
            gbest[i] = INT_MAX;
        }

        for(i = 0; i < dim + 1; i++)               //Initialization
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

        // refgap_flag: A flag for counting the refreshing gap
        // ite_a: A iteration counting variable
        // cl_rec: A 2-dim recording array for comprehensive learning
        int refgap_flag = 0, ite_a, cl_rec[N];

        for(vgen = 0; vgen < gennum; vgen++)            //gennum��������
        {
            if(refgap_flag == refgap )
                refgap_flag = 0;

            float inert_w;
            inert_w = (wmax-wmin)*(gennum-vgen-1)/gennum +wmin;

            for(i = 0; i < dim + 1; i++)
            {
                if(refgap_flag == 0 )
                {
                    for(ite_a=0; ite_a<N; ite_a++)
                    {
                        float rand_p = (float) rand() / RAND_MAX * 1;
                        if(rand_p<=Prob_c)
                        {
                            // Tournament
                            int learnee1=ceil(rand() / RAND_MAX * N) - 1;  // Get a integer number from 0 to N-1
                            int learnee2=ceil(rand() / RAND_MAX * N) - 1;  // Get a integer number from 0 to N-1
                            if(pbest[(dim+1)*learnee1+dim]< pbest[(dim+1)*learnee2+dim])
                                cl_rec[ite_a] = learnee1;
                            else
                                cl_rec[ite_a] = learnee2;
                            // End of tournament
                        }

                        else
                            cl_rec[ite_a] = ite_a;
                    } // End of for(ite_a=0; ite_a<N; ite_a++; )
                } // End of if(refgap_flag == 0 )

                for(j = 0; j < N; j++)
                {
                    if(i < dim)
                    {
                        double rand_1 = (double)rand() / RAND_MAX * 1;
                        //double r2 = (double)rand() / RAND_MAX * 1;

                        vel[(dim+1)*j+i] = inert_w * vel[(dim+1)*j+i] + para_c * rand_1 * (pbest[(dim+1)*cl_rec[j]+i] - coor[(dim+1)*j+i]);
                        //vel[(dim+1)*j+i] = getMin( getMax(vel[(dim+1)*j+i],(-Vmax) ),Vmax);

                        coor[(dim+1)*j+i] += vel[(dim+1)*j+i]; // Moving
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
                    } // End of if/else (i < dim)
                } // End of for(j = 0; j < N; j++)
            } // End of (i = 0; i < dim + 1; i++)

            refgap_flag++;

            gb = gbest[dim];
            for(i = 1; i < N; i++)
                if(gbest[(dim+1)*i+dim] < gb)
                    gb = gbest[(dim+1)*i+dim];

            if(gb < goal && sum != 0)
            {
                sum = 0; // Change the flag
                //ngoal++;
                sucCounter += 1;
                AveOptGen += vgen;
                printf("suc gb:%1.30f\n",gb);
            }

            if(vgen == 1000)
                recordInd[vtest]=gb;

            if(vgen == gennum - 1)
            {
                printf("last gb:%1.30f\n",gb);

                AveEndFit += gb;

                if(sum == 0) // Succeeded
                    AveSucFit += gb;
            }

        } // End of for VGEN!
    } // End of vtest

    FIndex=countFIndex(recordInd);

    AveEndFit = AveEndFit/testnum;

    if(sucCounter!=0) // If succeeded
    {
        AveSucFit = AveSucFit/sucCounter;
        AveOptGen = AveOptGen/sucCounter;
    }
    else
    {
        AveSucFit = -1; // A flag to show the failure
        AveOptGen = -1; // A flag to show the failure
    }


    AveOptRate=(double) sucCounter / testnum * 100;

    printf("%1.30f,%1.30f,%f,%f,%f\n",AveEndFit,AveSucFit,AveOptGen,AveOptRate,FIndex);
    fprintf(pso,"%1.30f,%1.30f,%f,%f,%f\n",AveEndFit,AveSucFit,AveOptGen,AveOptRate,FIndex);
    // AveSucFit is the Average Endfit for only successful runs

    fclose(pso);

    return 0;
}
