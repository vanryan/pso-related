/*
PSO with  Changing Perceptive Radius
Version 2
2013.11.07

Van Ryan
*/

#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"
#include"string.h"

#include"../hyb-settings.h"

#include<fstream>
using namespace std;

#include"ortho_rot_matrix.h"
#include"basicfuncs.h"
#include"hybridfuncs.h"


FILE* pso;

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
    int i, j, k, l, vfun, vtest, vgen; // iteration vars
    double edge, Vmax, goal, gb;  // gb is the global gbest
    double vel[coornum], pbest[coornum], gbest[coornum]; // Each bird's got their own gbest

    /*
        Vars for Radius
    */
    double absRadius; // Absolute perceptive radius
    double relRadius; // A iteration var for radius

    /*
        Vars for Real-time Recording
    */
    int adjmat[N][N]; // Adjacent Matrix
    int radRec; // A var to count the index of the radius
    double res = 0, avggennew=0, results[4*testnum];

    /*
        Vars for Storage
    */
    double AveOptGen[relRange];  // FuncNumber*RadiusNumber: Every 19[radius] with a func
    double AveOptRate[relRange];  // FuncNumber*RadiusNumber: Every 19[radius] with a func
    //So the truth is, we're getting the average(a matrix) of average(a test)
    //double AveOptFit[relRange][gennum/aveDegInterval+1]; // The average of the gbest fitness of the interval generations
    double AveEndFit[relRange]; // End fit
    double AveOptFit[relRange][gennum/aveDegInterval+1]; // The average of the gbest fitness of the interval generations
    double AveSucFit[relRange];// SUCCESSFUL EndFit The average of the gbest fitness in the end of each SUCCESSFUL test (last generation)
    double FIndex[relRange], recordInd[testnum];
    /* End of Storage Section */


    char buffer[19]; // A buffer to write the file, 15 is the estimated char numbers
    char storage[19*relRange];


    for(i=0; i<relRange; i++)
    {
        AveOptGen[i]=0;
        AveOptRate[i]=0;
        AveEndFit[i]=0;
        AveSucFit[i]=0;
        for(j=0; j<gennum/aveDegInterval; j++)
        {
            AveOptFit[i][j]=0;
        }

    }

    for(i=0; i<19*relRange; i++)
    {
        storage[i]=' ';
    }

    /*
        Some settings
    */
    pso = fopen("../../vicpso.txt", "a+");

    vfun=4;

    /*
     Initiate the hybrid composition function
    */
    hyb_func4_init(fmax,or_matrix,opt_arr,opt_arr_old);

    edge = hyb_func_goalnedgenvmax[vfun-1][1];
    Vmax = hyb_func_goalnedgenvmax[vfun-1][2];
    goal = hyb_func_goalnedgenvmax[vfun-1][0];


    radRec=0; // Set the radius index counting number

    for(relRadius=relStart; relRadius<relEnd; relRadius+=relStep) // Loop layer 1 (relRadius)
    {

                    int ngoal = 0, avggen = 0; // Goal achieving flag & Average generations

                    for(i = 0; i < testnum; i++)
                    {
                        results[4*i] = gennum;
                        results[4*i+1] = 0;
                        results[4*i+2] = gennum;
                        results[4*i+3] = 0;
                    }


                    for(vtest = 0; vtest < testnum; vtest++)  // Run the experiment
                    {
                        int sum = 1;


                        for(i = 0; i < coornum; i++)
                        {
                            pbest[i] = INT_MAX;
                            gbest[i] = INT_MAX;
                        }

                        //Initialization:

                        for(i = 0; i < dim + 1; i++)
                        {
                            for(j = 0; j < N; j++)
                            {
                                if(i < dim)
                                {
                                    coor[(dim+1)*j+i] = (double)rand() / RAND_MAX *  2 * edge - edge;
                                    vel[(dim+1)*j+i] = (double)rand() / RAND_MAX * 2 * Vmax - Vmax;
                                    pbest[(dim+1)*j+i] = coor[(dim+1)*j+i];
                                    gbest[(dim+1)*j+i] = coor[(dim+1)*j+i];

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
                                    // Initializing pbest
                                    pbest[(dim+1)*j+i] = coor[(dim+1)*j+i];
                                    gbest[(dim+1)*j+i] = coor[(dim+1)*j+i];

                                } // End of the condition i = dim
                            } // for birds



                        } // for dims

                        //neoDiag=0;
                        // Producing the initiating Adjacent Matrix
                        /*for(j=0; j<N; j++)
                            for(k=0; k<N; k++) // Avoid Relapse
                            {
                                dist=getDistance(j,k);
                                if(dist>neoDiag)
                                    neoDiag = dist;
                            }
                            */

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
                        // End of producing the initiating Adjacent Matrix

                     /*  Average Degree
                      for( j=0; j<N ; j++ )
                            for( k=0; k<N ; k++ )
                                if( adjmat[j][k] != 0) // k Can be perceived by j
                                    temp+=1;
                        temp/=N;
                        AveDeg[radRec][0]+=temp; // Got it
                    */

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
                        AveOptFit[radRec][0]+=gb;

                        // End of Updating each bird's gbest in initialization


                        for(vgen = 0; vgen < gennum; vgen++)         // Evolution starts
                        {
                            for(i = 0; i < dim + 1; i++)
                            {
                                for(j = 0; j < N; j++)
                                {
                                    if(i < dim)
                                    {
                                        double r1 = (double)rand() / RAND_MAX * 1;
                                        double r2 = (double)rand() / RAND_MAX * 1;
                                        vel[(dim+1)*j+i] = w * vel[(dim+1)*j+i] + c1 * r1 * (pbest[(dim+1)*j+i] - coor[(dim+1)*j+i])\
                                                           + c2 * r2 * (gbest[(dim+1)*j+i] - coor[(dim+1)*j+i]);

                                        coor[(dim+1)*j+i] += vel[(dim+1)*j+i];

                                    }
                                    else  // i = dim
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
                                            for(l = 0; l < dim + 1; l++)
                                                gbest[(dim+1)*j+l] = coor[(dim+1)*j+l];


                                    } // End of else
                                } // for birds
                            } // for dims

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
                            //if(neoDiag!=0)
                            //printf("%d,%f\n",vgen,neoDiag);
                            // End of producing the initiating Adjacent Matrix

                          /*  //Caculating Average Degree of the generation
                            if( (vgen+1) % aveDegInterval == 0 )
                            {
                                for( j=0; j<N ; j++ )
                                    for( k=0; k<N ; k++ )
                                        if( adjmat[j][k] != 0) // k Can be perceived by j
                                            aveDegree+=1;
                                aveDegree/=N; // Got it
                                AveDeg[radRec][(vgen+1) / aveDegInterval] += aveDegree; // Storing Average Degree of this particular generation
                            }
                        */


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

                            if( (vgen+1) % aveDegInterval == 0 )
                            {
                                AveOptFit[radRec][(vgen+1) / aveDegInterval]+=gb;
                            }

                            if(gb < goal && sum != 0)
                            {
                                results[4*vtest] = (double)vgen+1;
                                results[4*vtest+1] = gb;
                                sum = 0;
                                printf("\ngb:%f,%d\t",gb,vgen);
                                //ngoal++;
                                //break;
                            }

                            if(vgen == 1000)
                                recordInd[vtest]=gb;

                            if(vgen == gennum - 1)
                            {
                                results[4*vtest+2] = gennum;
                                results[4*vtest+3] = gb;
                                printf("lastgb:%f\n",gb);

                                AveEndFit[radRec]+=gb;
                            }


                        } // for vgens


                        printf("\nTest%d for Function%d of relRadius%f is done. The next one's coming at ya.\n\n",vtest,vfun,relRadius);

                    } // for vtest

                    FIndex[radRec]=countFIndex(recordInd);

                    AveEndFit[radRec] = AveEndFit[radRec]/testnum;

                    for(i = 0; i < testnum; i++)
                    {
                        if(results[4*i] != gennum)
                        {
                            ngoal++; //Suc Counter
                            avggen += results[4*i];
                            res += results[4*i+1];
                            AveSucFit[radRec]+=results[4*i+3];
                        }
                    }
                    //printf("%d,%d", avggen, ngoal);
                    if(ngoal == 0)
                    {
                        AveSucFit[radRec] = -1;
                        avggennew = (double)gennum;
                    }
                    else
                    {
                        AveSucFit[radRec] /= ngoal;
                        avggennew = (double)avggen / ngoal;
                    }


                    //printf("AVGEN:%f\t",avggennew);
                    AveOptGen[radRec]=avggennew;
                    //printf("AVRATE:%f",(double)ngoal / testnum * 100);
                    AveOptRate[radRec]=(double)ngoal / testnum * 100;

                    //printf("Function%d is tested under this radius.\n", vfun);


        printf("%1.15f,%1.15f,%f,%f,%f\n",AveEndFit[radRec],AveSucFit[radRec],AveOptGen[radRec],AveOptRate[radRec],FIndex[radRec]); // Record a radius in a row, literally
        fprintf(pso,"%1.15f,%1.15f,%f,%f,%f\n",AveEndFit[radRec],AveSucFit[radRec],AveOptGen[radRec],AveOptRate[radRec],FIndex[radRec]); // Record a radius in a row, literally

        radRec++;

        printf("\n\nRadius %d is tested\n\n",radRec);



    } // End of Loop layer 1 (relRadius)
    printf("\n\n Done! You got it, dude!");

    fclose(pso);

    for(i=0; i<relRange; i++)
        for(j=0; j<=(gennum/aveDegInterval); j++)
        {
            AveOptFit[i][j]=(double)AveOptFit[i][j]/(testnum);
        }

// Get da data in da file!
    sprintf(storage,"Gen,");
    for(i=0; i<relRange; i++)
    {
        sprintf(buffer,"ItemRel%f,",0.24+i*0.02);
        strcat(storage,buffer);
    }


/*
    pso = fopen("vicpsoFunc-v4.4-AveDeg.txt", "w");
    fprintf(pso, "%s\n",storage );
    for(j=0; j<=(gennum/aveDegInterval); j++)
    {
        int inIte;
        fprintf(pso,"%d,",j*10);
        for(inIte=0; inIte<relRange; inIte++)
        {
            fprintf(pso, "%f,",AveDeg[inIte][j]);
        }
        fprintf(pso, "\n");
    }
    fclose(pso)*/

    pso = fopen("../../vicpso-OptFit.txt", "a+");
    fprintf(pso, "%s\n",storage );

    for(j=0; j<=(gennum/aveDegInterval); j++)
    {
        int inIte;
        fprintf(pso,"%d,",j*10);
        for(inIte=0; inIte<relRange; inIte++)
        {
            fprintf(pso, "%f,",AveOptFit[inIte][j]);
        }
        fprintf(pso, "\n");
    }
    fclose(pso);


    return 0;

} // End of main
