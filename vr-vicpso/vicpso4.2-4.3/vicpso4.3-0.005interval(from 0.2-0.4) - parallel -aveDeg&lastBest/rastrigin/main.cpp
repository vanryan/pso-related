/*
PSO with perceptive radius
Version 1
2013.10.22

Changing Probing Radius
from 0.2-0.45 stepped by interval 0.05

Investigating the leap of stagnant average degree during this range

Van Ryan
*/

#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"
#include"string.h"

//const float relRadius = 0.3; // Perceptive radius (Relative)
const int dim = 30; // Number of dimensions
//const int N = 15;
const int N = 50; // Population; This will cause a waste of ram
const int testnum = 25; // Running times for experiments
const int gennum = 5000;  // Number of generations
const int coornum = (dim + 1) * N; // That's a tricky way of storing it

const int aveDegInterval=10; // The interval of recording average degree

float coor[coornum];

FILE* pso;

float getDistance(int a,int b){ // a & b stands for the index of the bird
    int i;
    float dis,sum=0;
    for(i=0;i<dim;i++){
        sum += pow(coor[(dim+1)*a+i]-coor[(dim+1)*b+i],2);
    }
    dis = sqrt(sum);
    //printf("dis:%f\tsum:%f\n",dis,sum);
    return dis;
}



void Fun3(int a)                             //Rastrigin
{
    coor[(dim+1)*a+dim] = 0;
    int i;
    for(i = 0; i < dim; i++)
        coor[(dim+1)*a+dim] += pow(coor[(dim+1)*a+i], 2) - 10 * cos(2 * 3.14 * coor[(dim+1)*a+i]) + 10;
}


int main()
{
    srand((unsigned)time(NULL));

    int i, j, k, l, varw, varc1, varc2, vfun=3, vtest, vgen, temp; // iteration vars
    float w, c1, c2, edge, Vmax, goal, gb, res = 0, avggennew=0;  // gb is the global gbest
    float vel[coornum], pbest[coornum], gbest[coornum], results[4*testnum]; // Each bird's got their own gbest
    float absRadius; // Absolute perceptive radius
    float relRadius; // A iteration var for radius
    int radRec; // A var to count the index of the radius
    int adjmat[N][N]; // Adjacent Matrix

    char buffer[19]; // A buffer to write the file, 15 is the estimated char numbers
    char storage[19*51];


    /* Storage Section */
//The strategy is: we store the 'sum' into 'aves'.
//When the time comes they got loaded into the file, we'll have them become real 'average'

    float AveOptGen[51];
    float AveOptRate[51];
    float AveDeg[51][gennum/aveDegInterval+1];  // [Every 10 gen]
                                            //So the truth is, we're getting the average(a matrix) of average(a test)
    float AveOptFit[51][gennum/aveDegInterval+1]; // The average of the gbest fitness in the end of each test (last generation)
    /* End of Storage Section */


    for(i=0;i<51;i++){
        AveOptGen[i]=0;
        AveOptRate[i]=0;
        for(j=0;j<gennum/aveDegInterval;j++){
            AveDeg[i][j]=0;
            AveOptFit[i][j]=0;
        }
    }

    for(i=0;i<19*14*5;i++){
        storage[i]=' ';
    }


    radRec=0; // Set the radius index counting number

for(relRadius=0.1;relRadius<1.41;relRadius+=0.02){ // Loop layer 1 (relRadius)

    w = 0.6;
    for(varw = 0; varw < 1; varw++)
    {
        c1 = 1.7;
        for(varc1 = 0; varc1 < 1; varc1++)
        {
            c2 = 1.7;
            for(varc2 = 0; varc2 < 1 ; varc2++)
            {


                        edge = 5.12;
                        Vmax = edge;
                        goal = 100;


                    absRadius = relRadius * sqrt(dim * pow(2 * edge , 2));
                    //printf("absRadius:%f\n",absRadius);

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

                        /*for(i = 0; i < N; i++)
                            for(j = 0; j < N; j++)
                                adjmat[i][j]=0;
                        // Finished initializing adjacent matrix
                        */

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
                                { // i=dim

                                        Fun3(j);

                                    // Initializing pbest
                                    pbest[(dim+1)*j+i] = coor[(dim+1)*j+i];
                                    gbest[(dim+1)*j+i] = coor[(dim+1)*j+i];
                                    /*
                                    if(coor[(dim+1)*j+i] < gbest[(dim+1)*j+i])
                                    {
                                        for(l = 0; l < dim + 1; l++)
                                            gbest[(dim+1)*j+l] = coor[(dim+1)*j+l];

                                    }
                                    for(k = 0; k < N; k++)
                                        if( gbest[(dim+1)*j+i] < gbest[(dim+1)*k+i])
                                        {

                                            for(l = 0; l < dim + 1; l++)
                                                gbest[(dim+1)*k+l] = gbest[(dim+1)*j+l];
                                        }
                                    */
                                } // End of the condition i = dim
                            } // for birds



                        } // for dims

                        // Producing the initiating Adjacent Matrix
                        for(j=0;j<N;j++)
                            for(k=0;k<N;k++) // Avoid Relapse
                            {
                                if( getDistance(j,k) <= absRadius && j!=k ) // k Can be perceived by j
                                    {
                                        adjmat[j][k]=1;
                                    }
                                else{
                                    adjmat[j][k]=0;
                                }
                            }
                        // End of producing the initiating Adjacent Matrix

                        for( j=0; j<N ; j++ )
                                for( k=0; k<N ; k++ )
                                    if( adjmat[j][k] != 0) // k Can be perceived by j
                                        temp+=1;
                        temp/=N;
                        AveDeg[radRec][0]+=temp; // Got it

                        // Updating each bird's gbest:
                        for(j=0;j<N;j++){
                             for(k=0;k<N;k++)
                                if( adjmat[j][k] != 0 ){// k Can be perceived by j
                                        // Updating the gbest of j
                                        //printf("%f\t",getDistance(j,k));
                                        if( coor[(dim+1)*k+dim] < gbest[(dim+1)*j+dim] ){
                                            for(l=0;l<=dim;l++)
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


                        for(vgen = 0; vgen < gennum; vgen++)         // Evolution
                        {
                            float aveDegree=0; // Average Degree


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

                                        coor[(dim+1)*j+i] += vel[(dim+1)*j+i];

                                    }
                                    else  // i = dim
                                    {

                                            Fun3(j);

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
                                        int dist=getDistance(j,k);
                                        if( dist <= absRadius && j!=k ) // k Can be perceived by j
                                        {
                                            adjmat[j][k]=1;
                                        }
                                        else{
                                                adjmat[j][k]=0;
                                        }
                                    }
                            // End of producing the initiating Adjacent Matrix

                            //Caculating Average Degree of the generation
                            if( (vgen+1) % aveDegInterval == 0 ){
                                for( j=0; j<N ;j++ )
                                    for( k=0; k<N ;k++ )
                                        if( adjmat[j][k] != 0) // k Can be perceived by j
                                            aveDegree+=1;
                                aveDegree/=N; // Got it
                                AveDeg[radRec][(vgen+1) / aveDegInterval] += aveDegree; // Storing Average Degree of this particular generation
                            }


                            // Updating each bird's gbest:
                            for( j=0; j<N ;j++ ){
                                for( k=0; k<N ;k++ )
                                    if( adjmat[j][k] != 0 ) // k Can be perceived by j
                                        // Updating the gbest of j
                                        if( coor[(dim+1)*k+dim] < gbest[(dim+1)*j+dim] )
                                            for(l=0;l<=dim;l++)
                                                gbest[(dim+1)*j+l] = coor[(dim+1)*k+l];
                                //printf("gbest of %d:%f\nfitness of %d:%f\n",j,gbest[(dim+1)*j+dim],j,coor[(dim+1)*j+dim]);
                            }
                            // End of Updating each bird's gbest

                            gb = gbest[dim]; // Updating the global gbest
                            for(i = 1; i < N; i++)
                                if(gbest[(dim+1)*i+dim] < gb)
                                    gb = gbest[(dim+1)*i+dim];

                            if( (vgen+1) % aveDegInterval == 0 ){
                                AveOptFit[radRec][(vgen+1) / aveDegInterval]+=gb;
                            }

                            if(gb < goal && sum != 0)
                            {
                                results[4*vtest] = (double)vgen+1;
                                results[4*vtest+1] = gb;
                                sum = 0;
                                //ngoal++;
                                //break;
                            }
                            if(vgen == gennum - 1)
                            {
                                results[4*vtest+2] = gennum;
                                results[4*vtest+3] = gb;

                            }


                        } // for vgens


                        printf("\nA test for Function%d is done. The next one's coming at ya.\n\n",vfun);

                    } // for vtest

                    for(i = 0; i < testnum; i++)
                    {
                        if(results[4*i] != gennum)
                        {
                            ngoal++;
                            avggen += results[4*i];
                            res += results[4*i+1];
                        }
                    }
                    printf("%d,%d", avggen, ngoal);
                    if(ngoal == 0)
                        avggennew = (double)gennum;
                    else
                        avggennew = (double)avggen / ngoal;

                    printf("AVGEN:%f\t",avggennew);
                    AveOptGen[radRec]=avggennew;
                    printf("AVRATE:%f",(double)ngoal / testnum * 100);
                    AveOptRate[radRec]=(double)ngoal / testnum * 100;


                    printf("Function%d is tested under this radius.\n", vfun);

                c2 += 0.1;
            }
            c1 += 0.1;
        }
        w += 0.1;
    }
//} // End of Loop layer 2 (population)

    radRec++;

    printf("\n\nRadius %d is tested\n\n",radRec);

} // End of Loop layer 1 (relRadius)
    printf("\n\n Done! You got it, dude!");

    for(i=0;i<51;i++)
         for(j=0;j<=(gennum/aveDegInterval);j++){
            AveDeg[i][j]=(double)AveDeg[i][j]/testnum;
            AveOptFit[i][j]=(double)AveOptFit[i][j]/testnum;
        }

// Get da data in da file!
    for(i=0;i<51;i++){
        sprintf(buffer,"Gen%d,AvDeg%d,",i,i);
        strcat(storage,buffer);
    }



    pso = fopen("vicpsoFunc3-adjmat.txt", "w");
    fprintf(pso, "RelRadius,AveOptGen,AveOptRate,AveDeg\n");
    for(i=0;i<51;i++){
        fprintf(pso, "%f,%f,%f\n",0.2+i*0.005,AveOptGen[i],AveOptRate[i]);
    }
    fclose(pso);

    pso = fopen("vicpsoFunc3-adjmat-AveDeg.txt", "w");
    fprintf(pso, "%s\n",storage );

        for(j=0;j<=(gennum/aveDegInterval);j++){
            int inIte;
            for(inIte=0;inIte<51;inIte++){
                fprintf(pso, "%d,%f,",j*10,AveDeg[inIte][j]);
            }
            fprintf(pso, "\n");
        }


    fclose(pso);

    pso = fopen("vicpsoFunc3-OptFit.txt", "w");
    fprintf(pso, "%s\n",storage );

        for(j=0;j<=(gennum/aveDegInterval);j++){
            int inIte;
            for(inIte=0;inIte<51;inIte++){
                fprintf(pso, "%d,%f,",j*10,AveOptFit[inIte][j]);
            }
            fprintf(pso, "\n");
        }
    fclose(pso);


    return 0;

} // End of main
