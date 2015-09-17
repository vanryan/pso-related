/*
PSO with Changing Perceptive Radius
Version 4.6
2013.11.15


Van Ryan
*/

#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"
#include"string.h"


const int relRange = 17; // The number of values of radius we chose
const int changeInte = 100; // The interval we change the perceptive radius
const int dim = 30; // Number of dimensions
const int N = 50; // Population; This will cause a waste of ram
const int testnum = 1; // Running times for experiments
const int gennum = 5000;  // Number of generations
const int coornum = (dim + 1) * N; // That's a tricky way of storing it
const double kstep = 0.01; // The step of radius parameter change taken by sorted particles
//const double kUpperLim = 1.0;
const double kLowerLim = 0;

const int aveDegInterval=10; // The interval of recording average degree

float coor[coornum];

FILE* pso;

float getDistance(int a,int b)  // a & b stands for the index of the bird
{
    int i;
    float dis,summup=0;
    for(i=0; i<dim; i++)
    {
        summup += pow(coor[(dim+1)*a+i]-coor[(dim+1)*b+i],2);
    }
    dis = sqrt(summup);
    //printf("dis:%f\tsum:%f\n",dis,summup);
    return dis;
}


void Func(int a)
{
    coor[(dim+1)*a+dim] = 0;
    int i;
    for(i = 0; i < dim; i++)
        coor[(dim+1)*a+dim] += pow(coor[(dim+1)*a+i], 2) - 10 * cos(2 * 3.14 * coor[(dim+1)*a+i]) + 10;

}
int partion(double a[2][N],int p,int r){
    //rand
    time_t t;
    srand((unsigned) time(&t));
    int e=rand()%(r-p+1)+p;
    int tem,tem2;
    tem=a[0][e];
    tem2=a[1][e];
    a[0][e]=a[0][r];
    a[1][e]=a[1][r];
    a[0][r]=tem;
    a[1][r]=tem2;
    int x=a[0][r], i=p-1;
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
        int q=partion(a,p,r);
        QuickSort(a,p,q-1);
        QuickSort(a,q+1,r);
    }
}





int main()
{
    srand((unsigned)time(NULL));

    int i, j, k, l, varw, varc1, varc2, vfun=3, vtest, vgen; // iteration vars
    float w, c1, c2, edge, Vmax, goal, gb, res = 0, avggennew=0;  // gb is the global gbest
    float vel[coornum], pbest[coornum], gbest[coornum], results[4*testnum]; // Each bird's got their own gbest
    double absRadius[N],paramK[N]; // Absolute perceptive radius and the radius parameter matrix K
    double relRadius; // A iteration var for radius
    int radRec; // A var to count the index of the radius
    int adjmat[N][N]; // Adjacent Matrix

    char buffer[19]; // A buffer to write the file, 15 is the estimated char numbers
    char storage[19*relRange];


    /* Storage Section */
//The strategy is: we store the 'sum' into 'aves'.
//When the time comes they got loaded into the file, we'll have them become real 'average'

    float AveOptGen[relRange];  // FuncNumber*RadiusNumber: Every 19[radius] with a func
    float AveOptRate[relRange];  // FuncNumber*RadiusNumber: Every 19[radius] with a func
    float AveOptFit[relRange][gennum/aveDegInterval+1]; // The average of the gbest fitness of the interval generations
    float AveEndFit[relRange];// The average of the gbest fitness in the end of each successful test (last generation)
    /* End of Storage Section */


    for(i=0; i<relRange; i++)
    {
        AveOptGen[i]=0;
        AveOptRate[i]=0;
        AveEndFit[i]=0;
        for(j=0; j<gennum/aveDegInterval; j++)
        {
  //          AveDeg[i][j]=0;
            AveOptFit[i][j]=0;
        }

    }

    for(i=0; i<19*relRange; i++)
    {
        storage[i]=' ';
    }


    radRec=0; // Set the radius index counting number

    for(relRadius=0.2; relRadius<1.01; relRadius+=0.05) // Loop layer 1 (relRadius)
    {

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
                        double neoDiag=0;
                        double dist=0;
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

                                    Func(j);
                                    // Initializing pbest
                                    pbest[(dim+1)*j+i] = coor[(dim+1)*j+i];
                                    gbest[(dim+1)*j+i] = coor[(dim+1)*j+i];

                                } // End of the condition i = dim
                            } // for birds



                        } // for dims


                        absRadius[0] = relRadius * sqrt(dim * pow(2 * edge , 2));
                        for(j=0;j<N;j++){
                            paramK[j]=relRadius;
                            absRadius[j]=absRadius[0];
                        }


                        for(j=0; j<N; j++)
                            for(k=0; k<N; k++) // Avoid Relapse
                            {
                                dist=getDistance(j,k);
                                if( dist <= absRadius[j] && j!=k ) // k Can be perceived by j
                                {
                                    adjmat[j][k]=1;
                                }
                                else
                                {
                                    adjmat[j][k]=0;
                                }


                            }
                        // End of producing the initiating Adjacent Matrix


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
                                        float r1 = (double)rand() / RAND_MAX * 1;
                                        float r2 = (double)rand() / RAND_MAX * 1;
                                        vel[(dim+1)*j+i] = w * vel[(dim+1)*j+i] + c1 * r1 * (pbest[(dim+1)*j+i] - coor[(dim+1)*j+i])\
                                                           + c2 * r2 * (gbest[(dim+1)*j+i] - coor[(dim+1)*j+i]);

                                        coor[(dim+1)*j+i] += vel[(dim+1)*j+i];

                                    }
                                    else  // i = dim
                                    {

                                        Func(j);
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

                            if( (vgen+1) % (changeInte) == 0 && vgen!=0)
                            {
                                double fit[2][N];
                                for(j=0;j<N;j++){
                                     fit[0][j]=coor[(dim+1)*j+dim];  // Fitness
                                     fit[1][j]=j;  // The index of particle
                                }
                                QuickSort(fit,0,N-1);  // Done sorting the fitness while keeping the index
                                // Changing the parameter k
                                // The quicksort gives a small to large sequence!
                                for(j=0;j<N;j++){
                                    if(j>=35){
                                        paramK[(int)fit[1][j]]-=kstep;
                                        if( paramK[(int)fit[1][j]] < kLowerLim)
                                             paramK[(int)fit[1][j]] = kLowerLim;
                                    }
                                    if(j<15){
                                        paramK[(int)fit[1][j]]+=kstep;
                                    }
                                    absRadius[(int)fit[1][j]]=paramK[(int)fit[1][j]]*neoDiag;
                                    //printf("%d,",(int)fit[1][j]);
                                }
                                //printf("\n");
                                //printf("%d,%f\n",vgen,neoDiag);
                            }

                            neoDiag=0;

                            // Producing the initiating Adjacent Matrix
                            for( j=0; j<N ; j++ )
                                for( k=0; k<N; k++) // Avoid Relapse
                                {
                                    dist=getDistance(j,k);
                                    if( dist <= absRadius[j] && j!=k ) // k Can be perceived by j
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
                            if(vgen == gennum - 1)
                            {
                                results[4*vtest+2] = gennum;
                                results[4*vtest+3] = gb;
                                printf("lastgb:%f\n",gb);
                            }


                        } // for vgens


                        printf("\nTest%d for Function%d of relRadius%f is done. The next one's coming at ya.\n\n",vtest,vfun,relRadius);

                    } // for vtest

                    for(i = 0; i < testnum; i++)
                    {
                        if(results[4*i] != gennum)
                        {
                            ngoal++;
                            avggen += results[4*i];
                            res += results[4*i+1];
                            AveEndFit[radRec]+=results[4*i+3];
                        }
                    }
                    printf("%d,%d", avggen, ngoal);
                    if(ngoal == 0)
                    {
                        AveEndFit[radRec]=-1;
                        avggennew = (double)gennum;
                    }
                    else
                    {
                        AveEndFit[radRec] /= ngoal;
                        avggennew = (double)avggen / ngoal;
                    }


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

    for(i=0; i<relRange; i++)
        for(j=0; j<=(gennum/aveDegInterval); j++)
        {
            AveOptFit[i][j]=(double)AveOptFit[i][j]/(testnum);
        }

// Get da data in da file!
    sprintf(storage,"Gen,");
    for(i=0; i<relRange; i++)
    {
        sprintf(buffer,"ItemRel%f,",0.2+i*0.05);
        strcat(storage,buffer);
    }


    pso = fopen("c:/pso/dynaRadKSort1.5/100/vicpsoFunc3-main.txt", "w");
    fprintf(pso, "RelRadius,AveOptGen,AveOptRate\n");
    for(i=0; i<relRange; i++)
    {
        fprintf(pso, "%f,%f,%f\n",0.2+i*0.05,AveOptGen[i],AveOptRate[i]);
    }
    fclose(pso);


    pso = fopen("c:/pso/dynaRadKSort1.5/100/vicpsoFunc3-OptFit.txt", "w");
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

    pso = fopen("c:/pso/dynaRadKSort1.5/100/vicpsoFunc3-EndSucFit.txt","w");
    for(i=0; i<relRange;i++){
        fprintf(pso,"%f,%f,\n",0.2+i*0.05,AveEndFit[i]);
    }
    fclose(pso);

    return 0;

} // End of main
