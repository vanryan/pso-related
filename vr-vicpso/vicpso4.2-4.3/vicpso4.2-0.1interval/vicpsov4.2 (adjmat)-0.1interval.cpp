/*
PSO with perceptive radius
Version 1
2013.10.11
*/

#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

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

void Fun3(int a)                             //Rastrigin
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
        b += pow(coor[(dim+1)*a+i],2) / 4000;
        c *= cos(coor[(dim+1)*a+i]/ sqrt(i+1));
    }
    coor[(dim+1)*a+dim] = 1 + b - c;
}

void Fun5(int a)
{
    coor[(dim+1)*a+dim] = 0;
    // Not written

}


int main()
{
    srand((unsigned)time(NULL));

    int i, j, k, l, varw, varc1, varc2, vfun, vtest, vgen; // iteration vars
    float w, c1, c2, edge, Vmax, goal, gb, res = 0, avggennew=0;  // gb is the global gbest
    float vel[coornum], pbest[coornum], gbest[coornum], results[4*testnum]; // Each bird's got their own gbest
    float absRadius; // Absolute perceptive radius
    float relRadius; // A iteration var for radius
    int radRec; // A var to count the index of the radius
    int adjmat[N][N]; // Adjacent Matrix


    /* Storage Section */
//The strategy is: we store the 'sum' into 'aves'.
//When the time comes they got loaded into the file, we'll have them become real 'average'

    float AveOptGen[52];  // FuncNumber*RadiusNumber: Every 19[radius] with a func
    float AveOptRate[52];  // FuncNumber*RadiusNumber: Every 19[radius] with a func
    float AveDeg[52][gennum/aveDegInterval];  // [FuncNumber*RadiusNumber][Every 10 gen]
                                            //So the truth is, we're getting the average(a matrix) of average(a test)
    /* End of Storage Section */


    for(i=0;i<52;i++){
        AveOptGen[i]=0;
        AveOptRate[i]=0;
        for(j=0;j<gennum/aveDegInterval;j++)
            AveDeg[i][j]=0;
    }


    radRec=0; // Set the radius index counting number

for(relRadius=0.1;relRadius<1.4;relRadius+=0.1){ // Loop layer 1 (relRadius)


    w = 0.6;
    for(varw = 0; varw < 1; varw++)
    {
        c1 = 1.7;
        for(varc1 = 0; varc1 < 1; varc1++)
        {
            c2 = 1.7;
            for(varc2 = 0; varc2 < 1 ; varc2++)
            {
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
                                    if(vfun == 1)
                                        Fun1(j);
                                    if(vfun == 2)
                                        Fun2(j);
                                    if(vfun == 3)
                                        Fun3(j);
                                    if(vfun == 4)
                                        Fun4(j);

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

                        AveDeg[radRec*4+vfun-1][0]=0;
                            for( j=0; j<N ; j++ )
                                for( k=0; k<N ; k++ )
                                    if( adjmat[j][k] != 0) // k Can be perceived by j
                                        AveDeg[radRec*4+vfun-1][0]+=1;
                        AveDeg[radRec*4+vfun-1][0]/=N; // Got it

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
                                        if(vfun == 1)
                                            Fun1(j);
                                        if(vfun == 2)
                                            Fun2(j);
                                        if(vfun == 3)
                                            Fun3(j);
                                        if(vfun == 4)
                                            Fun4(j);

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
                                AveDeg[radRec*4+vfun-1][(vgen+1) / aveDegInterval] += aveDegree; // Storing Average Degree of this particular generation
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

                        printf("Generations taken:%f\n",results[4*vtest] );
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
                    AveOptGen[radRec*4+vfun-1]=avggennew;
                    printf("AVRATE:%f",(double)ngoal / testnum * 100);
                    AveOptRate[radRec*4+vfun-1]=(double)ngoal / testnum * 100;

                    printf("Function%d is tested under this radius.\n", vfun);
                } // End of loop for functions
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

    for(i=0;i<13;i++)
        for(j=0;j<gennum/aveDegInterval;j++)
            AveDeg[i][j]=(double)AveDeg[i][j]/testnum;

// Get da data in da file!
    pso = fopen("vicpsoFunc1-adjmat.txt", "w");
    fprintf(pso, "RelRadius,AveOptGen,AveOptRate\n");
    for(i=0;i<13;i++){
        fprintf(pso, "%f,%f,%f\n",0.1+i*0.1,AveOptGen[i*4],AveOptRate[i*4]);
    }
    fclose(pso);

    pso = fopen("vicpsoFunc1-adjmat-AveDeg.txt", "w");
    fprintf(pso, "Gen1,AveDeg1,Gen2,AveDeg2,Gen3,AveDeg3,Gen4,AveDeg4,Gen5,AveDeg5,Gen6,AveDeg6,Gen7,AveDeg7,Gen8,AveDeg8,Gen9,AveDeg9,Gen10,AveDeg10,Gen11,AveDeg11,Gen12,AveDeg12,Gen13,AveDeg13\n");
    for(j=0;j<=gennum/aveDegInterval;j++)
        fprintf(pso, "%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f\n",\
j*10,AveDeg[0][j],j*10,AveDeg[4][j],j*10,AveDeg[8][j],j*10,AveDeg[12][j],j*10,AveDeg[16][j],j*10,AveDeg[20][j],j*10,AveDeg[24][j],j*10,AveDeg[28][j],j*10,AveDeg[32][j],j*10,AveDeg[36][j],j*10,AveDeg[40][j],j*10,AveDeg[44][j],j*10,AveDeg[48][j]);
    fclose(pso);

    pso = fopen("vicpsoFunc2-adjmat.txt", "w");
    fprintf(pso, "RelRadius,AveOptGen,AveOptRate,AveDeg\n");
    for(i=0;i<13;i++){
        fprintf(pso, "%f,%f,%f\n",0.1+i*0.1,AveOptGen[i*4+1],AveOptRate[i*4+1]);
    }
    fclose(pso);

    pso = fopen("vicpsoFunc2-adjmat-AveDeg.txt", "w");
    fprintf(pso, "Gen1,AveDeg1,Gen2,AveDeg2,Gen3,AveDeg3,Gen4,AveDeg4,Gen5,AveDeg5,Gen6,AveDeg6,Gen7,AveDeg7,Gen8,AveDeg8,Gen9,AveDeg9,Gen10,AveDeg10,Gen11,AveDeg11,Gen12,AveDeg12,Gen13,AveDeg13\n");
    for(j=0;j<=gennum/aveDegInterval;j++)
       fprintf(pso, "%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f\n",\
         j*10,AveDeg[1][j],j*10,AveDeg[5][j],j*10,AveDeg[9][j],j*10,
         AveDeg[13][j],j*10,AveDeg[17][j],j*10,AveDeg[21][j],j*10,AveDeg[25][j],j*10,AveDeg[29][j],
         j*10,AveDeg[33][j],j*10,AveDeg[37][j],j*10,AveDeg[41][j],j*10,AveDeg[45][j],j*10,AveDeg[49][j]);
    pso = fopen("vicpsoFunc3-adjmat.txt", "w");
    fprintf(pso, "RelRadius,AveOptGen,AveOptRate,AveDeg\n");
    for(i=0;i<13;i++){
        fprintf(pso, "%f,%f,%f\n",0.1+i*0.1,AveOptGen[i*4+2],AveOptRate[i*4+2]);
    }
    fclose(pso);

    pso = fopen("vicpsoFunc3-adjmat-AveDeg.txt", "w");
    fprintf(pso, "Gen1,AveDeg1,Gen2,AveDeg2,Gen3,AveDeg3,Gen4,AveDeg4,Gen5,AveDeg5,Gen6,AveDeg6,Gen7,AveDeg7,Gen8,AveDeg8,Gen9,AveDeg9,Gen10,AveDeg10,Gen11,AveDeg11,Gen12,AveDeg12,Gen13,AveDeg13\n");
    for(j=0;j<=gennum/aveDegInterval;j++)
        fprintf(pso, "%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f\n",\
         j*10,AveDeg[2][j],j*10,AveDeg[6][j],j*10,AveDeg[10][j],j*10,AveDeg[14][j],j*10,AveDeg[18][j],j*10,
         AveDeg[22][j],j*10,AveDeg[26][j],j*10,AveDeg[30][j],j*10,AveDeg[34][j],j*10,AveDeg[38][j],j*10,AveDeg[42][j],j*10,AveDeg[46][j],j*10,AveDeg[50][j]);
    fclose(pso);

    pso = fopen("vicpsoFunc4-adjmat.txt", "w");
    fprintf(pso, "RelRadius,AveOptGen,AveOptRate,AveDeg\n");
    for(i=0;i<13;i++){
        fprintf(pso, "%f,%f,%f\n",0.1+i*0.1,AveOptGen[i*4+3],AveOptRate[i*4+3]);
    }
    fclose(pso);

    pso = fopen("vicpsoFunc4-adjmat-AveDeg.txt", "w");
    fprintf(pso, "Gen1,AveDeg1,Gen2,AveDeg2,Gen3,AveDeg3,Gen4,AveDeg4,Gen5,AveDeg5,Gen6,AveDeg6,Gen7,AveDeg7,Gen8,AveDeg8,Gen9,AveDeg9,Gen10,AveDeg10,Gen11,AveDeg11,Gen12,AveDeg12,Gen13,AveDeg13\n");
    for(j=0;j<=gennum/aveDegInterval;j++)
        fprintf(pso, "%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f,%d,%f\n",\
         j*10,AveDeg[3][j],j*10,AveDeg[7][j],j*10,AveDeg[11][j],j*10,AveDeg[15][j],j*10,AveDeg[19][j],j*10,AveDeg[23][j],
         j*10,AveDeg[27][j],j*10,AveDeg[31][j],j*10,AveDeg[35][j],j*10,AveDeg[39][j],j*10,AveDeg[43][j],j*10,AveDeg[47][j],j*10,AveDeg[51][j]);
    fclose(pso);

    return 0;

} // End of main
