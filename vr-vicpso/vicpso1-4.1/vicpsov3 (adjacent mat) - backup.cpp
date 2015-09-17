/*
PSO with perceptive radius
Version 1
2013.09.22
*/

#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

//const float relRadius = 0.3; // Perceptive radius (Relative)
const int dim = 30; // Number of dimensions
//const int N = 15;
int N = 60; // Population
const int testnum = 30; // Running times for experiments
const int gennum = 5000;  // Number of generations
const int coornum = (dim + 1) * N; // That's a tricky way of storing it

float coor[coornum];

FILE* pso;

float getDistance(int a,int b){ // a & b stands for the index of the bird
    int i;
    float dis,sum=0;
    for(i=0;i<dim;i++){
        sum+=pow(coor[(dim+1)*a+i]-coor[(dim+1)*b+i],2);
    }
    dis=sqrt(sum);
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

    int i, j, k, l, varw, varc1, varc2, vfun, vtest, vgen ,n; // iteration vars
    float w, c1, c2, edge, Vmax, goal, gb, res = 0, avggennew=0;  // gb is the global gbest
    float vel[coornum], pbest[coornum], gbest[coornum], results[2*testnum]; // Each bird's got their own gbest
    float absRadius; // Absolute perceptive radius
    float relRadius; // a iteration var for radius
    int adjmat[N][N]; // Adjacent Matrix


    for(i = 0; i < testnum; i++)
    {
        results[2*i] = gennum;
        results[2*i+1] = 0;
    }

    pso = fopen("vicpso1-4.txt", "a+");
    fprintf(pso, "Vicpso:\n");

for(relRadius=0;relRadius<=1.3;relRadius+=0.1){ // Loop layer 1 (relRadius)

    for(n=0;n<3;n++){ // Loop layer 2 (population)

        switch(n){
        case 0:N=15; break;
        case 1:N=30; break;
        case 2:N=60; break;
        default:N=15; 
        }


    w = 0.6;
    for(varw = 0; varw < 1; varw++)
    {
        c1 = 1.7;
        for(varc1 = 0; varc1 < 1; varc1++)
        {
            c2 = 1.7;
            for(varc2 = 0; varc2 < 1 ; varc2++)
            {
                for(vfun = 1; vfun < 4; vfun++)
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

                    absRadius = relRadius*sqrt(dim * pow(edge,2));
                    //printf("absRadius:%f\n",absRadius);

                    int ngoal = 0, avggen = 0; // Goal achieving flag & Average generations

                    for(vtest = 0; vtest < testnum; vtest++)  // Run the experiment
                    {

                        for(i = 0; i < n; i++)
                            for(j = 0; j < n; j++)
                                adjmat[i][j]=0;
                        // Finished initializing adjacent matrix

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
                            for(k=j+1;k<N;k++) // Avoid Relapse
                                if( getDistance(j,k) <= absRadius) // k Can be perceived by j
                                    {
                                        adjmat[j][k]=1;
                                        adjmat[k][j]=1;
                                    }
                        // End of producing the initiating Adjacent Matrix

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
                        // End of Updating each bird's gbest


                        for(vgen = 0; vgen < gennum; vgen++)         // Evolution
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

                                        /*
                                        // Updating each bird's own gbest:
                                        if(coor[(dim+1)*j+i] < gbest[(dim+1)*j+i])
                                        {
                                            for(l = 0; l < dim + 1; l++)
                                                gbest[(dim+1)*j+l] = coor[(dim+1)*j+l];
                                        }
                                        for(k = 0; k < N; k++)
                                            if( gbest[(dim+1)*j+i] < gbest[(dim+1)*k+i])
                                                for(l = 0; l < dim + 1; l++)
                                                    gbest[(dim+1)*k+l] = gbest[(dim+1)*j+l];
                                        */
                                    } // End of else
                                } // for birds
                            } // for dims

                            // Producing the initiating Adjacent Matrix
                            for( j=0; j<N ; j++ )
                                for( k=j+1; k<N; k++) // Avoid Relapse
                                    if( getDistance(j,k) <= absRadius) // k Can be perceived by j
                                        {
                                            adjmat[j][k]=1;
                                            adjmat[k][j]=1;
                                        }
                            // End of producing the initiating Adjacent Matrix

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

                        } // for vgens
                    } // for vtest

                    for(i = 0; i < testnum; i++)
                    {
                        if(results[2*i] != gennum - 1)
                        {
                            ngoal++;
                            avggen += results[2*i];
                            res += results[2*i+1];
                        }
                    }
                    printf("%d,%d", avggen, ngoal);
                    if(ngoal == 0)
                        avggennew = (double)gennum - 1;
                    else
                        avggennew = (double)avggen / ngoal;
                    fprintf(pso, "\nfunction:%d\nrelRadius:%f\nPopulation:%d\n平均代数为：%f\n达优率为：%f%%",vfun,relRadius,N, avggennew, (double)ngoal / testnum * 100);
                }
                c2 += 0.1;
            }
            c1 += 0.1;
        }
        w += 0.1;
    }
} // End of Loop layer 2 (population)
} // End of Loop layer 1 (relRadius)    
    printf("\n\n Done! You've got it, dude!");
    fclose(pso);
}
