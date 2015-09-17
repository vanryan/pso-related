/*
 FIPSO v2.0
 Van Ryan

 A particle learns from random neighbors ( the number of the learnees is set [we still use vlneighs] )

 2014.02.20
*/
#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

const int dim = 30;
const int N = 50;
const int testnum = 20;
const int gennum = 5000;
const int coornum = (dim + 1) * N;

double pbest[coornum];
double update[coornum];

int randlneighs[N];

double coor[coornum], stdindmemo[testnum];

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

void Fun5(int a)
{
    coor[(dim+1)*a+dim] = 0;

}

void choose_rand_neighbors(int lneighs)
{
    int counter=0;
    int i,flag=1,inttemp;
    double tmp;
    for(i=0; i<lneighs; i++)
        randlneighs[i]=0; // initializing the indices of learnees
    while(counter<lneighs)
    {
        tmp = (double)rand() / RAND_MAX;
        tmp*= (N-1);
        inttemp=(int)tmp;
        flag=1;

        for(i=0; i<counter; i++) //repeat check
        {
            if(randlneighs[i]==inttemp) // We already have it
            {
                flag=0;
                break;
            }
        }

        if(flag!=0)// succeeded in the repeat check
        {
            randlneighs[counter]=inttemp;
            counter++;
        }
        // else if failed the repeat check: go on producing and checking
    }
}

void get_vel_update(int lneighs, int gbindex)
{
    int i,j,k;
    double coef[N];
    double phi = 4.1;

    double PBNINDEX[2][N]; // Pbest & Index

    for(j = 0; j < N; j++)
    {
        choose_rand_neighbors(lneighs);

        for(k=0; k<lneighs; k++)
        {
            PBNINDEX[0][k] = pbest[(dim+1)*randlneighs[k]+dim];  // pbest
            PBNINDEX[1][k] = randlneighs[k];  // The index of particle
        }

        for(i = 0; i < dim; i++)
        {
            double sum_coef = 0;
            double pos_ave = 0;

            for( k = 0; k < lneighs; k++) // One particle can only learn from a certain number of neighbors
            {
                coef[(int)PBNINDEX[1][k]] = (double)rand() / RAND_MAX;
                sum_coef += coef[(int)PBNINDEX[1][k]];
                pos_ave += coef[(int)PBNINDEX[1][k]] * pbest[(dim+1)*((int)PBNINDEX[1][k])+i];
            }

            // The particle learns from gbest
            double temp = (double)rand() / RAND_MAX;
            sum_coef += temp;
            pos_ave += temp * pbest[(dim+1)*gbindex+i];

            // The particle should also learn from itself
            temp = (double)rand() / RAND_MAX;
            sum_coef += temp;
            pos_ave += temp * pbest[(dim+1)*j+i];

            // After learning
            pos_ave /= sum_coef;
            double var_phi = sum_coef * phi / (lneighs+2);

            double social_central_tendency = pos_ave - coor[(dim+1)*j+i];
            update[(dim+1)*j+i] = var_phi * social_central_tendency ;
        }
    }
}

int main()
{
    //printf("a");
    srand((unsigned)time(NULL));
    int i, j, l, vfun, vtest, vgen, vlneighs;
    float x, edge, Vmax, goal, gb, res = 0,avggennew=0;
    float vel[coornum], results[2*testnum], stdindex[N-1];
    int sucCounter;

// vlneighs:The number of neighbors a paricle learns from (self-excluded)

    float AveEndFit[4];
    float AveAllFit[4];

    //Storage
    float Sto_SucEndFit[N-1]; // N-1 circumstances
    float Sto_SucGens[N-1];
    float Sto_SucRate[N-1];

    for(i=0;i<N-1;i++){
        Sto_SucEndFit[i]=-1;
        Sto_SucGens[i]=-1;
        Sto_SucRate[i]=-1;
    }

    pso = fopen("lipso2.3-all-f1.txt", "a+");
    //fprintf(pso, "num\tw\tc1\tc2\tfuncnum\tavggen\tratio\tfinalresult\n");

    x = 0.7298;

    for(vlneighs =0; vlneighs< N-1; vlneighs++)
    {
        fprintf(pso, "lneighs:%d\n",vlneighs);
        for(i = 0; i < 4; i++){
            AveAllFit[i]=0;
            AveEndFit[i]=0;
        }

        for(vfun = 1; vfun < 2; vfun++)
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

                for(i = 0; i < dim + 1; i++)               //Á£×Ó³õÊ¼»¯£»
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
/*
                int gbp;
                gb = pbest[dim];
                for(i = 1; i < N; i++)
                    if(pbest[(dim+1)*i+dim] < gb)
                    {
                        gb = pbest[(dim+1)*i+dim];
                        gbp = i;
                    }
*/


                for(vgen = 0; vgen < gennum; vgen++)
                {

                    // Finding gbest
                    double tempgb=pbest[dim];
                    int gbindex=0;
                    for( i = 1; i < N; i++)
                        if(pbest[(dim+1)*i+dim]<tempgb ){
                            tempgb=pbest[(dim+1)*i+dim];
                            gbindex=i;
                        }

                    get_vel_update(vlneighs, gbindex);

                    for(i = 0; i < dim + 1; i++)
                    {
                        for(j = 0; j < N; j++)
                        {
                            if(i < dim)
                            {
                                //double r1=(double)rand() / RAND_MAX;
                                //double r2=(double)rand() / RAND_MAX;
                                //vel[(dim+1)*j+i] = x * ( vel[(dim+1)*j+i] + phi*( r1*( pbest[(dim+1)*gbp+i]-coor[(dim+1)*j+i]) + r2*(pbest[(dim+1)*j+i]-coor[(dim+1)*j+i]) )/2 );
                                //vel[(dim+1)*j+i] = x * ( vel[(dim+1)*j+i] +  2.05*r1*(pbest[(dim+1)*j+i]-coor[(dim+1)*j+i])+2.05*r2*(pbest[(dim+1)*gbp+i]-coor[(dim+1)*j+i]));

                                vel[(dim+1)*j+i] = x * (vel[(dim+1)*j+i]+update[(dim+1)*j+i]);
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
                        {
                            gb = pbest[(dim+1)*i+dim];
                            //gbp = i;
                        }


                    //printf("\ngb:%f\n",gb);

                    if(gb < goal && sum != 0)
                    {
                        results[2*vtest] = (double)vgen;
                        results[2*vtest+1] = gb;
                        sum = 0;
                        //ngoal++;
                        printf("\ngb(E):%e,gb:%1.30f,vgen:%d\t",gb,gb,vgen);
                        sucCounter+=1;
                    }
                    if(vgen == 1000)
                    {
                        stdindmemo[vtest]=gb;
                        printf("\nvtest:%d,1000gen's gb:%1.30f\n",vtest,gb);
                    }
                    if(vgen == gennum - 1)
                    {
                        printf("lastgb:%e\n",gb);
                        AveAllFit[vfun-1] += gb;

                        if(sum == 0) // Succeeded
                            AveEndFit[vfun-1] += gb;
                        else
                            results[2*vtest] = vgen;
                    }

                } // End of for VGEN!
            } // End of for VTEST   !

            //Calculating the Standard Index
            double avememofit=0;
            for(i=0; i<testnum; i++)
                avememofit+=stdindmemo[i];
            avememofit/=testnum;
            double delta=0;
            for(i=0; i<testnum; i++)
                delta += pow((stdindmemo[i]-avememofit),2);
             delta = sqrt(delta/testnum);

            double tempmemeo, median;
            for(i=0; i<testnum-1; i++)
                for(j=i+1; j<testnum; j++)
                    if(stdindmemo[i]>stdindmemo[j])
                    {
                        tempmemeo = stdindmemo[i];
                        stdindmemo[i] = stdindmemo[j];
                        stdindmemo[j] = tempmemeo;
                    }
            if(testnum%2 == 0)
            {
                median = (stdindmemo[testnum/2]+stdindmemo[testnum/2+1])/2;
                printf("\nmedian:%f;  %f,%f",median,stdindmemo[testnum/2],stdindmemo[testnum/2+1]);
            }
            else
                median = stdindmemo[testnum/2+1];


            stdindex[vlneighs] = (median - avememofit)/delta;

            printf("\nstdindex[vlneighs]%f",stdindex[vlneighs]);

            if(sucCounter!=0){
                Sto_SucEndFit[vlneighs]=AveEndFit[vfun-1]/sucCounter;
               // fprintf(pso,"\n\n\nLast succeeded average gb fitness:%1.30f\n",AveEndFit[vfun-1]/sucCounter);
            }
            else
                printf("\n\n\nLast succeeded average gb fitness: U gotta be kidding, man!\n");

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
            if(ngoal == 0) // Unsuccessful
                avggennew = (double)gennum - 1;
            else
                avggennew = (double)avggen / ngoal;
            fprintf(pso, "AveGen,%f,SucRate,%f,SucEndFit,%f,LastFit,%f,StdIndex,%f\n", avggennew, (double)sucCounter/testnum*100,Sto_SucEndFit[vlneighs],AveAllFit[vfun-1]/testnum,stdindex[vlneighs]);
            printf("\nAveGen,%f,SucRate,%f,SucEndFit,%f,LastFit,%f,StdIndex,%f\n", avggennew, (double)sucCounter/testnum*100,Sto_SucEndFit[vlneighs],AveAllFit[vfun-1]/testnum,stdindex[vlneighs]);

            Sto_SucGens[vlneighs]=avggennew;
            Sto_SucRate[vlneighs]=(double)sucCounter / testnum * 100;
        } // End of for VFUN!
    } // End of for vlneighs!

    fclose(pso);

    pso = fopen("lipso2.3-all-f1-SucEndFit.txt", "a+");
    for(i=0;i<N;i++)
        fprintf(pso,"%d,%1.30f\n",i+1,Sto_SucEndFit[i]);
    fclose(pso);
    pso = fopen("lipso2.3-all-f1-SucGens.txt", "a+");
    for(i=0;i<N;i++)
        fprintf(pso,"%d,%f\n",i+1,Sto_SucGens[i]);
    fclose(pso);
    pso = fopen("lipso2.3-all-f1-SucRate.txt", "a+");
    for(i=0;i<N;i++)
        fprintf(pso,"%d,%f\n",i+1,Sto_SucRate[i]);
    fclose(pso);
}
