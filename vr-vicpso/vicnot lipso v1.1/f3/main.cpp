/*
 FIPSO v1.1
 Van Ryan

 A particle itself within its neighborhood
 2013.12.10
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

int partion(double a[2][N],int p,int r)
{
    //rand
    time_t t;
    srand((unsigned) time(&t));
    int e=rand()%(r-p+1)+p;
    double tem;
    int tem2;
    tem=a[0][e];
    tem2=a[1][e];
    a[0][e]=a[0][r];
    a[1][e]=a[1][r];
    a[0][r]=tem;
    a[1][r]=tem2;
    int x=a[0][r], i=p-1;
    for (int j=p; j<r; j++)
    {
        if (a[0][j]<=x)
        {
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

void QuickSort(double a[2][N],int p,int r)
{
    if (p<r)
    {
        int q=partion(a,p,r);
        QuickSort(a,p,q-1);
        QuickSort(a,q+1,r);
    }
}

void BubbleSort(double arr[2][N],int length){
    double tmp;
    int i, j, tmp2;
    for (i = 0; i < length; i++) {
        for (j = length - 1; j > i; j--) {
            if (arr[0][j] < arr[0][j-1]) {
                tmp = arr[0][j-1];
                arr[0][j-1] =  arr[0][j];
                arr[0][j] = tmp;
                tmp2 = arr[1][j-1];
                arr[1][j-1] =  arr[1][j];
                arr[1][j] = tmp2;
            }
        }
    }
}
void get_vel_update(int lneighs)
{
    int i,j,k;
    double coef[N];
    double phi = 4.1;

    double PBNINDEX[2][N]; // Pbest & Index
    for(k=0; k<N; k++)
    {
        PBNINDEX[0][k] = pbest[(dim+1)*k+dim];  // pbest
        PBNINDEX[1][k] = k;  // The index of particle
    }
    //QuickSort(PBNINDEX,0,N-1);  // Done sorting the fitness while keeping the index
    // The quicksort gives a small to large sequence! (And fitness is better when smaller

    BubbleSort(PBNINDEX,N);
    // The bubblesort gives a small to large sequence! (And fitness is better when smaller

    //for(k=0; k<N; k++)
    //printf("%f,[%f],",PBNINDEX[0][k],PBNINDEX[1][k]);

    //PBNINDEX[1][0]=gbp;

    for(j = 0; j < N; j++)
    {
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
            // The particle should also learn from itself
            double temp = (double)rand() / RAND_MAX;
            sum_coef += temp;
            pos_ave += temp * pbest[(dim+1)*j+i];

            // After learning
            pos_ave /= sum_coef;
            double var_phi = sum_coef * phi / (lneighs+1);

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
        Sto_SucEndFit[i]=0;
        Sto_SucGens[i]=-1;
        Sto_SucRate[i]=-1;
    }

    pso = fopen("../../fipso1.1-all-f3.txt", "a+");
    //fprintf(pso, "num\tw\tc1\tc2\tfuncnum\tavggen\tratio\tfinalresult\n");

    x = 0.729;

    for(vlneighs =1; vlneighs< N; vlneighs++)
    {
        fprintf(pso, "lneighs:%d\n",vlneighs);
        for(i = 0; i < 4; i++){
            AveAllFit[i]=0;
            AveEndFit[i]=0;
        }

        for(vfun = 3; vfun < 4; vfun++)
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


                for(vgen = 0; vgen < gennum; vgen++)            //gennum代进化；
                {
                    /*
                    int k;
                    double phi = 4.1;
                    double PBNINDEX[2][N]; // Pbest & Index
                    for(k=0; k<N; k++)
                    {
                        PBNINDEX[0][k]=pbest[(dim+1)*k+dim];  // pbest
                        PBNINDEX[1][k]=k;  // The index of particle
                    }
                     for(k=0; k<N; k++){
                         printf("<%f,%d> ",PBNINDEX[0][k],(int)PBNINDEX[1][k]);
                     }
                     printf("\n");*/
                    //QuickSort(PBNINDEX,0,N-1);  // Done sorting the fitness while keeping the index
                    // The quicksort gives a small to large sequence! (And fitness is better when smaller
                    //for(k=0; k<N; k++)
                    //printf("%f,[%f],",PBNINDEX[0][k],PBNINDEX[1][k]);
                    /*for(k=0; k<N; k++){
                        printf("<%f,%d> ",PBNINDEX[0][k],(int)PBNINDEX[1][k]);
                    }*/
                    //printf("\n");
                    /*
                    for(j = 0; j < N; j++)
                    {
                        for(i = 0; i < dim; i++)
                        {
                            double r1=(double)rand() / RAND_MAX, r2=(double)rand() / RAND_MAX;
                            int kkk=(int)PBNINDEX[1][0];
                            update[(dim+1)*j+i] = phi*( r1*( pbest[(dim+1)*kkk+i]-coor[(dim+1)*j+i]) + r2*(pbest[(dim+1)*j+i]-coor[(dim+1)*j+i]) )/2;
                        }
                    }
                    */

                    get_vel_update(vlneighs);

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
                        printf("\ngb:%f,%d\t",gb,vgen);
                        sucCounter+=1;
                    }
                    if(vgen == 1000)
                    {
                        stdindmemo[vtest]=gb;
                        printf("\nvtest:%d,1000gen's gb:%1.30f\n",vtest,gb);
                    }
                    if(vgen == gennum - 1)
                    {
                        printf("lastgb:%f\n",gb);
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


            stdindex[vlneighs-1] = (median - avememofit)/delta;

            printf("\nstdindex[vlneighs-1]%f",stdindex[vlneighs-1]);

           if(sucCounter!=0){
                Sto_SucEndFit[vlneighs-1]=AveEndFit[vfun-1]/sucCounter;
                //fprintf(pso,"\n\n\nLast succeeded average gb fitness:%1.30f\n",AveEndFit[vfun-1]/sucCounter);
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

            if(sucCounter==0)
                Sto_SucEndFit[vlneighs-1]=-1;

            fprintf(pso, "AveGen,%f,SucRate,%f,SucEndFit,%1.30f,LastFit,%1.30f,StdIndex,%f\n", avggennew, (double)sucCounter/testnum*100,Sto_SucEndFit[vlneighs-1],AveAllFit[vfun-1]/testnum,stdindex[vlneighs-1]);
            printf("\nAveGen,%f,SucRate,%f,SucEndFit,%f,LastFit,%f,StdIndex,%f\n", avggennew, (double)sucCounter/testnum*100,Sto_SucEndFit[vlneighs-1],AveAllFit[vfun-1]/testnum,stdindex[vlneighs-1]);
            Sto_SucGens[vlneighs-1]=avggennew;
            Sto_SucRate[vlneighs-1]=(double)sucCounter / testnum * 100;
        } // End of for VFUN!
    } // End of for vlneighs!

    fclose(pso);

     pso = fopen("fipso1.1-all-f3-SucEndFit.txt", "a+");
    for(i=0;i<N;i++)
        fprintf(pso,"%d,%1.30f\n",i+1,Sto_SucEndFit[i]);
    fclose(pso);
    pso = fopen("fipso1.1-all-f3-SucGens.txt", "a+");
    for(i=0;i<N;i++)
        fprintf(pso,"%d,%f\n",i+1,Sto_SucGens[i]);
    fclose(pso);
    pso = fopen("fipso1.1-all-f3-SucRate.txt", "a+");
    for(i=0;i<N;i++)
        fprintf(pso,"%d,%f\n",i+1,Sto_SucRate[i]);
    fclose(pso);
}
