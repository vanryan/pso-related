/*
 FIPSO v1.1
 Van Ryan

 BA network
 2013.12.10

*/
#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

const int dim = 30;
const int N = 50;
const int testnum = 5;
const int gennum = 5000;
const int coornum = (dim + 1) * N;
int t[N][N]; //Adjacent Matrix

float coor[coornum];
FILE* pso;
FILE* bafp;

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

void BubbleSort(double arr[2][N],int length)
{
    double tmp;
    int i, j, tmp2;
    for (i = 0; i < length; i++)
    {
        for (j = length - 1; j > i; j--)
        {
            if (arr[0][j] < arr[0][j-1])
            {
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

void get_vel_update(float update[coornum], float pbest[coornum],int lneighs)
{
    int i,j,k;
    double coef[N];
    double phi = 4.1;

    double PBNINDEX[2][N]; // Pbest & Index
    for(k=0; k<N; k++)
    {
        PBNINDEX[0][k]=pbest[(dim+1)*k+dim];  // pbest
        PBNINDEX[1][k]=k;  // The index of particle
    }
    BubbleSort(PBNINDEX,N);  // Done sorting the fitness while keeping the index
    // The bubblesort gives a small to large sequence! (And fitness is better when smaller

    //for(k=0;k<N;k++)
     //          printf("<%f,%d>",PBNINDEX[0][k],(int)PBNINDEX[1][k]);
      //      printf("\n");

    for(j = 0; j < N; j++)
    {
        int con_neighs_counter = 0; // the counter of connected neighbors
        for(i = 0; i < N; i++)
            if(i!=j && t[j][i]!=INT_MAX)
                con_neighs_counter++;

        for(i = 0; i < dim; i++)
        {
            double sum_coef = 0;
            double pos_ave = 0;
            double var_phi;
            int lcountdown; // Max of neighbors you can learn from
            int llockdown;

            //Comparing connected_neighbors' number and lneighs
            if(con_neighs_counter>=lneighs)
                lcountdown=lneighs;
            else lcountdown=con_neighs_counter;

            llockdown=lcountdown;

            for(k=0;k<N;k++){
                if(lcountdown>0){
                    if(t[j][(int)PBNINDEX[1][k]]!=INT_MAX && (int)PBNINDEX[1][k]!=j){
                        coef[(int)PBNINDEX[1][k]] = (double)rand() / RAND_MAX;
                        sum_coef += coef[(int)PBNINDEX[1][k]];
                        pos_ave += coef[(int)PBNINDEX[1][k]] * pbest[(dim+1)*((int)PBNINDEX[1][k])+i];
                        lcountdown--; //Cut one off
                    }
                }
                else break;
            }


/*
            double temp = (double)rand() / RAND_MAX;
            sum_coef += temp;
            pos_ave += temp * pbest[(dim+1)*j+i];
            var_phi = sum_coef * phi / (llockdown-lcountdown+1);

            pos_ave /= sum_coef;

            double social_central_tendency = pos_ave - coor[(dim+1)*j+i];
            update[(dim+1)*j+i] = var_phi * social_central_tendency ;
*/

            for(k=0;k<N;k++){if(t[j][(int)PBNINDEX[1][k]]!=INT_MAX && (int)PBNINDEX[1][k]!=j) break;}
            double r1=(double)rand() / RAND_MAX,r2=(double)rand() / RAND_MAX;
            update[(dim+1)*j+i] = phi/2*(r1*(pbest[(dim+1)*j+i]- coor[(dim+1)*j+i]) +r2*(pbest[(dim+1)*((int)PBNINDEX[1][k])+i]- coor[(dim+1)*j+i]) );
        }
    }
}

void BAnet_produce()
{
    if((bafp = fopen("fipso-v1-ba-dis.txt", "rt")) == NULL)
    {
        printf("ERROR!");
        return ;
    }
    int i, j;
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
        {
            fscanf(bafp, "%d,", &t[i][j]);
            if(i == 0 && j == 0)
                t[i][j] = 0;
            if(t[i][j] == 2)
                t[i][j] = INT_MAX;

        }
    fclose(bafp);

}


int main()
{
    //printf("a");
    srand((unsigned)time(NULL));
    int i, j, l, vfun, vtest, vgen, vlneighs;
    float x, edge, Vmax, goal, gb, res = 0,avggennew=0;
    float vel[coornum], pbest[coornum], update[coornum], results[2*testnum];
    int sucCounter;

// vlneighs:The number of neighbors a paricle learns from (self-excluded)

    float AveEndFit[4];

    BAnet_produce(); // Produce the particle's BA distribution

    pso = fopen("fipso1.1-ba-f1-f4.txt", "a+");
    fprintf(pso, "num\tw\tc1\tc2\tfuncnum\tavggen\tratio\tfinalresult\n");

    x = 0.7298;

    for(vlneighs =1; vlneighs< 2; vlneighs++)
    {
        fprintf(pso, "lneighs:%d\n",vlneighs);
        for(i = 0; i < 4; i++)
            AveEndFit[i]=0;

        for(vfun = 1; vfun < 5; vfun++)
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

                for(vgen = 0; vgen < gennum; vgen++)            //gennum代进化；
                {
                    get_vel_update(update, pbest,vlneighs);

                    for(i = 0; i < dim + 1; i++)
                    {
                        for(j = 0; j < N; j++)
                        {
                            if(i < dim)
                            {
                                vel[(dim+1)*j+i] = x * ( vel[(dim+1)*j+i] + update[(dim+1)*j+i] );

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
                            gb = pbest[(dim+1)*i+dim];
                    if(gb < goal && sum != 0)
                    {
                        results[2*vtest] = (double)vgen;
                        results[2*vtest+1] = gb;
                        sum = 0;
                        //ngoal++;
                        printf("\ngb:%f,%d\t",gb,vgen);
                        sucCounter+=1;
                    }
                    if(vgen == gennum - 1)
                    {
                        printf("lastgb:%f\n",gb);

                        if(sum == 0) // Succeeded
                            AveEndFit[vfun-1] += gb;
                        else
                            results[2*vtest] = vgen;
                    }

                } // End of for VGEN!
            } // End of for VTEST   !

            if(sucCounter!=0)
                fprintf(pso,"\n\n\nLast succeeded average gb fitness:%f\n",AveEndFit[vfun-1]/sucCounter);
            else
                fprintf(pso,"\n\n\nLast succeeded average gb fitness: U gotta be kidding, man!\n");

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
            fprintf(pso, "\n平均代数为：%f\n达优率为：%f%%", avggennew, (double)sucCounter / testnum * 100);
        } // End of for VFUN!
    } // End of for vlneighs!

    fclose(pso);
}
