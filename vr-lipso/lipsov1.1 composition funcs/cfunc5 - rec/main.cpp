/*
 FIPSO v1.1
 Van Ryan

 A particle itself within its neighborhood
 2013.12.10

 Using composite functions
 2014.01.20
*/
#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

#include"ortho_rot_matrix.h"
#include"basicfuncs.h"
#include"hybridfuncs.h"

const int N = 50;
const int testnum = 20;
const int gennum = 5000;
const int coornum = (dim + 1) * N;

double pbest[coornum];
double update[coornum];

double coor[coornum];

FILE* pso;

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
    float vel[coornum], results[2*testnum];
    int sucCounter;
    double fmax[basic_func_num], or_matrix[dim][dim], opt_arr[basic_func_num][dim], opt_arr_old[basic_func_num][dim];

// vlneighs:The number of neighbors a paricle learns from (self-excluded)
// AveEndFit stores the ending fitnesses of only the successful tests, meanwhile AveAllEndFit stores ending fitnesses regardless of success or failure
    float AveEndFit[4],AveAllEndFit[4];

    hyb_func5_init(fmax,or_matrix,opt_arr,opt_arr_old);

    pso = fopen("lipso1.1-rec-all-cf5.txt", "a+");
    fprintf(pso,">>>TEST RECORD START<<<\n");
    //Recording random stuff: fmax,or_matrix,opt_arr
    //fmax
    fprintf(pso,"\nfmax:\n");
    for(i=0;i<basic_func_num;i++){
        fprintf(pso,"fmax,%f\n",fmax[i]);
    }
    //opt_arr
    fprintf(pso,"\nopt_arr:\n");
    for(i=0;i<basic_func_num;i++){
        fprintf(pso,"f%d\n",i+1);
        for(j=0;j<dim;j++)
            fprintf(pso,"opt_arr,%f\n",opt_arr[i][j]);
        fprintf(pso,"\n");
    }
    fprintf(pso,"\nD*D rotate orthodox matrix:\n");
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++)
            fprintf(pso,"%f,",or_matrix[i][j]);
        fprintf(pso,"\n");
    }
    fprintf(pso,">>>TEST RECORD OVER<<<\n");


    fclose(pso);

    pso = fopen("lipso1.1-all-cf5.txt", "a+");

    x = 0.7298;

    for(vlneighs =1; vlneighs< N; vlneighs++)
    {
        for(i = 0; i < 4; i++){
            AveAllEndFit[i]=0;
            AveEndFit[i]=0;
        }

        for(vfun = 5; vfun < 6; vfun++)
        {
            for(i = 0; i < testnum; i++)
            {
                results[2*i] = gennum;
                results[2*i+1] = 0;
            }

            edge = hyb_func_goalnedgenvmax[vfun-1][1];
            Vmax = hyb_func_goalnedgenvmax[vfun-1][2];
            goal = hyb_func_goalnedgenvmax[vfun-1][0];

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
                        printf("\nReached gb:%f,at %d\t",gb,vgen);
                        sucCounter+=1;
                    }
                    if(vgen == gennum - 1)
                    {
                        printf("lastgb:%f\n",gb);
                        AveAllEndFit[vfun-1] += gb;
                        if(sum == 0) // Succeeded
                            AveEndFit[vfun-1] += gb;
                        else
                            results[2*vtest] = vgen;
                    }

                } // End of for VGEN!
            } // End of for VTEST   !

            if(sucCounter!=0)
                printf("\n\n\nLast succeeded average gb fitness:%f\nLast average gb fitness:%f\n",AveEndFit[vfun-1]/sucCounter,AveAllEndFit[vfun-1]/testnum);
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
            if(ngoal == 0)
                avggennew = (double)gennum - 1;
            else
                avggennew = (double)avggen / ngoal;

            fprintf(pso, "平均代数为,%f,达优率,%f%%,LastSucGb,%f,LastAllGb,%f\n", avggennew, (double)sucCounter / testnum * 100, AveEndFit[vfun-1]/sucCounter,AveAllEndFit[vfun-1]/testnum);

        } // End of for VFUN!
    } // End of for vlneighs!

    fclose(pso);
}
