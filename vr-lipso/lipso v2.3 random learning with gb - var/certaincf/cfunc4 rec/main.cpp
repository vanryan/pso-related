/*
 FIPSO v2.3
 Van Ryan

 A particle learns from random neighbors ( the number of the learnees is set [we still use vlneighs] )

 Using composite functions
 2014.02.20
*/
#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"
#include<fstream>
using namespace std;

#include"ortho_rot_matrix.h"
#include"basicfuncs.h"
#include"hybridfuncs.h"

const int N = 50;
const int testnum = 20;
const int gennum = 5000;
const int coornum = (dim + 1) * N;

double pbest[coornum];
double update[coornum];

int randlneighs[N];

double coor[coornum], stdindmemo[testnum];

FILE* pso;

void choose_rand_neighbors(int lneighs)
{
    int counter=0;
    int i,flag=1,inttemp;
    double tmp;
    for(i=0; i<N; i++)
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
    double fmax[basic_func_num], or_matrix[dim][dim], opt_arr[basic_func_num][dim], opt_arr_old[basic_func_num][dim];

// vlneighs:The number of neighbors a paricle learns from (self-excluded)
// AveEndFit stores the ending fitnesses of only the successful tests, meanwhile AveAllEndFit stores ending fitnesses regardless of success or failure
    float AveEndFit[4],AveAllEndFit[4];

    hyb_func4_init(fmax,or_matrix,opt_arr,opt_arr_old);

    pso = fopen("lipso2.3-rec-all-cf4.txt", "a+");
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

    pso = fopen("lipso2.3-all-cf4.txt", "a+");

    x = 0.7298;

    for(vlneighs =0; vlneighs< N-1; vlneighs++)
    {
        for(i = 0; i < 4; i++){
            AveAllEndFit[i]=0;
            AveEndFit[i]=0;
        }

        for(vfun = 4; vfun < 5; vfun++)
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
                    if(vgen == 1000)
                    {
                        stdindmemo[vtest]=gb;
                        printf("\nvtest:%d,1000gen's gb:%1.30f\n",vtest,gb);
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

            fprintf(pso, "平均代数为,%f,达优率,%f%%,LastSucGb,%f,LastAllGb,%f,StdIndex,%f\n", avggennew, (double)sucCounter / testnum * 100, AveEndFit[vfun-1]/sucCounter,AveAllEndFit[vfun-1]/testnum, stdindex[vlneighs]);
        } // End of for VFUN!
    } // End of for vlneighs!

    fclose(pso);
}
