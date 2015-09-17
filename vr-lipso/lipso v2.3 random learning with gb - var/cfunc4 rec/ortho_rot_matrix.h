#include <iostream>
#include"math.h"
#include"time.h"
#include"limits.h"
#include"stdlib.h"

using namespace std;

const int DIM_NUM=10;
const double PI=3.1415926;

void diag_emat(double me[DIM_NUM][DIM_NUM])
{
    int i,j;
    for(i=0; i<DIM_NUM; i++)
        for(j=0; j<DIM_NUM; j++)
        {
            if(i==j)
                me[i][j]=1;
            else
                me[i][j]=0;
        }

}
void add_mat(double ma[DIM_NUM][DIM_NUM],double mb[DIM_NUM][DIM_NUM]){
    int i,j;
    for(i=0;i<DIM_NUM;i++)
        for(j=0; j<DIM_NUM; j++)
            ma[i][j]+=mb[i][j];
}
void minus_mat(double ma[DIM_NUM][DIM_NUM],double mb[DIM_NUM][DIM_NUM]){
    int i,j;
    for(i=0;i<DIM_NUM;i++)
        for(j=0; j<DIM_NUM; j++)
            ma[i][j]-=mb[i][j];
}

void multp_mat(double ma[DIM_NUM][DIM_NUM],double mb[DIM_NUM][DIM_NUM])
{
    int i,j,k;
    double mresult[DIM_NUM][DIM_NUM];
    for(i=0;i<DIM_NUM;i++)
        for(j=0; j<DIM_NUM; j++)
            mresult[i][j]=0;
    for(i=0; i<DIM_NUM; i++)
        for(j=0; j<DIM_NUM; j++)
            for(k=0; k<DIM_NUM; k++)
               mresult[i][j]+=ma[i][k]*mb[k][j];
    for(i=0; i<DIM_NUM; i++)
        for(j=0; j<DIM_NUM; j++)
            ma[i][j]=mresult[i][j];
}
void multp_mat2(double ma[DIM_NUM],double mb[DIM_NUM][DIM_NUM]){
    int i,k;
    double mresult[DIM_NUM];
    for(i=0;i<DIM_NUM;i++)
        mresult[i]=0;
    for(i=0; i<DIM_NUM; i++)
        for(k=0; k<DIM_NUM; k++)
            mresult[i]+=ma[k]*mb[k][i];
    for(i=0; i<DIM_NUM; i++)
        ma[i]=mresult[i];
}

void rot_mat_create(double mturn[DIM_NUM][DIM_NUM],int i,int j) // Create a single rotation matrix
{
    double alpha;

    diag_emat(mturn);

    alpha=((double)rand() / RAND_MAX -0.5)*PI*0.5;

    mturn[i][i]=cos(alpha);
    mturn[j][j]=cos(alpha);
    mturn[i][j]=sin(alpha);
    mturn[j][i]=-sin(alpha);
}
void ort_rot_mat_create(double mresult[DIM_NUM][DIM_NUM])
{
    diag_emat(mresult);//initialize mresult
    srand((unsigned)time(NULL));

    int i;
    double mturn[DIM_NUM][DIM_NUM];
    for(i=1;i<DIM_NUM;i++){
        rot_mat_create(mturn,0,i);
        multp_mat(mresult,mturn);
    }
    for(i=1;i<DIM_NUM-1;i++){
        rot_mat_create(mturn,i,DIM_NUM-1);
        multp_mat(mresult,mturn);
    }
}
