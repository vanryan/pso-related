#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"

const int dim = 30;
const int N = 50;
const int testnum = 1000;
const int gennum = 5000;
const int coornum = (dim + 1) * N;

FILE* pso;
FILE* fp;

int t[N][N];
int Ft[N][N];
int Ct[N][N];
int Dt[N][N];
int pointdeg[N], gbpoint[testnum][3*(gennum+1)];
double vel[gennum+1][coornum], pbest[gennum+1][coornum], gbest[gennum+1][coornum], coor[gennum+1][coornum],edgert[gennum+1][2500];

void Fun1(int g, int a)                             //Sphere
{
    coor[g][(dim+1)*a+dim] = 0;
    int i;
    for(i = 0; i < dim; i++)
        coor[g][(dim+1)*a+dim] += pow(coor[g][(dim+1)*a+i], 2);
}

void Fun2(int g, int a)                            //Rosenbrock
{
    coor[g][(dim+1)*a+dim] = 0;
    int i;
    for(i = 0; i < dim; i++)
    {
        if(i < dim - 1)
            coor[g][(dim+1)*a+dim] += 100 * pow(( coor[g][(dim+1)*a+i+1]  -  pow(coor[g][(dim+1)*a+i] , 2) ) , 2)  + pow(coor[g][(dim+1)*a+i] - 1 , 2);
        else
            coor[g][(dim+1)*a+dim] += 0;//100*pow(( coor[(dim+1)*a]  -  pow(coor[(dim+1)*a+i] ,2) ) ,2)  + pow(coor[(dim+1)*a+i] - 1 ,2);
    }

}

void Fun3(int g, int a)                            //Rastrigrin
{
    coor[g][(dim+1)*a+dim] = 0;
    int i;
    for(i = 0; i < dim; i++)
        coor[g][(dim+1)*a+dim] += pow(coor[g][(dim+1)*a+i], 2) - 10 * cos(2 * 3.14 * coor[g][(dim+1)*a+i]) + 10;
}

void Fun4(int g, int a)                          //Griewank
{
    coor[g][(dim+1)*a+dim] = 0;
    int i;
    double b = 0, c = 1;
    for(i = 0; i < dim; i++)
    {
        b +=  pow(coor[g][(dim+1)*a+i]  , 2) / 4000;
        c  *=  cos(coor[g][(dim+1)*a+i]  / sqrt(i + 1));
    }
    coor[g][(dim+1)*a+dim] = 1 + b - c;
}

void Distribution()
{
    pso = fopen("7-1p=6distribution.txt", "a+");
    int degree[N], dis[N];
    int i, j, l;
    for(i = 0; i < N; i++)
    {
        degree[i] = 0;
        for(j = 0; j < N; j++)
            if(Dt[i][j] == 1)
                degree[i] += 1;
        //fprintf(divide,"%d : %d\n",i,degree[i]);
    }

    for(i = 0; i < N; i++)
    {
        pointdeg[i] = degree[i];
        fprintf(pso, "pointdeg[%d]:%d\n", i, pointdeg[i]);
    }

    for(l = 0; l < N; l++)
    {
        dis[l] = 0;
        for(i = 0; i < N; i++)
            if(degree[i] == l)
                dis[l] += 1;
        fprintf(pso, "%d:%d\n", l, dis[l]);
    }
    fclose(pso);
}

void Banet(int a)
{

    if((fp = fopen("a.txt", "rt")) == NULL)
    {
        printf("ERROR!");
        return ;
    }
    int i, j;
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
        {
            fscanf(fp, "%d,", &t[i][j]);
            if(i == 0 && j == 0)
                t[i][j] = 0;
            if(t[i][j] == 2)
                t[i][j] = INT_MAX;

        }
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
        {
            //printf("%d;", t[i][j]);
            Ft[i][j] = t[i][j];
            Ct[i][j] = t[i][j];
            Dt[i][j] = t[i][j];
        }
    fclose(fp);
    Distribution();
}


int main()
{
    srand((unsigned)time(NULL));
    int i, j, k, l, varw, varc1, varc2, vfun, vtest, vgen;
    double w, c1, c2, edge, Vmax, goal, gb, res = 0, avggennew = 0;
    double results[4*testnum];
    Banet(2);

    pso = fopen("sf-edge-fun1-dsinout.txt", "a+");
    //fprintf(pso, "num\tw\tc1\tc2\tfuncnum\tavggen\tratio\tfinalresult\n");
    /*for(i = 0; i < N; i++)
    {
        for(j = 0; j < N; j++)
            printf("%d,", t[i][j]);
        printf("\n");
    }*/

    w=0.6;
    for(varw=0; varw<1; varw++)
    {
        c1 = 1.7;
        for(varc1 = 0; varc1 < 1; varc1++)
        {
            c2 = 1.7;
            for(varc2 = 0; varc2 < 1; varc2++)
            {
                //fprintf(pso,"\n\nw=%f,c1=%f,c2=%f",w,c1,c2);
                for(vfun = 1; vfun < 2; vfun++)
                {
                    for(i=0; i<gennum+1; i++)
                        for(j=0; j<N*N; j++)
                            edgert[i][j]=0;
                    for(i = 0; i < testnum; i++)
                    {
                        results[4*i] = gennum;
                        results[4*i+1] = -1;
                        results[4*i] = gennum;
                        results[4*i+1] = -1;
                    }
                    int gbdeg[3*N];
                    for(i = 0; i < 3 * N; i++)
                    {
                        gbdeg[i] = 0;
                    }
                    res = 0;
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
                        edge = 500;
                        Vmax = edge;
                        goal = 0.05;
                    }
                    int ngoal = 0, avggen = 0;

                    for(vtest = 0; vtest < testnum; vtest++)                 //pso开始；
                    {
                        printf("a");
                        int sum = 1;
                        /*for(i = 0; i < gennum + 1; i++)                          //度记录数组gbpoint[][]初始化；
                        {
                            gbpoint[vtest][3*i] = N + 1;
                            gbpoint[vtest][3*i+1] = 0;
                            gbpoint[vtest][3*i+2] = INT_MAX;
                        }*/

                        for(i = 0; i < gennum + 1; i++)                      //粒子历史最优、邻域最优数组初始化；
                            for(j = 0; j < coornum; j++)
                            {
                                pbest[i][j] = INT_MAX;
                                gbest[i][j] = INT_MAX;
                            }

                        for(i = 0; i < dim + 1; i++)               //粒子初始化；
                        {
                            for(j = 0; j < N; j++)
                            {
                                if(i < dim)
                                {
                                    coor[0][(dim+1)*j+i] = (double)rand() / RAND_MAX * 2 * edge - edge;
                                    vel[0][(dim+1)*j+i] = (double)rand() / RAND_MAX * 2 * Vmax - Vmax;
                                    pbest[0][(dim+1)*j+i] = coor[0][(dim+1)*j+i];
                                    //printf("%f\t", coor[(dim+1)*j+i]);
                                }
                                else
                                {
                                    if(vfun == 1)
                                        Fun1(0, j);
                                    if(vfun == 2)
                                        Fun2(0, j);
                                    if(vfun == 3)
                                        Fun3(0, j);
                                    if(vfun == 4)
                                        Fun4(0, j);
                                    //printf("\n%d:%f,", j, coor[(dim+1)*j+i]);
                                    pbest[0][(dim+1)*j+i] = coor[0][(dim+1)*j+i];
                                    if(coor[0][(dim+1)*j+i] < gbest[0][(dim+1)*j+i])
                                    {
                                        for(l = 0; l < dim + 1; l++)
                                            gbest[0][(dim+1)*j+l] = coor[0][(dim+1)*j+l];
                                        //printf("gbest:%f\t", gbest[(dim+1)*j+i]);
                                    }

                                }
                            }

                        }

                        int edgetemp[N];
                        for(j = 0; j < N; j++)
                        {
                            edgetemp[j]=j;
                            for(k = 0; k < N; k++)
                                if( t[j][k] == 1 && pbest[0][(dim+1)*k+dim] < gbest[0][(dim+1)*j+dim])
                                {
                                    //printf("t[%d][%d],",j,k);
                                    for(l = 0; l < dim + 1; l++)
                                        gbest[0][(dim+1)*j+l] = pbest[0][(dim+1)*k+l];
                                    edgetemp[j] = k;
                                }
                            edgert[0][N*edgetemp[j]+j]+=1;
                            //printf("%d:%d:%d:%f\n",j,edgetemp[j],50*edgetemp[j]+j,edgert[0][50*edgetemp[j]+j]);
                        }

                        gb = gbest[0][dim];
                        for(i=1; i<N; i++)
                            if(gbest[0][(dim+1)*i+dim]<gb)
                                gb=gbest[0][(dim+1)*i+dim];
                        //fprintf(pso,"%d:%f\n",0,gb);

                        /*gbpoint[vtest][0] = 0;
                        gbpoint[vtest][1] = pointdeg[0];
                        gbpoint[vtest][2] = pbest[0][dim];
                        //float x =  pbest[0][dim];
                        for(i = 1; i < N; i++)
                            if(pbest[0][(dim+1)*i+dim] < gbpoint[vtest][2])
                            {
                                gbpoint[vtest][0] = i;
                                gbpoint[vtest][3] = i;
                                gbpoint[vtest][1] = pointdeg[i];
                                gbpoint[vtest][4] = pointdeg[i];
                                gbpoint[vtest][2] = pbest[0][(dim+1)*i+dim] ;
                                gbpoint[vtest][5] = pbest[0][(dim+1)*i+dim] ;
                            }*/
                        //printf("0:%f\n",x);

                        /*for(i = 0; i < coornum; i++)
                        {
                            printf("\nabc:%f,%f,%f", coor[i], pbest[i], gbest[i]);
                            if(i % (dim + 1) == dim)
                                printf("\n");
                        }*/
                        for(vgen = 0; vgen < gennum; vgen++)            //gennum代进化；
                        {
                            //w = 1 - 0.55 * (double)vgen / (gennum - 1);
                            for(i = 0; i < dim + 1; i++)
                            {
                                for(j = 0; j < N; j++)
                                {
                                    if(i < dim)
                                    {
                                        double r1 = (double)rand() / RAND_MAX * 1;
                                        double r2 = (double)rand() / RAND_MAX * 1;
                                        vel[vgen+1][(dim+1)*j+i] = w * vel[vgen][(dim+1)*j+i] + c1 * r1 * (pbest[vgen][(dim+1)*j+i] - coor[vgen][(dim+1)*j+i])\
                                                                   + c2 * r2 * (gbest[vgen][(dim+1)*j+i] - coor[vgen][(dim+1)*j+i]);
                                        /*if(vel[vgen+1][(dim+1)*j+i] > Vmax)
                                            vel[vgen+1][(dim+1)*j+i] = Vmax;
                                        if(vel[vgen+1][(dim+1)*j+i] < -Vmax)
                                            vel[vgen+1][(dim+1)*j+i] = -Vmax;*/

                                        coor[vgen+1][(dim+1)*j+i] = coor[vgen][(dim+1)*j+i] + vel[vgen+1][(dim+1)*j+i];
                                        /*if(coor[vgen+1][(dim+1)*j+i] > edge)
                                        {
                                            coor[vgen+1][(dim+1)*j+i] = edge;
                                            vel[vgen+1][(dim+1)*j+i] = 0;
                                        }

                                        if(coor[vgen+1][(dim+1)*j+i] < -edge)
                                        {
                                            coor[vgen+1][(dim+1)*j+i] = -edge;
                                            vel[vgen+1][(dim+1)*j+i] = 0;
                                        }*/
                                    }
                                    else
                                    {
                                        if(vfun == 1)
                                            Fun1(vgen + 1, j);
                                        if(vfun == 2)
                                            Fun2(vgen + 1, j);
                                        if(vfun == 3)
                                            Fun3(vgen + 1, j);
                                        if(vfun == 4)
                                            Fun4(vgen + 1, j);

                                        if(coor[vgen+1][(dim+1)*j+i] < pbest[vgen][(dim+1)*j+i])            //pbest更新；
                                            for(l = 0; l < dim + 1; l++)
                                                pbest[vgen+1][(dim+1)*j+l] = coor[vgen+1][(dim+1)*j+l];
                                        else
                                            for(l = 0; l < dim + 1; l++)
                                                pbest[vgen+1][(dim+1)*j+l] = pbest[vgen][(dim+1)*j+l];

                                        if(coor[vgen+1][(dim+1)*j+i] < gbest[vgen][(dim+1)*j+i])              //gbest更新；
                                        {
                                            edgetemp[j]=j;
                                            for(l = 0; l < dim + 1; l++)
                                                gbest[vgen+1][(dim+1)*j+l] = coor[vgen+1][(dim+1)*j+l];
                                        }
                                        else
                                            for(l = 0; l < dim + 1; l++)
                                                gbest[vgen+1][(dim+1)*j+l] = gbest[vgen][(dim+1)*j+l];


                                    }
                                }
                            }
                            for(j = 0; j < N; j++)
                            {
                                for(k = 0; k < N; k++)
                                    if( t[j][k] == 1 )
                                    {
                                        if(pbest[vgen+1][(dim+1)*k+dim] < gbest[vgen+1][(dim+1)*j+dim])
                                        {
                                            edgetemp[j]=k;
                                            for(l = 0; l < dim + 1; l++)
                                                gbest[vgen+1][(dim+1)*j+l] = pbest[vgen+1][(dim+1)*k+l];
                                        }
                                        /*if(gbest[vgen+1][(dim+1)*j+i] > gbest[vgen+1][(dim+1)*k+i])
                                            for(l = 0; l < dim + 1; l++)
                                                gbest[vgen+1][(dim+1)*j+l] = gbest[vgen+1][(dim+1)*k+l];*/
                                    }

                                edgert[vgen+1][N*edgetemp[j]+j]+=1;
                            }

                            /* if(vgen < gennum - 1)
                             {
                                 for(i = 0; i < N; i++)
                                     if(pbest[vgen+1][(dim+1)*i+dim]  < gbpoint[vtest][(vgen+1)*3+2] )
                                     {
                                         gbpoint[vtest][(vgen+1)*3] = i;

                                         gbpoint[vtest][(vgen+1)*3+1] = pointdeg[i];

                                         gbpoint[vtest][(vgen+1)*3+2] = pbest[vgen+1][(dim+1)*i+dim] ;
                                     }
                                 gbpoint[vtest][(vgen+1)*3+3] = gbpoint[vtest][(vgen+1)*3];
                                 gbpoint[vtest][(vgen+1)*3+4] = gbpoint[vtest][(vgen+1)*3+1];
                                 gbpoint[vtest][(vgen+1)*3+5] = gbpoint[vtest][(vgen+1)*3+2];
                                 //if(vg)
                             }

                             if(vgen == gennum - 1)
                             {
                                 printf("a");//printf("%d_\n",gbpoint[vtest][(vgen+1)*2]);
                                 for(i = 0; i < N; i++)
                                     if(pbest[vgen+1][(dim+1)*i+dim]  < gbpoint[vtest][(vgen+1)*3+2])
                                     {
                                         gbpoint[vtest][(vgen+1)*3] = i;//printf("%d_",gbpoint[vtest][(vgen+1)*2]);
                                         //gbpoint[vtest][(vgen+1)*2+2] = i;
                                         gbpoint[vtest][(vgen+1)*3+1] = pointdeg[i];
                                         //gbpoint[vtest][(vgen+1)*2+3] = pointdeg[i];
                                         gbpoint[vtest][(vgen+1)*3+2] = pbest[vgen+1][(dim+1)*i+dim] ;
                                     }
                             }*/


                            gb = gbest[vgen+1][dim];
                            for(i = 1; i < N; i++)
                                if(gbest[vgen+1][(dim+1)*i+dim] < gb)
                                    gb = gbest[vgen+1][(dim+1)*i+dim];
                            //fprintf(pso,"%d:%f\n",vgen+1,gb);
                            if(gb < goal && sum != 0)                                //达标点取得
                            {
                                results[4*vtest] = (double)vgen;
                                results[4*vtest+1] = gb;
                                sum = 0;
                                //ngoal++;
                                //break;
                            }
                            if(vgen == gennum - 1)
                            {
                                results[4*vtest+2] = gennum - 1;
                                results[4*vtest+3] = gb;
                            }
                        }                       //end of gennum;
                    }                           //end of testnum;

                    double sin[2*N];
                    double sout[2*N];
                    double eravg[196];
                    for(i=0; i<196; i++)
                        eravg[i]=0;
                    for(i=0; i<2*N; i++)
                    {
                        sin[i]=0;
                        sout[i]=0;
                    }
                    int eravgtemp=0;
                    for(i=0; i<N; i++)
                        for(j=0; j<N; j++)
                        {
                            if(t[i][j]==1)
                                fprintf(pso,",%d->%d",i,j);
                        }
                    for(i=0; i<gennum+1; i++)
                    {
                        eravgtemp =0;
                        fprintf(pso,"\n%d,",i);
                        for(j=0; j<N*N; j++)
                        {
                            int le=j/N;
                            int ri=j%N;
                            edgert[i][j] /= testnum;
                            if(t[le][ri]==1)
                            {
                                sin[2*pointdeg[ri]]+=edgert[i][j];
                                sin[2*pointdeg[ri]+1]+=1;
                                sout[2*pointdeg[le]]+=edgert[i][j];
                                sout[2*pointdeg[le]+1]+=1;
                                eravg[eravgtemp] +=edgert[i][j];
                                fprintf(pso,"%.16f,",edgert[i][j]);
                                // printf("%d,%d,%d,%f\n",le,ri,t[le][ri],edgert[i][j]);
                                eravgtemp++;
                            }
                        }

                    }
                    fprintf(pso,"\n");
                    for(i=0; i<196; i++)
                    {
                        eravg[i] /=5001;
                        fprintf(pso,",%.16f",eravg[i]);
                    }
                    fprintf(pso,"\n");
                    for(i=0; i<100; i++)
                    {
                        int ps=0;
                        for(j=0; j<196; j++)
                        {
                            if(eravg[j]>(double)i*0.01)
                            {
                                ps++;
                            }
                        }
                        fprintf(pso,"边权为,%.16f,%d\n",(double)i*0.01,ps);
                    }
                    for(i=0; i<N; i++)
                    {
                        if(sin[2*i+1]!=0)
                            fprintf(pso,"point right degree,%d,in,%.16f,",i,sin[2*i]/sin[2*i+1]);
                        if(sout[2*i+1]!=0)
                            fprintf(pso,"out,%.16f\n",sout[2*i]/sout[2*i+1]);
                    }

                    double iaoedgert[2500];               //p(s)-s-deg distribution;
                    for(j=0; j<2500; j++)
                        iaoedgert[j]=0;
                    for(i=0; i<gennum+1; i++)
                    {
                        for(j=0; j<2500; j++)
                        {
                            iaoedgert[j] +=edgert[i][j];
                        }
                    }
                    for(i=0; i<2500; i++)
                        iaoedgert[i] /= (gennum+1);

                    for(i=0; i<N; i++)
                    {
                        fprintf(pso,"%d",i);
                        for(k=0; k<100; k++)
                        {
                            int dsin=0;
                            int dsout=0;
                            for(j=0; j<2500; j++)
                            {
                                int le = j/N;
                                int ri = j%N;
                                if( pointdeg[le]==i && t[le][ri] ==1 && iaoedgert[j]> k*0.01)
                                {
                                    dsout++;
                                }
                                if(pointdeg[ri]==i && t[le][ri]==1 && iaoedgert[j]> k*0.01)
                                {
                                    dsin++;
                                }
                            }
                            fprintf(pso,",%.16f,dsin,%d,dsout,%d\n",(double)k*0.01,dsin,dsout);
                        }
                    }



                    /*printf("%f,\n", res);
                    for(i = 0; i < testnum; i++)
                    {
                        if(results[2*i+1] < goal) //avggen += results[2*i];
                        {
                            res += results[2*i+1];
                            ngoal++;
                        }
                    }
                    if(ngoal == 0)
                        printf("\nngoal=0!!\n");
                    else
                    {
                        printf("AAA");
                        printf("\nngoal=%d", ngoal);
                        res = res / ngoal;
                        printf("res:%f,\n", res);
                    }*/

                    //printf("%d,%d", avggen, ngoal);
                    /*if(ngoal == 0)
                        avggennew = (double)gennum - 1;
                    else
                        avggennew = (double)avggen / ngoal;*/
                    //fprintf(pso, "\nasdfasdfasdfafdasdfs\n平均代数为：%f\n达优率为：%f%%", avggennew, (double)ngoal / testnum * 100);

                    //fprintf(pso, "c1=%f,c2=%f,vfun=%d\n", c1, c2, vfun);

                    /*for(i = 0; i < testnum; i++)
                    {
                        int x = N;
                        for(k = 0; k < gennum; k++)
                        {
                            if(gbpoint[i][3*k+2] < goal)
                            {
                                //gbdeg[3*j] = k;
                                x = gbpoint[i][3*k+1];
                                gbdeg[3*x+1] += 1;
                                break;
                            }
                        }
                    }
                    for(j = 0; j < N; j++)
                    {
                        for(i = 0; i < testnum; i++)
                        {
                            //for(j=0;j<gennum+1;j++)
                             //fprintf(pso,"NO.%d gen:  num_%d, degree_%d\t",j,gbpoint[i][2*j],gbpoint[i][2*j+1]);
                            //fprintf(pso, "testnum=%d, num_%d, degree_%d", i, , gbpoint[i][2*gennum+1]);
                            //fprintf(pso, "\n");

                            if(gbpoint[i][3*gennum+1] == j)
                                gbdeg[3*j+2] += 1;

                        }
                        fprintf(pso, "degree_%d_num1_%d_num2_%d\n", j,gbdeg[3*j+1], gbdeg[3*j+2]);
                    }*/
                    double zhi=0;
                    for(i = 0; i < testnum; i++)
                    {
                        if(results[4*i] != gennum)
                        {
                            ngoal++;
                            avggen += results[4*i];
                            res += results[4*i+1];
                            zhi +=results[4*i+3];
                        }
                    }
                    printf("%d,%d", avggen, ngoal);
                    if(ngoal == 0)
                    {
                        avggennew = (double)gennum - 1;
                        zhi=-1;
                    }
                    else
                    {
                        avggennew = (double)avggen / ngoal;
                        zhi =zhi/ngoal;
                    }
                    fprintf(pso, "平均代数为:%.16f:达优率为:%.16f%%:最优值:%.16f", avggennew, (double)ngoal / testnum * 100,zhi);

                }             //end of vfun1-4;
                c2 += 0.1;
            }                //end of c2;
            c1 += 0.1;
        }                      //end of c1;
        w+=0.03;

    }
    fclose(pso);
}
