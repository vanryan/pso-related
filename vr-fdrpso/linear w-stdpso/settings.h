/*
 Parameters for PSO
*/
const double wmax = 0.9;
const double wmin = 0.4;
const double c1 = 1;
const double c2 = 1;
const double c3 = 2;

const int dim = 20;
const int N = 10;
const int testnum = 20;
const int gennum = 1000;
const int coornum = (dim + 1) * N;

double coor[coornum];
double vel[coornum], pbest[coornum], gbest[coornum];

FILE* pso;

/*
 Useful Functions
*/
double getMax(double a, double b)
{
    if(a>b)
        return a;
    else
        return b;
}

double getMin(double a, double b)
{
    if(a<b)
        return a;
    else
        return b;
}

double countFIndex(double record[testnum]) //Getting the standardized index
{
    double findex;
    double avememo=0,delta=0;
    int i,j;
    for(i=0; i<testnum; i++)
        avememo+=record[i];
    avememo /= testnum;
    for(i=0; i<testnum; i++) // Calculating Standard Deviation
        delta += pow((record[i]-avememo),2);
    delta = sqrt(delta/testnum);

    double tempmemeo, median;
    for(i=0; i<testnum-1; i++)
        for(j=i+1; j<testnum; j++)
            if(record[i]>record[j])
            {
                tempmemeo = record[i];
                record[i] = record[j];
                record[j] = tempmemeo;
            }
    if(testnum%2 == 0)
        median = (record[testnum/2]+record[testnum/2+1])/2;
    else
        median = record[testnum/2+1];

    findex=(median - avememo)/delta;

    return findex;
}

int findNeigh_fdr(int a, int d) // Find the neighbor to learn from based on
{
    /*
        @params
        int a: the number of the particle
        int d: the number of the dimension
        @End of params
    */
    int i,tmpmaxid,init_flag;
    double tmpmax,fdr[N][2]; // The second dimension is for both storing values and flags

    for(i=0; i<N; i++)
        if(fabs(pbest[(dim+1)*i+d]-coor[(dim+1)*a+d])!=0)
        {
            fdr[i][0]=1;
            fdr[i][1]=(coor[(dim+1)*a+dim]-pbest[(dim+1)*i+dim])/fabs(pbest[(dim+1)*i+d]-coor[(dim+1)*a+d]);
        }
        else
        {
            fdr[a][0]=0;
            fdr[a][1]=0;
        }

    /* Get a value to start comparison*/
    init_flag=0;
    i=0;
    while(init_flag==0)
    {
        if(fdr[i][0]!=0)
        {
            tmpmax=fdr[i][1];
            tmpmaxid=i;
            init_flag=1;
        }
        else
        {
            i++;
        }
    }


    for(i=0; i<N; i++)
        if (fdr[i][0]!=0)
            if(fdr[i][1]>tmpmax)
            {
                tmpmax=fdr[i][1];
                tmpmaxid=i;
            }

   /* Test
    printf(" NO:%d\t",a);

    for(i=0; i<N; i++)
        printf("[%f,%f]",fdr[i][0],fdr[i][1]);

    printf("\n");
    printf("maxid:%d\n",tmpmaxid);

    */


    return tmpmaxid;
}


/*
 Basic Test Functions
*/
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
    double b = 0, c = 1;
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
