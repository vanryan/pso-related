/*
 Parameters for relative radius
*/
const double relStart = 0.24;
const double relEnd = 1.01;
const double relStep = 0.02;
const int relRange = 39;

/*
 Parameters for PSO
*/
const float x = 0.729;
const float c1 = 1.494;
const float c2 = 1.494;

const int dim = 10; // Number of dimensions
const int N = 50; // Population; This will cause a waste of ram
const int testnum = 50; // Running times for experiments
const int gennum = 5000;  // Number of generations
const int coornum = (dim + 1) * N; // That's a tricky way of storing it

//const int aveDegInterval=100; // The interval of recording average degree

double coor[coornum];

/*
 Functions
*/
double getDistance(int a,int b)  // a & b stands for the index of the bird
{
    int i;
    double dis,summup=0;
    for(i=0; i<dim; i++)
    {
        summup += pow(coor[(dim+1)*a+i]-coor[(dim+1)*b+i],2);
    }
    dis = sqrt(summup);
    //printf("dis:%f\tsum:%f\n",dis,summup);
    return dis;
}

double countFIndex(double record[testnum])
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
