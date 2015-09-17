/*
 Parameters for PSO
*/
const double w = 0.729;
const double c1 = 1.494;
const double c2 = 1.494;

const int dim = 10;
const int N = 50;
const int testnum = 50;
const int gennum = 5000;
const int coornum = (dim + 1) * N;

double coor[coornum];
FILE* pso;

/*
 Useful Functions
*/
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
