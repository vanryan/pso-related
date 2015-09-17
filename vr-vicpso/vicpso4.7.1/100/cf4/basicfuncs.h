
//Basic Functions

double Func1(double coord[], int dimens)                              //Sphere
{
    double fit=0;
    int i;
    for(i = 0; i < dimens; i++)
        fit += pow(coord[i], 2);
    return fit;
}

double Func2(double coord[], int dimens)  //Weierstrass
{
    double fit=0;
    int i,j,kmax=20;
    double para_a=0.5,para_b=3,temp1=0,temp2=0;
    for(j=0; j<=kmax; j++)
        temp2+=pow(para_a,j)*cos(PI*pow(para_b,j));
    temp2*=dimens;
    for(i=0; i<dimens; i++)
    {
        for(j=0; j<=kmax; j++)
        {
            temp1+=pow(para_a,j)*cos(2*PI*pow(para_b,j)*(coord[i]+0.5));
        }
    }
    fit += temp1;
    fit -= temp2;
    return fit;
}

double Func3(double coord[], int dimens)                             //Rastrigrin
{
    double fit=0;
    int i;
    for(i = 0; i < dimens; i++)
        fit += pow(coord[i], 2) - 10 * cos(2 * 3.14 * coord[i]) + 10;
    return fit;
}

double Func4(double coord[], int dimens)                           //Griewank
{
    double fit=0;
    int i;
    float b = 0, c = 1;
    for(i = 0; i < dimens; i++)
    {
        b +=  pow(coord[i]  , 2) / 4000;
        c  *=  cos(coord[i]  / sqrt(i + 1));
    }
    fit = 1 + b - c;
    return fit;
}

double Func5(double coord[], int dimens) //Ackley
{
    double fit=0;
    int i;
    double fa=0,fb=0;
    for(i = 0; i < dimens; i++)
    {
        fa+=coord[i]*coord[i];
        fb+=cos(2*PI*coord[i]);
    }
    fit+=-20*exp( -0.2*sqrt(fa/dimens) ) - exp(fb/dimens) + 20 +exp(1);
    return fit;
}

double Func6(double coord[], int dimens)                             //Rosenbrock
{
    double fit=0;
    int i;
    for(i = 0; i < dimens; i++)
    {
        if(i < dimens - 1)
            fit += 100 * pow(( coord[i+1]  -  pow(coord[i] , 2) ) , 2)  + pow(coord[i] - 1 , 2);
        else
            fit += 0;//100*pow(( coord[(dimens+1)*a]  -  pow(coord[(dimens+1)*a+i] ,2) ) ,2)  + pow(coord[(dimens+1)*a+i] - 1 ,2);
    }
    return fit;
}

