int partion(double a[2][N],int p,int r){
    //rand
    time_t t;
    srand((unsigned) time(&t));
    int e=rand()%(r-p+1)+p;
    double tem,tem2;
    tem=a[0][e];
    tem2=a[1][e];
    a[0][e]=a[0][r];
    a[1][e]=a[1][r];
    a[0][r]=tem;
    a[1][r]=tem2;
    double x=a[0][r];
    int i=p-1;
    for (int j=p;j<r;j++){
        if (a[0][j]<=x){
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

void QuickSort(double a[2][N],int p,int r){
    if (p<r){
        int q=partion(a,p,r);
        QuickSort(a,p,q-1);
        QuickSort(a,q+1,r);
    }
}
