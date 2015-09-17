#include"stdio.h"
#include"math.h"
#include"stdlib.h"
#include"time.h"
#include"limits.h"
const int N=30;
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
int main() {

    srand((unsigned)time(NULL));
    double r[2][N];
    int k;
    for(k=0;k<N;k++){
        r[0][k]=(double)rand();
        r[1][k]=k;
        printf("<%f,%d>",r[0][k],(int)r[1][k]);
    }
    printf("\n");
    BubbleSort(r,N);
    for(k=0;k<N;k++)
        printf("<%f,%d>",r[0][k],(int)r[1][k]);
    return 0;
}
