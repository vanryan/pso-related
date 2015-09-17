
const int basic_func_num=10;
const int search_range=5;  // The search range of the composition function
const int f_bias=0;
const int bias[basic_func_num]={0,100,200,300,400,500,600,700,800,900};
const int pre_c=2000;

const double hyb_func_goalnedgenvmax[6][3]={{100,search_range,search_range},{100,search_range,search_range}\
            ,{100,search_range,search_range},{400,search_range,search_range},{100,search_range,search_range}\
            ,{100,search_range,search_range}}; // Row:function index; 1st rank--goal; 2nd rank--edge; 3rd rank--vmax

//srand((unsigned)time(NULL));


void update_opt_arr(double arr[basic_func_num][dim])
{
    int i,j;
    for(i=0;i<basic_func_num-1;i++)
        for(j=0;j<dim;j++)
            arr[i][j] =(double)rand() / RAND_MAX *  2 * search_range - search_range;
    for(j=0;j<dim;j++)
        arr[basic_func_num-1][j]=0;
}


template <class T>

int getArrayLen(T& array)
{
return (sizeof(array) / sizeof(array[0]));
}

double get_max_from_arr(double arr[],int arr_len){
    int i;
    double temp=arr[0];
    for(i=1;i<arr_len;i++){
        if(arr[i]>temp)
            temp=arr[i];
    }
    return temp;
}

void handle_weight(double arr[basic_func_num]){
    double max_wei,sum_wei=0;
    int i,j;
    j=getArrayLen(arr);
    max_wei=get_max_from_arr(arr,j);
    for(i=0;i<basic_func_num;i++){
        if(max_wei!=arr[i])
            arr[i]*=(1-pow(max_wei,10));
        sum_wei+=arr[i];
    }
    for(i=0;i<basic_func_num;i++)
        arr[i]/=sum_wei;
}
/*
void hyb_func1_init(double fmax[basic_func_num], double or_matrix[dim][dim], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;
    double zedd[dim],lamda[basic_func_num];
    double (*funcp[basic_func_num])(double coord[], int);

    for(i=0;i<basic_func_num;i++){
        funcp[i]=Func1;
        lamda[i]=0.05;
        for(j=0;j<dim;j++)
            opt_arr_old[i][j]=0;
    }

    //initializing
    update_opt_arr(opt_arr);
    //End of initializing

    ort_rot_mat_create(or_matrix);

    for(i=0;i<basic_func_num;i++){
        for(j=0;j<dim;j++)
            zedd[j]=search_range/lamda[i];
        multp_mat2(zedd,or_matrix);
        fmax[i]=funcp[i](zedd,dim);
    }
}
*/
void hyb_func1_init(double fmax[basic_func_num], double or_matrix[dim][dim], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;

    for(i=0;i<basic_func_num;i++){
        for(j=0;j<dim;j++)
            opt_arr_old[i][j]=0;
    }

    //initializing opt_arr
	ifstream infile1("cf1-opt.txt",ios::in);
	if(!infile1) {
            cerr<<"open file1 error!"<<endl;
            exit(1);
	}
	for(i=0;i<basic_func_num;i++)
		for(j=0;j<dim;j++)
			infile1>>opt_arr[i][j];
    infile1.close();
    //End of initializing opt_arr

    //Print them out
    printf("Opt_Arr:\n");
    for(i=0;i<basic_func_num;i++){
		for(j=0;j<dim;j++)
            printf("%f,",opt_arr[i][j]);
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");



	//initializing or_matrix
	FILE* read;
	read = fopen("cf1-ortmat.txt", "rt");
	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			fscanf(read,"%lf,",&or_matrix[i][j]);
	fclose(read);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++)
			printf("%f,",or_matrix[i][j]);
        printf("\n");
	}
	printf("\n");
    //End of initializing or_matrix

    //inputting fmax
    ifstream infilem("cf1-fmax.txt",ios::in);
	if(!infilem) {
            cerr<<"open file error!"<<endl;
            exit(1);
	}
	for(i=0;i<basic_func_num;i++)
        infilem>>fmax[i];
    infilem.close();

    for(i=0;i<basic_func_num;i++){
        printf("%f,",fmax[i]);
    }
    printf("\n");

}

void hyb_func1(double coord[], int a, double or_matrix[dim][dim], double fmax[basic_func_num], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){//Composition Function 1
    int i,j;
    int theta[basic_func_num];
    double lamda[basic_func_num],weight[basic_func_num],temp1;
    double fit[basic_func_num],fval[basic_func_num];
    double mult_mat[dim];

    double (*funcp[basic_func_num])(double coord[], int);
    for(i=0;i<basic_func_num;i++){
        funcp[i]=Func1;
        theta[i]=1;
        lamda[i]=0.05;
    }

    for(i=0;i<basic_func_num;i++){
        temp1=0;

        for(j=0;j<dim;j++)
            temp1+=pow(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j],2);
        weight[i]=temp1;
        // -------------------------------
        for(j=0;j<dim;j++)
            mult_mat[j]=(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j])/lamda[i];

        multp_mat2(mult_mat,or_matrix);

        fit[i]= funcp[i](mult_mat,dim);
        fval[i]=pre_c*fit[i]/fmax[i];
    }
     // -------------------------------

    for(i=0;i<basic_func_num;i++)
        weight[i]=exp(- weight[i]/( 2*dim*pow(theta[i],2) ));
    handle_weight(weight);

    coord[(dim+1)*a+dim]=f_bias;
    for(i=0;i<basic_func_num;i++)
        coord[(dim+1)*a+dim]+=weight[i]*(fval[i]+bias[i]);
}

void hyb_func2(double coord[], int a, double or_matrix[dim][dim], double fmax[basic_func_num], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;
    int theta[basic_func_num];
    double lamda[basic_func_num],weight[basic_func_num],temp1;
    double fit[basic_func_num],fval[basic_func_num];
    double mult_mat[dim];


    double (*funcp[basic_func_num])(double coord[], int);
    for(i=0;i<basic_func_num;i++){
        funcp[i]=Func4;
        theta[i]=1;
        lamda[i]=0.05;
        for(j=0;j<dim;j++)
            opt_arr_old[i][j]=0;
    }

     for(i=0;i<basic_func_num;i++){
        temp1=0;

        for(j=0;j<dim;j++)
            temp1+=pow(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j],2);
        weight[i]=temp1;
        // -------------------------------
        for(j=0;j<dim;j++)
            mult_mat[j]=(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j])/lamda[i];

        multp_mat2(mult_mat,or_matrix);

        fit[i]= funcp[i](mult_mat,dim);
        fval[i]=pre_c*fit[i]/fmax[i];
    }
     // -------------------------------

    for(i=0;i<basic_func_num;i++)
        weight[i]=exp(- weight[i]/( 2*dim*pow(theta[i],2) ));
    handle_weight(weight);

    coord[(dim+1)*a+dim]=f_bias;
    for(i=0;i<basic_func_num;i++)
        coord[(dim+1)*a+dim]+=weight[i]*(fval[i]+bias[i]);
}
void hyb_func3(double coord[], int a, double or_matrix[dim][dim], double fmax[basic_func_num], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;
    int theta[basic_func_num];
    double lamda[basic_func_num],weight[basic_func_num],temp1;
    double fit[basic_func_num],fval[basic_func_num];
    double mult_mat[dim];


    double (*funcp[basic_func_num])(double coord[], int);
    for(i=0;i<basic_func_num;i++){
        funcp[i]=Func4;
        theta[i]=1;
        lamda[i]=1;
        for(j=0;j<dim;j++)
            opt_arr_old[i][j]=0;
    }

     for(i=0;i<basic_func_num;i++){
        temp1=0;

        for(j=0;j<dim;j++)
            temp1+=pow(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j],2);
        weight[i]=temp1;
        // -------------------------------
        for(j=0;j<dim;j++)
            mult_mat[j]=(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j])/lamda[i];

        multp_mat2(mult_mat,or_matrix);

        fit[i]= funcp[i](mult_mat,dim);
        fval[i]=pre_c*fit[i]/fmax[i];
    }
     // -------------------------------

    for(i=0;i<basic_func_num;i++)
        weight[i]=exp(- weight[i]/( 2*dim*pow(theta[i],2) ));
    handle_weight(weight);

    coord[(dim+1)*a+dim]=f_bias;
    for(i=0;i<basic_func_num;i++)
        coord[(dim+1)*a+dim]+=weight[i]*(fval[i]+bias[i]);
}

void hyb_func4_init(double fmax[basic_func_num], double or_matrix[dim][dim], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;

    for(i=0;i<basic_func_num;i++){
        for(j=0;j<dim;j++)
            opt_arr_old[i][j]=0;
    }

    //initializing opt_arr
	ifstream infile4("cf4-opt.txt",ios::in);
	if(!infile4) {
            cerr<<"open file4 error!"<<endl;
            exit(1);
	}
	for(i=0;i<basic_func_num;i++)
		for(j=0;j<dim;j++)
			infile4>>opt_arr[i][j];
    infile4.close();
    //End of initializing opt_arr

    //Print them out
    printf("Opt_Arr:\n");
    for(i=0;i<basic_func_num;i++){
		for(j=0;j<dim;j++)
            printf("%f,",opt_arr[i][j]);
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");



	//initializing or_matrix
	FILE* read;
	read = fopen("cf4-ortmat.txt", "rt");
	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			fscanf(read,"%lf,",&or_matrix[i][j]);
	fclose(read);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++)
			printf("%f,",or_matrix[i][j]);
        printf("\n");
	}
	printf("\n");
    //End of initializing or_matrix

    //inputting fmax
    ifstream infilem("cf4-fmax.txt",ios::in);
	if(!infilem) {
            cerr<<"open file error!"<<endl;
            exit(1);
	}
	for(i=0;i<basic_func_num;i++)
        infilem>>fmax[i];
    infilem.close();

    for(i=0;i<basic_func_num;i++){
        printf("%f,",fmax[i]);
    }
    printf("\n");
}

void hyb_func4(double coord[], int a, double or_matrix[dim][dim], double fmax[basic_func_num], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;
    int theta[basic_func_num];
    double lamda[basic_func_num],weight[basic_func_num],temp1;
    double fit[basic_func_num],fval[basic_func_num];
    double mult_mat[dim];


    double (*funcp[basic_func_num])(double coord[], int);


    for(i=0;i<basic_func_num;i++){
        theta[i]=1;
        for(j=0;j<dim;j++)
            opt_arr_old[i][j]=0;
    }

    for(i=0;i<2;i++){
        funcp[i]=Func5;
        lamda[i]=0.15625;
    }
    for(i=2;i<4;i++){
        funcp[i]=Func3;
        lamda[i]=1;
    }
    for(i=4;i<6;i++){
        funcp[i]=Func2;
        lamda[i]=10;
    }
    for(i=6;i<8;i++){
        funcp[i]=Func4;
        lamda[i]=0.05;
    }
    for(i=8;i<10;i++){
        funcp[i]=Func1;
        lamda[i]=0.05;
    }

     for(i=0;i<basic_func_num;i++){
        temp1=0;

        for(j=0;j<dim;j++)
            temp1+=pow(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j],2);
        weight[i]=temp1;
        // -------------------------------
        for(j=0;j<dim;j++)
            mult_mat[j]=(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j])/lamda[i];

        multp_mat2(mult_mat,or_matrix);

        fit[i]= funcp[i](mult_mat,dim);
        fval[i]=pre_c*fit[i]/fmax[i];
    }
     // -------------------------------

    for(i=0;i<basic_func_num;i++)
        weight[i]=exp(- weight[i]/( 2*dim*pow(theta[i],2) ));
    handle_weight(weight);

    coord[(dim+1)*a+dim]=f_bias;
    for(i=0;i<basic_func_num;i++)
        coord[(dim+1)*a+dim]+=weight[i]*(fval[i]+bias[i]);
}

void hyb_func5_init(double fmax[basic_func_num], double or_matrix[dim][dim], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;

    for(i=0;i<basic_func_num;i++){
        for(j=0;j<dim;j++)
            opt_arr_old[i][j]=0;
    }

    //initializing opt_arr
	ifstream infile5("cf5-opt.txt",ios::in);
	if(!infile5) {
            cerr<<"open file5 error!"<<endl;
            exit(1);
	}
	for(i=0;i<basic_func_num;i++)
		for(j=0;j<dim;j++)
			infile5>>opt_arr[i][j];
    infile5.close();
    //End of initializing opt_arr

    //Print them out
    printf("Opt_Arr:\n");
    for(i=0;i<basic_func_num;i++){
		for(j=0;j<dim;j++)
            printf("%f,",opt_arr[i][j]);
        printf("\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");



	//initializing or_matrix
	FILE* read;
	read = fopen("cf5-ortmat.txt", "rt");
	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			fscanf(read,"%lf,",&or_matrix[i][j]);
	fclose(read);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++)
			printf("%f,",or_matrix[i][j]);
        printf("\n");
	}
	printf("\n");
    //End of initializing or_matrix

//inputting fmax
    ifstream infilem("cf5-fmax.txt",ios::in);
	if(!infilem) {
            cerr<<"open file error!"<<endl;
            exit(1);
	}
	for(i=0;i<basic_func_num;i++)
        infilem>>fmax[i];
    infilem.close();

    for(i=0;i<basic_func_num;i++){
        printf("%f,",fmax[i]);
    }
    printf("\n");
}

void hyb_func5(double coord[], int a, double or_matrix[dim][dim], double fmax[basic_func_num], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;
    int theta[basic_func_num];
    double lamda[basic_func_num],weight[basic_func_num],temp1;
    double fit[basic_func_num],fval[basic_func_num];
    double mult_mat[dim];

    double (*funcp[basic_func_num])(double coord[], int);

    for(i=0;i<basic_func_num;i++)
        theta[i]=1;

    for(i=0;i<2;i++){
        lamda[i]=0.2;
        funcp[i]=Func3;
    }
    for(i=2;i<4;i++){
        lamda[i]=10;
        funcp[i]=Func2;
    }
    for(i=4;i<6;i++){
        lamda[i]=0.05;
        funcp[i]=Func4;
    }
    for(i=6;i<8;i++){
        lamda[i]=0.15625;
        funcp[i]=Func5;
    }
    for(i=8;i<10;i++){
        lamda[i]=0.05;
        funcp[i]=Func1;
    }

    for(i=0;i<basic_func_num;i++){
        temp1=0;

        for(j=0;j<dim;j++)
            temp1+=pow(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j],2);
        weight[i]=temp1;
        // -------------------------------
        for(j=0;j<dim;j++)
            mult_mat[j]=(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j])/lamda[i];

        multp_mat2(mult_mat,or_matrix);

        fit[i]= funcp[i](mult_mat,dim);
        fval[i]=pre_c*fit[i]/fmax[i];
    }
     // -------------------------------

    for(i=0;i<basic_func_num;i++)
        weight[i]=exp(- weight[i]/( 2*dim*pow(theta[i],2) ));
    handle_weight(weight);

    coord[(dim+1)*a+dim]=f_bias;
    for(i=0;i<basic_func_num;i++)
        coord[(dim+1)*a+dim]+=weight[i]*(fval[i]+bias[i]);
}

void hyb_func6_init(double fmax[basic_func_num], double or_matrix[dim][dim], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;
    double zedd[dim],lamda[basic_func_num];
    double (*funcp[basic_func_num])(double coord[], int);

    for(i=0;i<basic_func_num;i++){
        for(j=0;j<dim;j++)
            opt_arr_old[i][j]=0;
    }

    lamda[0]=0.02;
    lamda[1]=0.04;
    lamda[2]=3;
    lamda[3]=4;
    lamda[4]=0.025;
    lamda[5]=0.03;
    lamda[6]=0.109375;
    lamda[7]=0.125;
    lamda[8]=0.045;
    lamda[9]=0.05;

    for(i=0;i<2;i++){
        funcp[i]=Func3;
    }
    for(i=2;i<4;i++){
        funcp[i]=Func2;
    }
    for(i=4;i<6;i++){
        funcp[i]=Func4;
    }
    for(i=6;i<8;i++){
        funcp[i]=Func5;
    }
    for(i=8;i<10;i++){
        funcp[i]=Func1;
    }

    //initializing
    update_opt_arr(opt_arr);
    //End of initializing

    ort_rot_mat_create(or_matrix);

    for(i=0;i<basic_func_num;i++){
        for(j=0;j<dim;j++)
            zedd[j]=search_range/lamda[i];
        multp_mat2(zedd,or_matrix);
        fmax[i]=funcp[i](zedd,dim);
    }
}

void hyb_func6(double coord[], int a, double or_matrix[dim][dim], double fmax[basic_func_num], double opt_arr[basic_func_num][dim], double opt_arr_old[basic_func_num][dim]){
    int i,j;
    int theta[basic_func_num];
    double lamda[basic_func_num],weight[basic_func_num],temp1;
    double fit[basic_func_num],fval[basic_func_num];
    double mult_mat[dim];


    double (*funcp[basic_func_num])(double coord[], int);

    double acccctor=0.1;
    for(i=0;i<basic_func_num;i++){
        theta[i]=acccctor;
        acccctor+=0.1;
        for(j=0;j<dim;j++)
            opt_arr_old[i][j]=0;
    }

    lamda[0]=0.02;
    lamda[1]=0.04;
    lamda[2]=3;
    lamda[3]=4;
    lamda[4]=0.025;
    lamda[5]=0.03;
    lamda[6]=0.109375;
    lamda[7]=0.125;
    lamda[8]=0.045;
    lamda[9]=0.05;

    for(i=0;i<2;i++){
        funcp[i]=Func3;
    }
    for(i=2;i<4;i++){
        funcp[i]=Func2;
    }
    for(i=4;i<6;i++){
        funcp[i]=Func4;
    }
    for(i=6;i<8;i++){
        funcp[i]=Func5;
    }
    for(i=8;i<10;i++){
        funcp[i]=Func1;
    }

     for(i=0;i<basic_func_num;i++){
        temp1=0;

        for(j=0;j<dim;j++)
            temp1+=pow(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j],2);
        weight[i]=temp1;
        // -------------------------------
        for(j=0;j<dim;j++)
            mult_mat[j]=(coord[(dim+1)*a+j]-opt_arr[i][j]+opt_arr_old[i][j])/lamda[i];

        multp_mat2(mult_mat,or_matrix);

        fit[i]= funcp[i](mult_mat,dim);
        fval[i]=pre_c*fit[i]/fmax[i];
    }
     // -------------------------------

    for(i=0;i<basic_func_num;i++)
        weight[i]=exp(- weight[i]/( 2*dim*pow(theta[i],2) ));
    handle_weight(weight);

    coord[(dim+1)*a+dim]=f_bias;
    for(i=0;i<basic_func_num;i++)
        coord[(dim+1)*a+dim]+=weight[i]*(fval[i]+bias[i]);
}
