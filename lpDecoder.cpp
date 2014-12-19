#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
//#include <stdlib.h>
#include <cstdlib>
#include <cstring>
#include <glpk.h>
using namespace std;

const int SIGMA = 1;
int M,N;
int MaxColDegree,MaxRowDegree;
int *initCode;
double *recvCode;
int **H;
int **A_V,**B_C;
//double var;
int **coef;
int *pColDegree;
int *pRowDegree;
int wCount = 0,rCount = 0;
int *temp_B;

void print_H();
void read_file(){
	int i,j,tmp,*temp_p;
	FILE *fp;
	fp = fopen("alist_arr.txt","r");
	if(!fp){
		printf("open file alisr_arr.txt err");
		return ;
	}
	fscanf(fp,"%d",&N);
	fscanf(fp,"%d",&M);
	fscanf(fp,"%d",&MaxColDegree);
	fscanf(fp,"%d",&MaxRowDegree);
	H = (int **)malloc(M*sizeof(int*));
	for(i=0;i<M;i++){
		H[i] = (int *)malloc(N*sizeof(int));
	}
	for(i=0;i<M;i++)
		for(j=0;j<N;j++)
			H[i][j] = 0;
	pColDegree = (int *)malloc(N*sizeof(int));
	for (temp_p=pColDegree; temp_p<pColDegree+N; temp_p++)
	{
		fscanf(fp, "%d", temp_p);	//读每一列中1的个数到矩阵pColDegree[]中
	}

	pRowDegree = (int *)malloc(M*sizeof(int));
	for (temp_p=pRowDegree; temp_p<pRowDegree+M; temp_p++)
	{
		fscanf(fp,"%d",temp_p);//读每一行中1的个数到矩阵pRowDegree[]中
	}
	A_V = (int **)malloc(N*sizeof(int *));
	for(i = 0;i < N;i ++)
		A_V[i] = (int *)malloc(MaxColDegree*sizeof(int));
	for(i = 0;i < N;i++){
		for(j = 0;j < MaxColDegree;j ++){
			A_V[i][j] = 0;
		}
	}
	
	for(i = 0;i < N;i ++){
		for(j = 0;j < MaxColDegree;j ++){
			fscanf(fp,"%d",&tmp);
			A_V[i][j] = tmp-1;
			if(tmp > 0)
				H[tmp-1][i] = 1;
		}
	}
	B_C = (int **)malloc(M*sizeof(int*));
	for(i = 0;i < M;i ++)
		B_C[i] = (int *)malloc(MaxRowDegree*sizeof(int));
	for(i = 0;i < M;i ++){
		for(j = 0;j < MaxRowDegree;j++){
			fscanf(fp,"%d",&tmp);
			B_C[i][j] = tmp-1;
		//	H[i][tmp-1] = 1;
		}
	}
	fclose(fp);
	print_H();
}

void print_H(){
	int i,j;
	FILE *fp;
	fp = fopen("H","w");
	if(! fp){
		printf("error when open H");
		return ;
	}
	for(i = 0;i < M;i++){
		for(j = 0;j < N;j++){
			fprintf(fp,"%d",H[i][j]);
			putc(' ',fp);
		}
		putc('\n',fp);
	}
	fclose(fp);
}

double Gauss(){
	double ret;
	double UA,UB;
	static double U1,U2;
	double s;
	static double fac;
	static int phase = 0;
	if(phase == 0){
		do{
			UA = (float)rand()/RAND_MAX;
			UB = (float)rand()/RAND_MAX;
			U1 = 1-2*UA;
			U2 = 1-2*UB;
			s = U1*U1+U2*U2;
		}while(s>=1 || s<=0);
		fac = sqrt(-2.0*SIGMA*SIGMA*log(s)/s);
		ret = U1*fac;
	}else{
		ret = U2*fac;
	}		
	phase = 1-phase;
	return ret;
}

void functionObjective(const double &stdev){
	int i;
	initCode = (int *)malloc(N*sizeof(int));
	recvCode = (double *)malloc(N*sizeof(double));
	for(i = 0;i < N;i ++)
		initCode[i] = 0;
	for(i = 0;i < N;i ++)
		recvCode[i] = pow(-1,initCode[i]);//Because all-0
	for(i = 0;i < N;i ++){
		recvCode[i] += Gauss()*stdev;
		recvCode[i] = 2*recvCode[i]/(stdev*stdev);
	}
    FILE *fp;
    fp = fopen("objCoef","w");

    fprintf(fp,"[");
    for (i=0; i<N-1; i++)
    {
        fprintf(fp,"%f , ",recvCode[i]);
    }
    fprintf(fp,"%f ]\n",recvCode[N-1]);

    fclose(fp);
	free(initCode);
//	free(recvCode);
}

long factorial(int n){
	long result = 1;
	while(n > 0){
		result *= n;
		n--;
	}
	return result;
}

long combination(const int n,int i){
	if(n<i){
		printf("error when call combination method");
		return 0;
	}
	if(i==0 || i==n)
		return 1;
	else if(i == 1)
		return n;
	else
	{
		if (i > (n/2))   // C(n,i) = C(n,(n-i))
		{
			i = n-i;
		}

		//最简单的算法，返回 n!/i!*(n-i)! ，但非常容易溢出
		//return factorial(n)/(factorial(i)*factorial(n-i));

		//第二种算法，返回 n*(n-1)*(n-2)*...*(n-i+1) / i! ，并逐级约分
		int *numerator; //存储分子的数组
		int *denominator; //存储分母的数组
		numerator = (int *)malloc(i*sizeof(int));
		denominator = (int *)malloc((i-1)*sizeof(int));

		int k,s;
		for (k=0; k<i; k++)
		{
			numerator[k] = n-k; 
		}//numerator[]存储n，n-1...n-i+1
		for (k=0; k<i-1; k++)
		{
			denominator[k] = i-k;
		}//denominator[]存储i，i-1...2 ,s[0]=i...s[i-2]=2
		for (k=0; k<i; k++)
		{
			for (s=(i-1)-1; s>=0; s--)
			{
				if (denominator[s]==1)
				{
					continue;
				}
				if (numerator[k]%denominator[s] == 0)
				{
					numerator[k] /= denominator[s];
					denominator[s] = 1;
				}
			}
		}//从分子的最低位，和分母逐个进行约分，若已约分，则相应分子除以分母，分母为1
		
		long num=1, den=1;
		for (k=0; k<i; k++)
		{
			num *= numerator[k];
		}
		for (k=0; k<i-1; k++)
		{
			den *= denominator[k];
		}
		
		free(numerator);
		free(denominator);

		return num/den;
	}
}

int combOddNum(const int n){
	int i,sum = 0;
	for(i = 0;i < n;i ++){
		if(i%2 != 0)
			sum += combination(n,i);
	}
	return sum;
}

void matrix(){
	int **A;
	int *B;
	int i,j,k,offset=0,combNum,assign;
	int max_sign_i;
	printf("MaxRowDe=%d\n",MaxRowDegree);
	printf("MaxColDe=%d\n",MaxColDegree);
 	wCount = combOddNum(MaxRowDegree);
	printf("wCount=%d\n",wCount);
	A = (int **)malloc((wCount)*sizeof(int*));
	for(i = 0;i < wCount;i ++)
		A[i] = (int *)malloc(MaxRowDegree*sizeof(int));
	for(i = 0;i < wCount;i++){
		for(j = 0;j < MaxRowDegree; j ++)
			A[i][j] = -1;
	}
	B = (int*) malloc (wCount*sizeof(int));
	temp_B = (int*)malloc(M*wCount*sizeof(int));
	k = 1;
  	while( k <= MaxRowDegree )
  	{
        int *c = (int *)malloc( k*sizeof(int) );
        for (i=0; i<k; i++)
        {
            c[i] = i;
        }

        combNum = combination( MaxRowDegree , k );

        for ( int r = 0; r < combNum; r++)
        {
           B[ r + offset] = k - 1;
        }

        for (assign = 0; assign <= combNum -1; assign++)
        {
            for (i = 0; i < k; i++)
            {
                A[assign + offset][c[i]] =  1;
            }

            max_sign_i = -1;
            for (i=0; i<k; i++)
            {
                if( c[i] < MaxRowDegree-k+i )
                    max_sign_i = i;
            }
            if (max_sign_i == -1)
            {
                if (assign != combNum-1)
                {
                    printf("error in function(enum_comb)");
                    system("pause");exit(0);
                }
            }
            if (max_sign_i != -1)
            {
                c[max_sign_i] = c[max_sign_i]+1;
                for (j=max_sign_i+1; j<k; j++)
                {
                    c[j] = c[j-1]+1;
                }
            }
        }
        offset += combNum;
        k += 2;
        free(c);
	}
	coef = (int **)malloc(rCount*N*sizeof(int*));
	for(i = 0;i < rCount;i ++){
		coef[i] = (int*)malloc(N*sizeof(int));
		for(k = 0;k < N;k++)
			coef[i][k] = 0;
	}
	for(i = 0;i < M; i++){
		for(j = 0;j < wCount;j++){
			for(k = 0;k < MaxRowDegree; k++){
				coef[i*wCount+j][B_C[i][k]] = A[j][k];
			}
		}
	}
	FILE* fp;
	fp = fopen("LPInput","w");
	fprintf(fp,"[\n");
    for ( i = 0; i < rCount; i ++)
    {
        fprintf(fp,"[");
        for ( j = 0; j < N; j++)
        {
            fprintf(fp, j==N-1 ? "%2d " : "%2d , ",coef[i][j]);
        }
        fprintf(fp,i==rCount-1 ? "]\n" : "] ,\n");
    }
    fprintf(fp,"]\n");

    fprintf(fp,"[");
	int tmp = 0;
    for ( i = 0; i < M ; i++)
   	{
       	for ( j = 0; j < wCount ; j++)
       	{
			fprintf(fp,j == wCount-1?"%d ":"%2d , ",B[j]);
			temp_B[tmp++] = B[j];
       	}
   	}
    fprintf(fp,"]\n");
    fclose(fp);
	
	
//	for(i = 0;i < rCount;i ++){
//		printf("%d,",temp_B[i]);
//		//std::cout << "value: " << temp_B[i] << std::endl;
//	}
/*
	for ( i = 0; i < rCount; i++)
    {
      free(coef[i]);
    }
    free(coef);
*/
}

int main(){
	double stdev;
	read_file();
	srand((unsigned)time(NULL));
	matrix();
//	printf("%d,%d\n",M,N);
	FILE *fp;
	fp = fopen("result","w");
	if(! fp){
		printf("err when open file result");
		exit(0);
	}
	for(double snr = 0.0;snr <= 0.0;snr += 0.1){
		int errRate = 0;
		long sumCount = 0;
		int fifty = 50;
		int BEC = 0;
		int WEC = 0;
		stdev=sqrt(1.0/((2.0*(N - M)/N)*pow(10,snr/10)));
		int *ia = (int*)malloc((1+N*rCount)*sizeof(int));
		int *ja = (int*)malloc((N*rCount+1)*sizeof(int));
		double *ar = (double*)malloc((1+N*rCount)*sizeof(double));
		double *x = (double *)malloc(sizeof(double)*(N+1));
		while(true){
			sumCount++;
			functionObjective(stdev);
			glp_prob *lp;
			int i,j,k=0;			
			lp = glp_create_prob();
			glp_set_prob_name(lp,"lpdecoder");
			glp_set_obj_dir(lp,GLP_MIN);
			glp_add_rows(lp,rCount);
			for(i = 0;i < rCount;i ++){
				std::stringstream strOs;
				strOs << i+1;
			//	glp_set_row_name(lp,i+1,strOs.str().c_str());//""+i
				glp_set_row_bnds(lp,i+1,GLP_UP,0.0,temp_B[i]);
			}
			glp_add_cols(lp,N);
			for(i = 1;i <= N;i ++){
//				std::stringstream strOs;
//				strOs << i;
//				string s;
//				sprintf(s,"%d",i);
	//			glp_set_col_name(lp,i,"x");//""+x[i]   
//				glp_set_col_bnds(lp,i,GLP_LO,0.0,0.0);
				glp_set_obj_coef(lp,i,recvCode[i-1]);
			}
			int tmp=1;
			for(i = 1;i <= rCount;i ++){
				for(j = 1;j <= N;j ++){
					ia[tmp] = i,ja[tmp] = j,ar[tmp++] = coef[i-1][j-1];
//				std::cout << coef[i-1][j-1] << " ";	
//				std::cout << endl;
				}
			}
			i--;j--;tmp--;
			//std::cout << i << j << tmp << std::endl;
			glp_load_matrix(lp,tmp,ia,ja,ar);
			glp_simplex(lp,NULL);
			bool wErr = false;
			for(i = 0;i < N;i ++){
				x[i+1] = glp_get_col_prim(lp,i+1);
//				std::cout << x[i+1] << " ";
				if(x[i+1] != 0){
					wErr = true;
					BEC += 1;
				}
			}
			glp_delete_prob(lp);
			if(wErr)
				WEC += 1;
			if(WEC >= 1){
				printf("snr = %f,errCodeRate = %f,errFrameRate = %f\n",snr,(float)BEC/(N*sumCount),(float)WEC/sumCount);
				fprintf(fp,"snr = %f,errCodeRate = %f,errFrameRate = %f\n",snr,(float)BEC/(N*sumCount),(float)WEC/sumCount);
				break;
			}
		}
		delete[] x;	
		delete[] ia;
		delete[] ja;
		delete[] ar;
	//	double frameErrRate = (double)errRate/sumCount;
	//	fprintf(fp,"when snr = %lf,frame error rate = %lf\n",snr,frameErrRate);
	//	printf("when snr=%lf,frame erroe rate=%lf\n",snr,frameErrRate);
	}
	fclose(fp);
	return 0;
}
