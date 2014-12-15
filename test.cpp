#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>
using namespace std;

const int SIGMA = 1;
int M,N;
int MaxColDegree,MaxRowDegree;
int *initCode,*recvCode;
int **H;
int **A_V,**B_C;
double var;
int **coef;

void print_H();
void read_file(){
	int i,j,tmp;
	FILE *fp;
	fp = fopen("alist_arr.txt","r");
	if(fp != NULL){
		printf("open file alisr_arr.txt err");
		return ;
	}
	fscanf(fp,"%d",&N);
	fscanf(fp,"%d",&M);
	fscanf(fp,"%d",&MaxColDegree);
	fscanf(fp,"%d",&MaxRowDegree);
	H = (int **)malloc(M*sizeof(int));
	for(i=0;i<M;i++){
		H[i] = (int *)malloc(N*sizeof(int));
	}
	for(i=0;i<M;i++)
		for(j=0;j<N;j++)
			H[i][j] = 0;
	A_V = (int **)malloc(N*sizeof(int *));
	for(i = 0;i < N;i ++)
		A_V[i] = (int *)malloc(MaxColDegree*sizeof(int));
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
	for(i = 0;i < M;i ++)
		for(j = 0;j < MaxRowDegree;j++){
			fscanf(fp,"%d",&tmp);
			B_C[i][j] = tmp-1;
		}
	fclose(fp);
	print_H();
}

void print_H(){
	int i,j;
	FILE *fp;
	fp = fopen("H_file.txt","w");
	if(! fp){
		printf("error when open H_file");
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

void functionObjective(){
	int i;
	initCode = (int *)malloc(N*sizeof(int));
	recvCode = (int *)malloc(M*sizeof(int));
	for(i = 0;i < N;i ++)
		initCode[i] = 0;
	for(i = 0;i < N;i ++)
		recvCode[i] = pow(-1,initCode[i]);//Because all-0
	for(i = 0;i < N;i ++){
		recvCode[i] += Gauss()*var;
		recvCode[i] = 2*recvCode[i]/(var*var);
	}
}

long factorial(int n){
	long result = 1;
	while(n > 0){
		result *= n;
		n--;
	}
	return result;
}

long combination(const int n,const int k){
	if(n<k){
		printf("error when call combination method");
		return 0;
	}
	if(k==0 || k==n)
		return 1;
	else if(k == 1)
		return n;
	else
		return factorial(n)/(factorial(k)*factorial(n-k));
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
	int wCount=0,rCount,i,j,k,offset=0,comNum,assign;
	int maxSignI;
 	wCount = combOddNum(MaxRowDegree);
	A = (int **)malloc((wCount)*sizeof(int*));
	for(i = 0;i < wCount;i ++)
		A[i] = (int *)malloc(MaxRowDegree*sizeof(int));
	for(i = 0;i < wCount;i++){
		for(j = 0;j < MaxRowDegree; j ++)
			A[i][j] = -1;
	}
	B = (int*) malloc (wCount*sizeof(int));
	k = 1;
	while(k <= MaxRowDegree){
		int *c = (int *)malloc(k*sizeof(int));
		for(i = 0;i < k;i++){
			c[i] = i;
		}
		comNum = combination(MaxRowDegree,k);
		for(int r = 0;r < comNum;r++)
			B[r+offset] = k-1;
		for(assign = 0;assign <= comNum-1;assign++){
			for(i = 0;i < k;i++)
				A[assign+offset][c[i]] = 1;
			maxSignI = -1;			
			for(i = 0;i < k;i ++)
				if(c[i] < MaxRowDegree-k+i)
					maxSignI = i;
			if(maxSignI == -1){
				if(assign != comNum-1){
					printf("error when maxSignI");
					return ;
				}
			}
			if(maxSignI != -1){
				c[maxSignI] += 1;
				for(j = maxSignI + 1;j < k;j ++)
					c[j] = c[j-1]+1;
			}
		}
		offset += comNum;
		k += 2;
		free(c);
	}
	printf("print A");
	for(i = 0;i < wCount;i ++){
		for(j = 0;j < MaxRowDegree;j ++){
			printf("%d ",A[i][j]);
		printf("\n");
		}
	}
	printf("print B");
	for(j = 0;j < wCount; j++)
		printf("%d ",B[j]);
	printf("\n");
	rCount = M*wCount;
	coef = (int **)malloc(rCount*sizeof(int*));
	for(i = 0;i < rCount;i ++){
		coef[i] = (int*)malloc(N*sizeof(int));
		for(k = 0;k < N;k++)
			coef[i][k] = 0;
	}
	for(i = 0;i < M; i++){
		for(j = 0;j < wCount;j++){
			for(k = 0;k < MaxRowDegree; k++){
				coef[i*wCount+j][B_C[i][j]] = A[j][k];
			}
		}
	}
}

int main(){
	read_file();
	srand((unsigned)time(NULL));
	matrix();
	for(double snr = 2.0;snr < 5.0;snr += 0.1){
		var = sqrt(1.0/((2.0*(N-M)/N)*pow(10,snr/10)));
		while(true){
			functionObjective();
		}
	}

	return 0;
}
