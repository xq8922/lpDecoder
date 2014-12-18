#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>
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
	H = (int **)malloc(M*sizeof(int));
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
	for(i = 0;i < N;i++)
		for(j = 0;j < MaxColDegree;j ++)
			A_V[i][j] = 0;
	for(i = 0;i < N;i ++){
		for(j = 0;j < MaxColDegree;j ++){
			fscanf(fp,"%d",&tmp);
			A_V[i][j] = tmp-1;
//			if(tmp > 0)
//				H[tmp-1][i] = 1;
		}
	}
	B_C = (int **)malloc(M*sizeof(int*));
	for(i = 0;i < M;i ++)
		B_C[i] = (int *)malloc(MaxRowDegree*sizeof(int));
	for(i = 0;i < M;i ++)
		for(j = 0;j < MaxRowDegree;j++){
			fscanf(fp,"%d",&tmp);
			B_C[i][j] = tmp-1;
			if(tmp > 0)
				H[i][tmp-1] = 1;
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

void functionObjective(const double var){
	int i;
	initCode = (int *)malloc(N*sizeof(int));
	recvCode = (double *)malloc(N*sizeof(double));
	for(i = 0;i < N;i ++)
		initCode[i] = 0;
	for(i = 0;i < N;i ++)
		recvCode[i] = pow(-1,initCode[i]);//Because all-0
	for(i = 0;i < N;i ++){
		recvCode[i] += Gauss()*var;
		recvCode[i] = 2*recvCode[i]/(var*var);
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
    free(recvCode);
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
	temp_B = (int*)malloc(wCount*sizeof(int));
	k = 1;
/*
	while(k <= MaxRowDegree){
		int *c = (int *)malloc(k*sizeof(int));
		for(i = 0;i < k;i++){
			c[i] = i;
		}
		comNum = combination(MaxRowDegree,k);
		for(int r = 0;r < comNum;r++)
			B[r+offset] = k-1;
		for(assign = 0;assign <= comNum-1;assign++){
			for(i = 0;i < k;i++){
				A[assign+offset][c[i]] = 1;
			}
			maxSignI = -1;			
			for(i = 0;i < k;i ++){
				if(c[i] < MaxRowDegree-k+i)
					maxSignI = i;
			}
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
	}*/

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

	printf("print A\n");
	for(i = 0;i < 10;i ++){
		for(j = 0;j < MaxRowDegree;j ++){
			printf("%d ",A[i][j]);
		}
		printf("\n");
	}
/*	printf("print B\n");
cout<<"[";
    for ( i = 0; i < M ; i++)
    {
        for ( j = 0; j < wCount ; j++)
        {
            cout<<B[j]<<" ";
        }
    }
    cout<<"]\n";
*/
//	for(j = 0;j < wCount; j++)
//		printf("%d ",B[j]);
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
/*
	for ( i = 0; i < rCount; i++)
    {
      free(coef[i]);
    }
    free(coef);
*/
}

int main(){
	double var;
	read_file();
	srand((unsigned)time(NULL));
	matrix();
	printf("%d,%d\n",M,N);
	FILE *fp;
	fp = fopen("result","w");
	if(! fp){
		printf("err when open file result");
		exit(0);
	}
	const float ESP = 0.000001;
	for(double snr = 2.0;snr <= 5.0;snr += 0.1){
		int errRate = 0;
		long sumCount = 0;
		int fifty = 50;
int BEC = 0,WEC = 0;
		var = sqrt(1.0/((2.0*(N-M)/N)*pow(10,snr/10)));
		while(true){
			
			functionObjective(var);
			//decoding . judge if right sumCount++;else errRate++;
			//if(resultCode == 0)sumCOunt++;else errRate++;
			glp_prob *lp;
			int *ia = (int*)malloc(rCount*sizeof(int));
			int *ja = (int*)malloc(N*sizeof(int));
			double *ar = (double*)malloc(N*rCount*sizeof(double));
			int i,j,k=1;
			double *x = (double *)malloc(sizeof(double)*N);
			lp = glp_create_prob();
		glp_set_prob_name(lp,"lpdecoder");
		glp_set_obj_dir(lp,GLP_MIN);
			glp_add_rows(lp,rCount);
			for(i = 1;i <= rCount;i += 2){
		glp_set_row_name(lp,i,"");//""+i
		glp_set_row_bnds(lp,i+1,GLP_UP,0.0,temp_B[i-1]);
			}
		glp_add_cols(lp,3);
			for(i = 1;i < N;i += 3){
		glp_set_col_name(lp,i,"");//""+x[i]
		glp_set_col_bnds(lp,i,GLP_LO,0.0,0.0);
		glp_set_obj_coef(lp,i,recvCode[i-1]);
			}
		int tmp=0;
		for(i = 0;i < rCount;i ++){
			for(j = 0;j < N;j ++){
			ia[i] = i,ja[j] = j,ar[tmp++] = coef[i][j];
			}
		}
		glp_load_matrix(lp,tmp-1,ia,ja,ar);
		glp_simplex(lp,NULL);
		for(i = 0;i < N;i ++){
		x[i] = glp_get_col_prim(lp,i+1);
		printf("\nx = %g",x[i]);
		}
		glp_delete_prob(lp);
			bool wErr = false;
			float *f = new float[N];
			for(int i = 0;i < N; i++){
				f[i] = 0;
			}
			for(int i = 0;i < N;i ++){
				if(f[i] > ESP){
					wErr = true;
					BEC += 1;
				}
			}
			delete[] f;
			if(wErr)
				WEC += 1;
			if(WEC >= 50){
				//fprintf
				break;
			}
break;
		}
	//	double frameErrRate = (double)errRate/sumCount;
	//	fprintf(fp,"when snr = %lf,frame error rate = %lf\n",snr,frameErrRate);
	//	printf("when snr=%lf,frame erroe rate=%lf\n",snr,frameErrRate);
	}

	return 0;
}
