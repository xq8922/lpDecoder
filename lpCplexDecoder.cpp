#include <ilcplex/ilocplex.h>
#include <stdio.h>
#include <time.h>

using namespace std;

ILOSTLBEGIN
//#define SNR   6.0 //2.6

#define ESP 0.000001  //��0�Ƚϵľ���
double p=0.20;
int **H; //У�����
int N; //��
int M; //��
int MaxColDegree;
int MaxRowDegree;
int *pColDegree;
int *pRowDegree;
int **A_V; //��Ϣ�ڵ�
int **B_C; //У��ڵ�

int *codeword;	

double *r_code;	 //��ʾ���Թ滮Ŀ�꺯���еĲ���

struct w_Enum
{
	int CNode;
	int vecSize;
//	vector<int> vec;
	int *vec;
};

struct wTofMatrix
{
   int ** A;
   int *B;
   int rNum;
   int cNum;
};

void read_H();
void print_H();
void read_code();
void LP_target_function_construct();
void LP_matrix_construct();
void enum_comb(int, int, int, int, int, w_Enum* );
int comb_even_number_sum(int n);		//����ÿ��У��ڵ�C(0,n)+C(2,n)+C(4,n)...
int comb_odd_number_sum(int n);		//����ÿ��У��ڵ�C(1,n)+C(3,n)+C(5,n)...
long factorial(int);
long combination(int, int);
wTofMatrix* wTof(int rDeg);
void freeWtoFMatrix(wTofMatrix *w);
void freeABH();

int main()
{
	long dCount = 0; //decode count
	int WEC = 0; // word error count
	int BEC = 0; // bit error count

  	srand((unsigned)time(NULL)); 
	read_H();
	float snrStart = 2.0,snrEnd = 6.0,snrInterval = 1.0;
	while(1)
	{
	   cout<<"��������ȷ�Χ�Ͳ�������ʽ: ��ʼֵ ��ֵֹ ��������"<<endl;
	   cin>>snrStart>>snrEnd>>snrInterval;
	   if(snrStart <= snrEnd){break;}
	}
	p = snrStart;
	ofstream out("result.txt",ofstream::app);
	out<<"************************************************************************"<<endl;
	out<<"LDPC code :"<<"("<<MaxColDegree<<","<<MaxRowDegree<<"), N="<<N<<",M="<<M<<",R="<<(float)(N-M)/N<<endl;
	out.close();

	LP_matrix_construct();
	freeABH();

	//return 0;
	clock_t start,end;
	while(true)
	{
		LP_target_function_construct();
		
		//return 0;
		IloEnv env;

	 try {
	   const char* filename = "LP_StandardInput.txt";
	   IloNumVar::Type type = ILOFLOAT;

	  
	  IloNumArray objCoef(env);
	  IloNumArray2 fcoef(env);
	  IloNumArray fb(env);

	  ifstream file(filename);
	  if(!file)
	  {
		cerr<<"ERROR: could not open file '"<<filename
			<<"' for reading"<<endl;
		return -1;
	  }

	  file >> fcoef >> fb;
	  file.close();
	  file.open("objectiveCoef.txt");
	  file>>objCoef;

	  IloInt fNum = objCoef.getSize();
	  IloInt rNum = fcoef.getSize();

	  IloModel mod(env);
	  IloNumVarArray fOptimal(env,fNum,0.0,1.0,type);
	  mod.add(IloMinimize(env,IloScalProd(fOptimal,objCoef)));

	  for(IloInt i = 0; i < rNum; i++)
	  {
		IloExpr expr(env);
		for(IloInt j = 0; j <fcoef[i].getSize(); j++)
		{
		  expr += fOptimal[j] * fcoef[i][j];
		}
		mod.add(expr <= fb[i]);
		expr.end();
		
	  }
	  //solve model
	  IloCplex cplex(mod);
	  cplex.solve();
	  float *f = new float[N];
	  for(IloInt i = 0; i < N; i++)
	  {
		  f[i] =  cplex.getValue(fOptimal[i]);
	  }
	  bool wError = false;
	  for(int i = 0; i < N; i++)
	  {
		 if(f[i] > ESP)
		 {
		   wError = true;
		   BEC += 1;
		 }
	  } 
	  delete[] f;
	  if(wError)
	  {
		WEC += 1;
	  }

	  dCount++;
	  if(WEC >= 50)
	  {
		
		ofstream out("result.txt",ofstream::app);
		out<<"p:"<<setprecision(1)<<setiosflags(ios::fixed)<<p<<endl<<setprecision(6);
		out<<"Count: "<<dCount << " , bit error : " << BEC << " , word error : " << WEC << endl;
		out<<"BER : " << (float)BEC/(dCount*N) << ",WER : " << (float)WEC/dCount <<endl<<endl;
		out.close();

		p += snrInterval;
		if(p > snrEnd)
		{
		  break;
		}
		else
		{
		   dCount = 0;
		   BEC = 0;
		   WEC = 0;
		}
	 }   
	 }catch (IloException& ex) {
		  cerr << "Error: " << ex << endl;
	   }
	   catch (...) {
		  cerr << "Error" << endl;
	   }
	  env.end();
  }//while true
  system("pause");
  return 0;
}

void read_H()
{
	int i,j,temp;
	int *temp_p;

	FILE *fp;
	fp = fopen("alist_arr.txt", "r");
	if (fp==NULL)
	{
		printf("open file alist_arr.txt error!");
		system("pause");exit(0);
	}

	fscanf(fp, "%d", &N);
	fscanf(fp, "%d", &M);

	fscanf(fp, "%d", &MaxColDegree);
	fscanf(fp, "%d", &MaxRowDegree);

	pColDegree = (int *)malloc(N*sizeof(int));
	for (temp_p=pColDegree; temp_p<pColDegree+N; temp_p++)
	{
		fscanf(fp, "%d", temp_p);	//��ÿһ����1�ĸ���������pColDegree[]��
	}

	pRowDegree = (int *)malloc(M*sizeof(int));
	for (temp_p=pRowDegree; temp_p<pRowDegree+M; temp_p++)
	{
		fscanf(fp,"%d",temp_p);//��ÿһ����1�ĸ���������pRowDegree[]��
	}

	H = (int **)malloc(M*sizeof(int *));
	for (i=0; i<M; i++)
	{
		H[i] = (int *)malloc(N*sizeof(int));
	}//����һ��M*N��H[]����,�洢�������

	for (i=0; i<M; i++)
	{
		for (j=0; j<N; j++)
		{
			H[i][j] = 0;
		}
	}

	A_V = (int **)malloc(N*sizeof(int *));
	for (i=0; i<N; i++)
	{
		A_V[i] = (int *)malloc(MaxColDegree*sizeof(int));
	}//A_V[]����洢ÿ��Ϊ1���������±꣬�������Ǵ�1��ʼ�������ʱ��Ҫ��1

	for (i=0; i<N; i++)
	{
		for (j=0; j<MaxColDegree; j++)
		{
			fscanf(fp, "%d", &temp);
			A_V[i][j] = temp-1;
			if (temp>0)
			{
				H[temp-1][i] = 1;
			}
		}
	}//H[]����������

	B_C = (int **)malloc(M*sizeof(int *));
	for (i=0; i<M; i++)
	{
		B_C[i] = (int *)malloc(MaxRowDegree*sizeof(int));
	}

	for (i=0; i<M; i++)
	{
		for (j=0; j<MaxRowDegree; j++)
		{
			fscanf(fp, "%d", &temp);
			B_C[i][j] = temp-1;
		}
	}//B_C[]����洢ÿ��Ϊ1���������±꣬�������Ǵ�1��ʼ�������ʱ��Ҫ��1

	fclose(fp);

	print_H();  //��H�������뵽�ļ���
	
}

void print_H()
{
	int i,j;
	FILE *fp;
	fp = fopen("H_file.txt","w");
	if (fp==NULL)
	{
		printf("open file H_file.txt error!");
		system("pause");exit(0);
	}
	for (i=0; i<M; i++)
	{
		for (j=0; j<N; j++)
		{
			fprintf(fp, "%d", H[i][j]);
			putc(' ',fp);
		}
		putc('\n',fp);
	}

	fclose(fp);	
}

void read_code()//������
{
	codeword = (int *)malloc(N*sizeof(int));//�洢��������

#if READ_CODE//���ⲿ������

	int *temp;
	FILE *fp;
	fp = fopen("code_file.txt","r");
	if (fp==NULL)
	{
		printf("open file code_file.txt error!");
		system("pause");exit(0);
	}
	
	for (temp=codeword; temp<codeword+N; temp++)
	{
		fscanf(fp, "%d", temp);
	}

/*	int i;
	for (i=0;i<N;i++)
	{
		printf("%d  ",codeword[i]);
	}
*/
	fclose(fp);

#else 
	
	for (int i=0; i<N; i++)
	{
		codeword[i] = 0;
	}

#endif
}

void LP_target_function_construct()//���ݽ������֣�����Ŀ�꺯��
{
	read_code();  //��������
	int i;
	double *y_signal;
	y_signal = (double *)malloc(N*sizeof(double));

	for (i=0; i<N; i++)
	{
		float  temp = ((float)(rand()%100))/100;
		if(temp - p < ESP)
		{
			y_signal[i] = (codeword[i]+1)%2;
		}
		else if (temp - p >= ESP)
		{
			y_signal[i] = codeword[i];
		}
		else printf("code word error!");		//δ���б��룡����
	}//�����������y
	FILE *fp;
	fp = fopen("objectiveCoef.txt","w");

	fprintf(fp,"[");
	for (i=0; i<N-1; i++)
	{
		fprintf(fp,"%d , ",y_signal[i]);
	}
	fprintf(fp,"%d ]\n",y_signal[N-1]);

	fclose(fp);
	free(y_signal);
	free(codeword);

}

void setArrayToZero(int *array, int len)
{
	int i;
	for( i = 0; i < len; i++)
	{
	  array[i] =0;
	}
}
void LP_matrix_construct()
{
	
	int i,j,k,s;
	int rCount;
	wTofMatrix *WF;
	int* R;
	int** coef;
	FILE *fp;


	fp = fopen("LP_StandardInput.txt","w");

	WF = wTof(MaxRowDegree);//���ɾ���A��B
	rCount = M*WF->rNum;
	coef = (int**) malloc(rCount*sizeof(int*));
	for( i = 0; i < rCount; i++)
	{
	   coef[i] = (int*) malloc (N * sizeof(int));
	   for( k = 0; k < N; k++)
	   {
	      coef[i][k] = 0;
	   }
	}//coef[]:rCount*N

	
	for( i = 0; i < M; i++)
	{
		for (j = 0; j < WF->rNum; j++)
		{
		   for ( k = 0; k < WF->cNum; k++)
		   {
			  coef[i*WF->rNum + j][B_C[i][k]] = (WF->A)[j][k];
		   }
		}
	}
	
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
	for ( i = 0; i < M; i++)
	{
		for ( j = 0; j < WF->rNum; j++)
		{
			fprintf(fp, (i == M-1) && (j==WF->rNum-1) ? "%d " : "%d , ",(WF->B)[j]);
		}
	}
	fprintf(fp,"]\n");
	fclose(fp);

	freeWtoFMatrix(WF);
	for ( i = 0; i < rCount; i++)
	{
	  free(coef[i]);
	}
	free(coef);
}

void enum_comb(int c_node, int k, int n, int j_begin, int j_end, w_Enum* wComb)
			  //�ڵ�c_node��У��ڵ�������n�������ڵ��У�ö�ٳ����� C(n, k) ����ϣ�
			  //����ṹ wComb[j_begin] �� wComb[j_end] ������vec��

//void enum_comb(int c_node, int k, int n, int j_begin, int j_end, vector<vector<int> >& vec)
			  //�ڵ�c_node��У��ڵ�������n�������ڵ��У�ö�ٳ����� C(n, k) ����ϣ�
			  //�������� vec[j_begin] �� ���� vec[j_end] ��
//enum_comb(i,k,pRowDegree[i],j_begin,j_end,w[i]);

{

	
/*	if( (j_end-j_begin+1) != combination(k,n) )		// ����ע�͵��Խ�Լ����ʱ��
	{
		printf("input function(enum_comb) variable error! 1");
		system("pause");exit(0);
	}
*/	
	
	int i,j;
	int max_sign_i,assign;
	int *arr = (int *)malloc(pRowDegree[c_node]*sizeof(int));
	for (i=0; i<pRowDegree[c_node]; i++)
	{
		arr[i] = B_C[c_node][i];
	}
	
	if ( k==0 || k==n )
	{
		if (j_begin != j_end)
		{
			printf("input function(enum_comb) variable error! ");
			system("pause");exit(0);
		}
	}

	if (k==0)
	{
		wComb[j_begin].vecSize = 1;
		wComb[j_begin].vec = (int *)malloc(1*sizeof(int));
		wComb[j_begin].vec[0] = -1;
	//	wComb[j_begin].vec.push_back(-1);	// -1����ռ�
	//	vec[j_begin].push_back(-1);	// -1����ռ�
	}
	if (k==n)
	{
/*�½�ƽע�� 2011-10-11		
		for (i=0; i<pRowDegree[c_node]; i++)
		{
*/
			wComb[j_begin].vecSize = k;
			wComb[j_begin].vec = (int *)malloc(k*sizeof(int));
			for (j=0; j<k; j++)
			{
				wComb[j_begin].vec[j] = arr[j];
			}
			//wComb[j_begin].vec.push_back(arr[i]);
			//vec[j_begin].push_back(arr[i]);
/*�½�ƽע�� 2011-10-11				
		}
*/
	}
	else if ( k>0 && k<n )
	{
		int *c = (int *)malloc(k*sizeof(int));
		for (i=0; i<k; i++)
		{
			c[i] = i;
		}
		
		for (assign=j_begin; assign<=j_end; assign++)
		{
			wComb[assign].vecSize = k;
			wComb[assign].vec = (int *)malloc(k*sizeof(int));

			for (i=0; i<k; i++)
			{
				wComb[assign].vec[i] = arr[c[i]];
			//	wComb[assign].vec.push_back(arr[c[i]]);
			//	vec[assign].push_back(arr[c[i]]);
			}
			
			max_sign_i = -1;
			for (i=0; i<k; i++)
			{
				if( c[i] < n-k+i ) 
					max_sign_i = i;
			}
			if (max_sign_i == -1)
			{
				if (assign != j_end)
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
		free(c);
	}
	free(arr);
}


int comb_even_number_sum(int n)		//����ÿ��У��ڵ�C(0,n)+C(2,n)+C(4,n)...
{
	int i;
	int sum = 0;
	for (i=0; i<=n; i++)
	{
		if ( i % 2 == 0 )
		{
			sum += combination(n , i);
		}
	}
	return sum;
}

int comb_odd_number_sum(int n)		//������1������������ĺͣ�����ÿ��У��ڵ�C(1,n)+C(3,n)+C(5,n)...
{
	int i;
	int sum = 0;
	for (i=0; i<=n; i++)
	{
		if ( i % 2 != 0 )
		{
			sum += combination(n , i);
		}
	}
	return sum;
}

long factorial(int n)		//����n�Ľ׳�
{
	if(n>16)
	{
		printf("can't compute the factorial!");
		system("pause");exit(0);
	}
	long result=1;
/*	if(n==1)
		result=1;
	else
		result*=factorial(n-1);
*/
	while(n>0)
	{
		result*=n;
		n--;
	}
	return result;
}

long combination(int n, int i)		//���������C(n, i)
{
	if (i>n)
	{
		printf("the combination input error!");
		system("pause");exit(0);
	}
	if (i==0 || i==n)
	{
		return 1;
	}
	else if (i==1)
	{
		return n;
	}
	else 
	{
		if (i > (n/2))   // C(n,i) = C(n,(n-i))
		{
			i = n-i;
		}

		//��򵥵��㷨������ n!/i!*(n-i)! �����ǳ��������
		//return factorial(n)/(factorial(i)*factorial(n-i));

		//�ڶ����㷨������ n*(n-1)*(n-2)*...*(n-i+1) / i! ������Լ��
		int *numerator; //�洢���ӵ�����
		int *denominator; //�洢��ĸ������
		numerator = (int *)malloc(i*sizeof(int));
		denominator = (int *)malloc((i-1)*sizeof(int));

		int k,s;
		for (k=0; k<i; k++)
		{
			numerator[k] = n-k; 
		}//numerator[]�洢n��n-1...n-i+1
		for (k=0; k<i-1; k++)
		{
			denominator[k] = i-k;
		}//denominator[]�洢i��i-1...2 ,s[0]=i...s[i-2]=2
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
		}//�ӷ��ӵ����λ���ͷ�ĸ�������Լ�֣�����Լ�֣�����Ӧ���ӳ��Է�ĸ����ĸΪ1
		
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

wTofMatrix* wTof(int rDeg)//����dc(m)
{
  int** A;
  int* B;
  int i,j,k;
  int assign;
  int combNum, offset = 0;
  int wCount = 0;
  int max_sign_i;
  
  wCount = comb_odd_number_sum(rDeg);//������1�����������
  
  A = (int**)malloc((wCount)*sizeof(int*));
  for(i = 0; i < wCount; i++)
  {
     A[i] = (int*) malloc (rDeg*sizeof(int));
  }//A[ wCount][rDeg]
  
  for( i = 0 ; i < wCount; i++)
  {
	 for( j = 0; j < rDeg; j++)
     {
        A[i][j] = -1;	
     } 
  }

  B = (int*) malloc (wCount*sizeof(int));
  
  k = 1;
  while( k <= rDeg )
  { 
  	    int *c = (int *)malloc( k*sizeof(int) );
		for (i=0; i<k; i++)
		{
			c[i] = i;
		}
		
		
		combNum = combination( rDeg , k );

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
				if( c[i] < rDeg-k+i ) 
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
 
  /*
  //print A,B
  printf("A:\n");
  for ( i = 0; i <  wCount; i++)
  {
     for ( j = 0; j < rDeg; j++)
     {
       printf("%2d ", A[i][j]);	
     }	
      printf("\n");
  }
  printf("B:\n");
  for( i = 0; i < wCount; i++)
  {
    printf("%2d ",B[i]);
  }
  printf("\n");
  */
  wTofMatrix *result = (wTofMatrix*) malloc (sizeof(wTofMatrix));
  result->A = A;
  result->B = B;
  result->rNum = wCount;
  result->cNum = rDeg;
  
  return result;	
}

void freeWtoFMatrix(wTofMatrix* w)
{
   int i;
   for( i = 0; i < w->rNum; i++)
   {
     free((w->A)[i]);
   }
   free(w->A);
   free(w->B);
   free(w);
}

void freeABH()
{
  int i , j;
  for (i = 0; i < N; i++)
  {
    free(A_V[i]);
  }
  free(A_V);

  for (j = 0; j < M ; j++)
  {
    free(H[j]);
	free(B_C[j]);
  }
  free(H);
  free(B_C);
}