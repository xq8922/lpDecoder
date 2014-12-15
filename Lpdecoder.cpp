#include<iostream>
//#include <ilcplex>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include "random.h"
using namespace std;
//ILOSTLBEGIN
int N;//��������
int M;//��������

int MaxColDegree;//�еĶ�
int MaxRowDegree;//�еĶ�

int *initcodeword;
int *receivercode;

int **H;//У�����

int **A_V; //��Ϣ�ڵ�
int **B_C; //У��ڵ�

void read_file();
void print_H();

double var;

int **coef;
//���ļ������ļ��е����ݶ������
void read_file()
{
	int i,j,temp;

	FILE *fp;
	fp = fopen("alist_arr.txt", "r");
	if (fp==NULL)
	{
		printf("open file alist_arr.txt error!");
	}

	fscanf(fp, "%d", &N);
	fscanf(fp, "%d", &M);

	fscanf(fp, "%d", &MaxColDegree);
	fscanf(fp, "%d", &MaxRowDegree);

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




	//���ж�У�����������
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



	//���ж�У����������䣬����䲽��ʡ�ԣ���Ϊ������������Ѿ����ж�У�������������
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
		system("pause");
		exit(0);
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


//Ŀ�꺯���Ĺ��죬��Ŀ�꺯���о������Ƽ���������
void function_construct()
{
	int i;
	initcodeword = (int *)malloc(N*sizeof(int));//�洢��������
	receivercode = (int *)malloc(N*sizeof(int));
	for (i=0; i<N; i++)
	{
		initcodeword[i]=receivercode[i]=0;
	}
	for(i=0;i<N;i++)
	{
		receivercode[i]=pow(-1, initcodeword[i]);
	}
	for(i=0; i<N; i++) 
	{
		receivercode[i]+=Gauss() * var ;
	    receivercode[i]=2 * receivercode[i] / (var * var);
	}
}


//����n�Ľ׳�
long factorial(int n)
{
	int result=1;
	while(n>0)
	{
		result*=n;
		n--;
	}
	return result;
}



//���������c(n,i)
long combination(int n,int i)
{
	if(n<i)
	{
		cout<<"the combination input error!"<<endl;
	    return -1;
	}
	if(i==0||i==n)
	{
		return 1;
	}
	else if(i==1)
	{
		return n;
	}
	else
	{
		return factorial(n)/(factorial(i)*factorial(n-i));
	}

}


//���������c(n,1)+c(n,3)+c(n,5)..............
int comb_odd_number(int n)
{
	int i;
	int sum=0;
	for(i=0;i<n;i++)
	{
		if(i%2!=0)
			sum+=combination(n,i);
	}
	return sum;
}

//LP��Լ��������������1�ľ���Ĺ���
void LP_matrix_construct()
{

	int **A;
	int *B;
	int wCount=0,rCount,i,j,k,offset=0,combNum,assign;
	int max_sign_i;

	wCount=comb_odd_number(MaxRowDegree);

	A = (int**)malloc((wCount)*sizeof(int*));
	for(i = 0; i < wCount; i++)
	{
		A[i] = (int*) malloc (MaxRowDegree*sizeof(int));
	}
	
	for( i = 0 ; i < wCount; i++)
	{
		for( j = 0; j < MaxRowDegree; j++)
		{
			A[i][j] = -1;	
		}
	}

	B = (int*) malloc (wCount*sizeof(int));
	
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

	cout<<"�������A:"<<endl;
	for ( i = 0; i <  wCount; i++)
	{
		for ( j = 0; j < MaxRowDegree; j++)
		{
			cout<<A[i][j]<<" ";
		}
		cout<<endl;
	}

	cout<<"�������B��"<<endl;
	for( i = 0; i < wCount; i++)
	{
		cout<<B[i]<<" ";
	}//�����һ��С����
	cout<<endl;



	//��ʼ����Լ����������Ҫ�ľ���coef
	rCount = M*wCount;
	coef = (int**) malloc(rCount*sizeof(int*));
	for( i = 0; i < rCount; i++)
	{
	   coef[i] = (int*) malloc (N * sizeof(int));
	   for( k = 0; k < N; k++)
	   {
	      coef[i][k] = 0;
	   }
	}
	for( i = 0; i < M; i++)
	{
		for (j = 0; j < wCount; j++)
		{
		   for ( k = 0; k < MaxRowDegree; k++)
		   {
			  coef[i*wCount + j][B_C[i][k]] = A[j][k];
		   }
		}
	}

/*
	cout<<"["<<endl;
	for ( i = 0; i < rCount; i ++)
	{
		cout<<"[";
		for ( j = 0; j < N; j++)
		{
			cout<<coef[i][j]<<" ";
		}
		if(i==rCount-1)
			cout<<"]\n";
		else
			cout<<"] ,\n";
	}
	cout<<"]\n";
	

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
	cout<<endl;
}

void main()
{


	read_file();

	srand((unsigned)time(NULL));
	LP_matrix_construct();

	for(double SNR=2.0;SNR<5.1;SNR+=0.1)
	{
		double No=0.0;
	//	cout<<SNR<<endl;
		No=pow(10,SNR/10);
		
		var=sqrt(1.0/((2.0*(N-M)/N)*No));
		
		while(true)
		{
			function_construct();
		}
	}
}
