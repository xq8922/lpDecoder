#ifndef LINEARPROGRAMMING_H
#define LINEARPROGRAMMING_H


#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;


//���Թ滮���õ��Ĵ洢�ṹ	 
struct LP_PIVOT_STRUCT 
{
	int* Nonbasic;
	int* BasicVar;
	double **A;
	double *b;
	double *c;
	double v;
};


// ��������
void MemoryAllocateLP(LP_PIVOT_STRUCT *LP_STR, int NonbasicSize);
void MemoryFreeLP(LP_PIVOT_STRUCT *LP_STR);
void MemoryCopy(LP_PIVOT_STRUCT *dest, LP_PIVOT_STRUCT *src, int NonbasicSize);
void PIVOT(LP_PIVOT_STRUCT* LPstruct, int l, int e, int NonbasicSize);
void PIVOT_for_INITIALIZE(LP_PIVOT_STRUCT* LPstruct, int l, int e, int NonbasicSize);
double LinearProgrammingWithStandardInput();



#endif