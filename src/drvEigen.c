/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa												*/
/*																		*/
/*	Versi�n : 															*/
/*		   		v 1.0		1-XI-1994			10-II-1995				*/
/*																		*/
/************************************************************************/

#include <stdlib.h>
#ifdef _DEBUG_
#include <stdio.h>
#endif

#include "defines.h"
#include "io.h"
#include "operator.h"
#include "convert.h"
#include "qzExplicit.h"

void DRdrivereig(string matA, string matB, bool flag, escalar S)
{
	indice i;
	unsigned int Dim;
	vector L,Y;
	matriz A,B,R,Q,X;

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOrdmat(&A,matA,&Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOrdmat(&B,matB,&Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOmemmat(&R,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOmemmat(&Q,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	FCqrgram(Q,R,B,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOfreemat(B,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	OPtransmat(Q,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	OPprodmat(A,Q,A,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOfreemat(Q,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	FChesstriorto2(A,R,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	QZexplicito(A,R,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOmemvect(&L,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	for(i = 0;i < Dim;i++)
	{
		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		L[i] = A[i][i] / R[i][i];
	}

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	OPsort(L,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOfreemat(A,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOfreemat(R,Dim),

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOrdmat(&A,matA,&Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOrdmat(&B,matB,&Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	if(flag)
	{
		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		i = 0;

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		while((L[i] < S) && (i < Dim))
		{
			#ifdef _DEBUG_
			printf("\nfile:%s (%d)",__FILE__,__LINE__);
			#endif
			i++;
		}

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		if(i == Dim)
		{
			#ifdef _DEBUG_
			printf("\nfile:%s (%d)",__FILE__,__LINE__);
			#endif
			L[0] = L[Dim-1];

			#ifdef _DEBUG_
			printf("\nfile:%s (%d)",__FILE__,__LINE__);
			#endif
		}
		else
		{
			#ifdef _DEBUG_
			printf("\nfile:%s (%d)",__FILE__,__LINE__);
			#endif
			if((i > 0) && ((L[i] - S) > (S - L[i-1])))
			{
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				L[0] = L[i-1];
			}

			#ifdef _DEBUG_
			printf("\nfile:%s (%d)",__FILE__,__LINE__);
			#endif
			if((i > 0) && ((S - L[i-1]) > (L[i] - S)))
			{
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				L[0] = L[i];
			}
		}

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		IOmemvect(&Y,Dim);

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		QZveplu1(Y,A,B,L[0],Dim);

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		IOwrvect(Y,"vep00",Dim);

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		IOwrvect(L,"vap00",1);

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		IOfreevect(Y,Dim);
	}
	else
	{
		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		IOmemmat(&X,Dim);

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		QZvepslu(X,A,B,L,Dim);

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		IOwrpar(X,L,"vep","vap",Dim);

		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		IOfreemat(X,Dim);
	}

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOfreemat(A,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOfreemat(B,Dim);

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	IOfreevect(L,Dim);
}
