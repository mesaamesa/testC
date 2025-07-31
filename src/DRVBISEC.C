/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa												*/
/*																		*/
/*	Versi¢n : 															*/
/*		   		v 1.0		1-XI-1994			10-II-1995				*/
/*																		*/
/************************************************************************/

#include "nuevo.h"
#include "bisec.h"
#include "iomatriz.h"
#include "fcmatriz.h"
#include "opmatriz.h"
#include <stdio.h>

void DRdriverbisec(string matA, string matB, escalar S)
{
	unsigned int Dim,k;
    escalar lambda,a,b;
	matriz A,B,R,Q;


	IOrdmat(&A,matA,&Dim);
	IOrdmat(&B,matB,&Dim);

	IOmemmat(&R,Dim);
	IOmemmat(&Q,Dim);

	FCqrgram(Q,R,B,Dim);
	IOfreemat(B,Dim);
	OPtransmat(Q,Dim);
	OPprodmat(A,Q,A,Dim);
    IOfreemat(Q,Dim);

	FChesstriorto2(A,R,Dim);

	OPsubdiaghess(A,R,Dim);

	BIbuscainterval(&k,&a,&b,A,R,S,Dim);
/*traza*/printf("\nvap = %d \na0 = %g  \nb0 = %g",k,a,b);getchar();
	BIbisec(&lambda,a,b,k,A,R,Dim);

	printf("\nlambda = %1.17g",lambda);
	IOfreemat(A,Dim);
    IOfreemat(R,Dim);
}