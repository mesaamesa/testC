/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa												*/
/*																		*/
/*	Versi�n : 															*/
/*		   		v 1.0		1-XI-1994			10-II-1995				*/
/*																		*/
/************************************************************************/

#include "defines.h"
#include "io.h"
#include "convert.h"
#include "operator.h"
#include "bisection.h"

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
	BIbisec(&lambda,a,b,k,A,R,Dim);

	IOfreemat(A,Dim);
    IOfreemat(R,Dim);
}