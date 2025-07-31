#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include "nuevo.h"
#include "iomatriz.h"
#include "opmatriz.h"
#include "fcmatriz.h"
#include "qzimplic.h"

void DRdriver(string matA, string matB)
{
	unsigned int i,j;
	unsigned int dim;
	vector L;
	matriz A,B,R,Q,X;


	IOrdmat(&A,matA,&dim);
	IOrdmat(&B,matB,&dim);

	IOmemmat(&R,dim);
	IOmemmat(&Q,dim);


/*traza*/puts("\nQR matriz B");
	FCqrgram(Q,R,B,dim);
	IOfreemat(B,dim);
	OPtransmat(Q,dim);
	OPprodmat(A,Q,A,dim);
    IOfreemat(Q,dim);

/*traza*/puts("\nHESS-TRI matrices A-B");
	FChesstriorto2(A,R,dim);

/*traza*/puts("\nQZ");
/*	EIqzdesp2(A,R,dim);*/
	while(getchar() != 'x')
    {
    	OPfiltramat(A,EPS*OPnormamaxmat(A,dim),dim);
        OPfiltramat(R,EPS*OPnormamaxmat(R,dim),dim);
    	EIqzimplicito(A,R,dim);
/*traza*/printf("\nA%d %d = %g",dim-1,dim-2,A[dim-1][dim-2]);
	}

/*traza*/puts("\nVAP's");
	IOmemvect(&L,dim);
	for(i = 0;i < dim;i++)
		L[i] = A[i][i] / R[i][i];
/*traza*/IOputvect(L,"l",'l',dim);

    IOfreemat(A,dim);
    IOfreemat(R,dim),
	IOfreevect(L,dim);
}