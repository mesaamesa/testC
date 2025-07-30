/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa												*/
/*																		*/
/*	Versi¢n : 															*/
/*		   		v 1.0		1-XI-1994			10-II-1995				*/
/*																		*/
/************************************************************************/

#include "nuevo.h"
#include "errormat.h"
#include "iomatriz.h"
#include "opmatriz.h"
#include "fcmatriz.h"
#include <math.h>
/*traza*/
#include <stdio.h>

void VEvepnRe(matriz X, matriz A, matriz B, vector L, unsigned int Dim)
{
	indice i,j;
	unsigned int iter,itermax;
    escalar errmed,normx;
    vector x,y,t,v;
    matriz Ad,M,U;

    IOmemvect(&x,Dim);
    IOmemvect(&y,Dim);
    IOmemvect(&t,Dim);
    IOmemvect(&v,Dim);
    IOmemmat(&Ad,Dim);
    IOmemmat(&M,Dim);
    IOmemmat(&U,Dim);

    itermax = KITERMAX * Dim;
    iter = 0;
    for(i = 0;i < Dim;i++)
	{
		if(L[i] == INFINITO)
			OPcerovect(x,Dim);
		else
        {
		   	iter = 0;
    	    errmed = COTA;
    		OPctemat(Ad,B,L[i],Dim);
	        OPrestamat(Ad,A,Ad,Dim);
    	    FClu(M,U,Ad,Dim);
			OPunovect(y,Dim);
	        while((errmed >= COTA) && (iter < itermax))
    	    {
        		iter++;
	            OPmatvect(t,B,y,Dim);
    	        OPsisteqinf(M,v,t,Dim);
        	    OPsisteqsup(U,x,v,Dim);
            	normx = OPnorma2(x,Dim);
	            OPctevect(x,x,(1 / normx),Dim);
    	        for(j = 0;j < Dim;j++)
        	    	t[j] = fabs(x[j]) - fabs(y[j]);
				errmed = OPnormamax(t,Dim);
				OPcpyvect(y,x,Dim);
    	    }
        }
        if(iter == itermax)
			OPcerovect(x,Dim);
        for(j = 0;j < Dim;j++)
        	X[j][i] = x[j];
    }

    IOfreevect(x,Dim);
    IOfreevect(y,Dim);
    IOfreevect(t,Dim);
    IOfreevect(v,Dim);
    IOfreemat(Ad,Dim);
    IOfreemat(M,Dim);
    IOfreemat(U,Dim);
}

void VEvep1Re(vector X, matriz A, matriz B, escalar L, unsigned int Dim)
{
	indice i;
	unsigned int iter,itermax;
    escalar errmed,normx;
    vector x,y,t,v;
    matriz Ad,M,U;

    IOmemvect(&x,Dim);
    IOmemvect(&y,Dim);
    IOmemvect(&t,Dim);
    IOmemvect(&v,Dim);
    IOmemmat(&Ad,Dim);
    IOmemmat(&M,Dim);
    IOmemmat(&U,Dim);

	itermax = KITERMAX * Dim;
	iter = 0;
	if(L == INFINITO)
		OPcerovect(x,Dim);
	else
	{
		errmed = COTA;
   		OPctemat(Ad,B,L,Dim);
		OPrestamat(Ad,A,Ad,Dim);
		FClu(M,U,Ad,Dim);
		OPunovect(y,Dim);
    	while((errmed >= COTA) && (iter < itermax))
	    {
    	 	iter++;
        	OPmatvect(t,B,y,Dim);
	        OPsisteqinf(M,v,t,Dim);
    	    OPsisteqsup(U,x,v,Dim);
        	normx = OPnorma2(x,Dim);
	        OPctevect(x,x,(1 / normx),Dim);
    	    for(i = 0;i < Dim;i++)
        		t[i] = fabs(x[i]) - fabs(y[i]);
			errmed = OPnormamax(t,Dim);
			OPcpyvect(y,x,Dim);
	    }
    }
    if(iter == itermax)
		OPcerovect(x,Dim);
    for(i = 0;i < Dim;i++)
   	   	X[i] = x[i];

    IOfreevect(x,Dim);
    IOfreevect(y,Dim);
    IOfreevect(t,Dim);
    IOfreevect(v,Dim);
    IOfreemat(Ad,Dim);
    IOfreemat(M,Dim);
    IOfreemat(U,Dim);
}
