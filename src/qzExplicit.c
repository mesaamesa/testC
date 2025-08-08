/**
    Autor : Antonio Mesa
    Version : v 1.0
        1-XI-1994
        10-II-1995
**/

#include "defines.h"

void QZexplicito(matriz A, matriz B, unsigned int Dim)
{
	indice i,k;
	escalar sigma,epsA,h,c,s;
	matriz T;

	IOmemmat(&T,Dim);

	for(k = 0;k < (Dim-1);k++)
	{
		epsA = EPS * fabs(A[Dim-1-k][Dim-1-k] + A[Dim-2-k][Dim-2-k]);
		while(fabs(A[Dim-1-k][Dim-2-k]) >= epsA)
		{
			sigma = A[Dim-1-k][Dim-1-k] / B[Dim-1-k][Dim-1-k];
			OPctemat(T,B,sigma,Dim-k);
			OPrestamat(A,A,T,Dim-k);
			OPidentidad(T,Dim);
			for(i = 1;i < (Dim-k);i++)
			{
				h = sqrt(A[i-1][i-1] * A[i-1][i-1] + A[i][i-1] * A[i][i-1]);
				c = A[i-1][i-1] / h;
				s = A[i][i-1] / h;
                OPgivenspre(A,i-1,i,c,s,Dim-k);
                OPgivenspre(B,i-1,i,c,s,Dim-k);
				h = sqrt(B[i][i] * B[i][i] + B[i][i-1] * B[i][i-1]);
				c = B[i][i] / h;
				s = B[i][i-1] / h;
                OPgivenspost(B,i-1,i,c,s,Dim-k);
                OPgivenspost(T,i-1,i,c,s,Dim-k);
			}
			OPprodmat(A,A,T,Dim-k);
			OPctemat(T,B,sigma,Dim-k);
			OPsumamat(A,A,T,Dim-k);
		}
	}
	IOfreemat(T,Dim);
}

void QZvepslu(matriz X, matriz A, matriz B, vector L, unsigned int Dim)
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
    for(i = 0;i < Dim;i++)
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
        if(iter == itermax)
			OPcerovect(x,Dim);
        for(j = 0;j < Dim;j++)
        	X[j][i] = x[j];
    }

    IOfreevect(x,Dim);
    IOfreevect(y,Dim);
    IOfreevect(t,Dim);
    IOfreemat(Ad,Dim);
    IOfreemat(M,Dim);
    IOfreemat(U,Dim);
}

void QZveplu1(vector X, matriz A, matriz B, escalar L, unsigned int Dim)
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
    if(iter == itermax)
		OPcerovect(x,Dim);
    for(i = 0;i < Dim;i++)
       	X[i] = x[i];

    IOfreevect(x,Dim);
    IOfreevect(y,Dim);
    IOfreevect(t,Dim);
    IOfreemat(Ad,Dim);
    IOfreemat(M,Dim);
    IOfreemat(U,Dim);
}

			        	


    	
	