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

void QZstep(matriz A, matriz B,	matriz Q, matriz Z, escalar eps,
			escalar epsA, escalar epsB, indice inf, indice sup, unsigned int Dim)
{
	indice i;
	escalar sigma,h,c,s;
	matriz T,TT;

	IOmemmat(&T,sup);

	sigma = A[sup-1][sup-1] / B[sup-1][sup-1];
	OPctesubmat(T,B,sigma,inf,sup);
	OPrestasubmat(A,A,T,inf,sup);
	OPidentidad(T,sup);
	for(i = (inf+1);i < sup;i++)
	{
		if(fabs(A[i][i-1]) > eps * (fabs(A[i][i]) + fabs(A[i-1][i-1])))
		{
			h = sqrt(A[i-1][i-1] * A[i-1][i-1] + A[i][i-1] * A[i][i-1]);
			c = A[i-1][i-1] / h;
			s = A[i][i-1] / h;
			OPgivenspre(A,i-1,i,i-1,c,s,epsA,sup);
			OPgivenspre(B,i-1,i,i-1,c,s,epsB,sup);
			OPgivenspre(Q,i-1,i,0,c,s,eps,Dim);
		}
		if(fabs(B[i][i-1]) > eps * (fabs(B[i][i]) + fabs(B[i-1][i-1])))
		{
			h = sqrt(B[i][i] * B[i][i] + B[i][i-1] * B[i][i-1]);
			c = B[i][i] / h;
			s = B[i][i-1] / h;
			OPgivenspost(B,i-1,i,inf,c,s,epsB,sup);
			OPgivenspost(T,i-1,i,inf,c,s,epsB,sup);
			OPgivenspost(Z,i-1,i,0,c,s,eps,Dim);
		}
	}
	IOmemmat(&TT,sup);
	OPprodsubmat(TT,A,T,inf,sup);
	OPctesubmat(T,B,sigma,inf,sup);
	OPsumasubmat(A,TT,T,inf,sup);

	IOfreemat(TT,sup);
	IOfreemat(T,sup);
}

void QZstep2(matriz A, matriz B, escalar eps, escalar epsA, escalar epsB, indice inf, indice sup)
{
	indice i;
	escalar sigma,h,c,s;
	matriz T,TT;

	IOmemmat(&T,sup);

	sigma = A[sup-1][sup-1] / B[sup-1][sup-1];
	OPctesubmat(T,B,sigma,inf,sup);
	OPrestasubmat(A,A,T,inf,sup);
	OPidentidad(T,sup);
	for(i = (inf+1);i < sup;i++)
	{
		if(fabs(A[i][i-1]) > eps * (fabs(A[i][i]) + fabs(A[i-1][i-1])))
        {
			h = sqrt(A[i-1][i-1] * A[i-1][i-1] + A[i][i-1] * A[i][i-1]);
			c = A[i-1][i-1] / h;
			s = A[i][i-1] / h;
			OPgivenspre(A,i-1,i,i-1,c,s,epsA,sup);
			OPgivenspre(B,i-1,i,i-1,c,s,epsB,sup);
		}
		if(fabs(B[i][i-1]) > eps * (fabs(B[i][i]) + fabs(B[i-1][i-1])))
		{
			h = sqrt(B[i][i] * B[i][i] + B[i][i-1] * B[i][i-1]);
			c = B[i][i] / h;
			s = B[i][i-1] / h;
			OPgivenspost(B,i-1,i,inf,c,s,epsB,sup);
			OPgivenspost(T,i-1,i,inf,c,s,epsB,sup);
		}
	}
	IOmemmat(&TT,sup);
	OPprodsubmat(TT,A,T,inf,sup);
	OPctesubmat(T,B,sigma,inf,sup);
	OPsumasubmat(A,TT,T,inf,sup);

	IOfreemat(TT,sup);
	IOfreemat(T,sup);
}

void QZdeflaccion(matriz A, matriz B, matriz Q, matriz Z, indice k,
			escalar eps, escalar epsA, escalar epsB, indice inf, indice sup, unsigned int Dim)
{
	indice i;
    escalar h,c,s;

	if(k != inf)
	{    	

		for(i = k;i < (sup-1);i++)
		{
			if(fabs(B[i+1][i+1]) > epsB)
			{
				h = sqrt(B[i][i+1] * B[i][i+1] + B[i+1][i+1] * B[i+1][i+1]);
				c = B[i][i+1] / h;
				s = B[i+1][i+1] / h;
				OPgivenspre(A,i,i+1,i-1,c,s,epsA,sup);
				OPgivenspre(B,i,i+1,i+1,c,s,epsB,sup);
				OPgivenspre(Q,i,i+1,0,c,s,eps,Dim);
			}
			if(fabs(A[i+1][i-1]) > epsA)
            {
				h = sqrt(A[i+1][i-1] * A[i+1][i-1] + A[i+1][i] * A[i+1][i]);
				c = A[i+1][i] / h;
				s = A[i+1][i-1] / h;
				OPgivenspost(A,i-1,i,inf,c,s,epsA,sup);
				OPgivenspost(B,i-1,i,inf,c,s,epsB,sup);
				OPgivenspost(Z,i-1,i,0,c,s,eps,Dim);
			}
		}
		if(fabs(A[sup-1][sup-2]) > epsA)
		{
			h = sqrt(A[sup-1][sup-2] * A[sup-1][sup-2] + A[sup-1][sup-1] * A[sup-1][sup-1]);
			c = A[sup-1][sup-1] / h;
			s = A[sup-1][sup-2] / h;
			OPgivenspost(A,sup-2,sup-1,inf,c,s,epsA,sup);
			OPgivenspost(B,sup-2,sup-1,inf,c,s,epsB,sup);
			OPgivenspost(Z,sup-2,sup-1,0,c,s,eps,Dim);
		}
	}
}

void QZdeflaccion2(matriz A, matriz B, indice k, escalar epsA, escalar epsB, indice inf, indice sup)
{
	indice i;
    escalar h,c,s;

	if(k != inf)
	{ 		
		for(i = k;i < (sup-1);i++)
		{
            if(fabs(B[i+1][i+1]) > epsB)
			{
				h = sqrt(B[i][i+1] * B[i][i+1] + B[i+1][i+1] * B[i+1][i+1]);
				c = B[i][i+1] / h;
				s = B[i+1][i+1] / h;
				OPgivenspre(A,i,i+1,i-1,c,s,epsA,sup);
				OPgivenspre(B,i,i+1,i+1,c,s,epsB,sup);
			}
            if(fabs(A[i+1][i-1]) > epsA)
            {
				h = sqrt(A[i+1][i-1] * A[i+1][i-1] + A[i+1][i] * A[i+1][i]);
				c = A[i+1][i] / h;
				s = A[i+1][i-1] / h;
				OPgivenspost(A,i-1,i,inf,c,s,epsA,sup);
				OPgivenspost(B,i-1,i,inf,c,s,epsB,sup);
			}
		}
        if(fabs(A[sup-1][sup-2]) > epsA)
		{
			h = sqrt(A[sup-1][sup-2] * A[sup-1][sup-2] + A[sup-1][sup-1] * A[sup-1][sup-1]);
			c = A[sup-1][sup-1] / h;
			s = A[sup-1][sup-2] / h;
			OPgivenspost(A,sup-2,sup-1,inf,c,s,epsA,sup);
			OPgivenspost(B,sup-2,sup-1,inf,c,s,epsB,sup);
		}
	}
}

void QZexplic(matriz A, matriz B, matriz Q, matriz Z, escalar eps,
			escalar epsA, escalar epsB, indice inf, indice sup, unsigned int Dim)
{
	indice d,k,i;

	d = sup - inf;
	while(d > 2)
	{
		for(i = inf;i < sup;i++)
		{
			if(fabs(B[i][i]) < epsB)
			{
				QZdeflaccion(A,B,Q,Z,i,eps,epsA,epsB,inf,sup,Dim);
				sup--;
			}
		}
		QZstep(A,B,Q,Z,eps,epsA,epsB,inf,sup,Dim);
		k = sup - 1;
		while(k > inf)
		{
			if(fabs(A[k][k-1]) <= eps*(fabs(A[k][k]) + fabs(A[k-1][k-1])))
			{
				if((sup - k) > 2)
				{
						QZexplic(A,B,Q,Z,eps,epsA,epsB,k,sup,Dim);
						d = k - inf;
				}
				else
					d -= (sup - k);
				sup = k;
			}
			k--;
		}
	}
}

void QZexplic2(matriz A, matriz B, escalar eps, escalar epsA, escalar epsB, indice inf, indice sup)
{
	indice d,k,i;

	d = sup - inf;
	while(d > 2)
	{
		for(i = inf;i < sup;i++)
		{
			if(fabs(B[i][i]) < epsB)
			{
				QZdeflaccion2(A,B,i,epsA,epsB,inf,sup);
				sup--;
			}
		}
		QZstep2(A,B,eps,epsA,epsB,inf,sup);
		k = sup - 1;
		while(k > inf)
		{
			if(fabs(A[k][k-1]) <= eps*(fabs(A[k][k]) + fabs(A[k-1][k-1])))
			{
				if((sup - k) > 2)
				{
						QZexplic2(A,B,eps,epsA,epsB,k,sup);
						d = k - inf;
				}
				else
					d -= (sup - k);
				sup = k;
			}
			k--;
		}
	}
}

