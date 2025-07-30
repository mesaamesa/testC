/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa												*/
/*																		*/
/*	Versi¢n : 															*/
/*		   		v 1.0		1-XI-1994		10-II-1995					*/
/*																		*/
/************************************************************************/

#include "nuevo.h"
#include "iomatriz.h"
#include "opmatriz.h"
#include <math.h>
/*traza*/
#include <stdio.h>

void FCldu(L,U,D,A,Dim)
matriz L,U,A;
vector D;
unsigned int Dim;
{
	indice i,j,k;
	escalar ll,uu,dd;

	D[0] = A[0][0];
	U[0][0] = 1.0;
	L[0][0] = 1.0;
	for(j = 1;j < Dim;j++)
		U[0][j] = A[0][j] / D[0];
	for(i = 1;i < Dim;i++)
		L[i][0] = A[i][0] / D[0];
	for(i = 1;i < Dim;i++)
		for(j = 1;j < Dim;j++)
		{
			if(i == j)
			{
				dd = A[i][i];
				for(k = 0;k < i;k++)
					dd -= D[k] * L[i][k] * U[k][j];
				D[i] = dd;
				U[i][i] = 1.0;
				L[i][i] = 1.0;
			}
			if(j > i)
			{
				uu = A[i][j];
				for(k = 0;k < i;k++)
					uu -= D[k] * L[i][k] * U[k][j];
				U[i][j] = uu / D[i];
			}
			if(i > j)
			{
				ll = A[i][j];
				for(k = 0;k < j;k++)
					ll -= D[k] * L[i][k] * U[k][j];
				L[i][j] = ll / D[j];
			}
		}
}

void FClu(L,U,A,Dim)
matriz L,U,A;
unsigned int Dim;
{
	indice i,j,k;
	escalar ll,uu;

	U[0][0] = A[0][0];
	L[0][0] = 1.0;
	for(j = 1;j < Dim;j++)
		U[0][j] = A[0][j];
	for(i = 1;i < Dim;i++)
		for(j = 0;j < Dim;j++)
		{
			if(j >= i)
			{
				uu = A[i][j];
				for(k = 0;k < i;k++)
					uu -= L[i][k] * U[k][j];
				U[i][j] = uu;
                L[i][i] = 1.0;
			}
			if(i > j)
			{
				ll = A[i][j];
				for(k = 0;k < j;k++)
					ll -= L[i][k] * U[k][j];
				L[i][j] = ll / U[j][j];
			}
		}
}

void FCchol(U,A,Dim)
matriz U,A;
unsigned int Dim;
{
	indice i,j,k;
	escalar uu;

	for(i = 0;i < Dim;i++)
		for(j = i;j < Dim;j++)
			if(i == j)
			{
				uu = 0.0;
				for(k = 0;k < i;k++)
					uu += U[k][i] * U[k][i];
				U[i][i] = sqrt(A[i][i] - uu);
			}
			else
			{
				uu = 0.0;
				for(k = 0;k < i;k++)
					uu += U[k][i] * U[k][j];
				U[i][j] = (A[i][j] - uu) / U[i][i];
			}
}
    
void FCqrgram(Q,R,A,Dim)
matriz Q,R,A;
unsigned int Dim;
{
	indice i,j,k,l;
	escalar rkk,rkl;
	vector v,al,q,qtemp;
	matriz AA;

	IOmemmat(&AA,Dim);
	IOmemvect(&v,Dim);
	IOmemvect(&al,Dim);
	IOmemvect(&q,Dim);
	IOmemvect(&qtemp,Dim);

	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
			AA[i][j] = A[i][j];
	for(k = 0;k < Dim;k++)
	{
		for(i = 0;i < Dim;i++)
			v[i] = AA[i][k];
		R[k][k] = rkk = OPnorma2(v,Dim);
		for(i = 0;i < Dim;i++)
			Q[i][k]  = q[i] = v[i] / rkk;
		for(l = k+1;l < Dim;l++)
		{
			for(i = 0;i < Dim;i++)
				al[i] = AA[i][l];
			R[k][l] = rkl = OPprodint(q,al,Dim);
			OPctevect(qtemp,q,rkl,Dim);
			OPrestavect(al,al,qtemp,Dim);
			for(i = 0;i < Dim;i++)
				AA[i][l] = al[i];
		}
	}
	IOfreemat(AA,Dim);
	IOfreevect(v,Dim);
	IOfreevect(al,Dim);
	IOfreevect(q,Dim);
	IOfreevect(qtemp,Dim);
}

void FCqrhouse(Q,R,A,Dim)
matriz Q,R,A;
unsigned int Dim;
{
	indice i,j,k,l;
	escalar sk,alfa,x;
	vector v;
	matriz P,T;

	IOmemmat(&P,Dim);
	IOmemmat(&T,Dim);
	IOmemvect(&v,Dim);

	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
		{
			R[i][j] = A[i][j];
			if(i == j)
				Q[i][i] = 1.0;
			else
				Q[i][j] = 0.0;
		}
	for(k = 0;k < Dim;k++)
	{
		for(i = 0;i < k;i++)
			v[i] = 0.0;
		for(i = k;i < Dim;i++)
			v[i] = R[i][k];
		sk = OPnorma2(v,Dim);
		if(R[k][k] > 0.0)
			sk = -sk;
		v[k] = v[k] - sk;
		alfa = 1 / (sk * (sk - R[k][k]));
		for(i = k;i < Dim;i++)
			for(j = k;j < Dim;j++)
				if(i == j)
					P[i][j] = 1 - alfa * v[i] * v[j];
				else
					P[i][j] = -alfa * v[i] * v[j];
		for(i = k;i < Dim;i++)
			for(j = k;j < Dim;j++)
			{
				x = 0.0;
				for(l = k;l < Dim;l++)
					x += P[i][l] * R[l][j];
				T[i][j] = x;
			}
		for(i = k;i < Dim;i++)
			for(j = k;j < Dim;j++)
				R[i][j] = T[i][j];
		for(i = 0;i < Dim;i++)
			for(j = k;j < Dim;j++)
			{
				x = 0.0;
				for(l = k;l < Dim;l++)
					x += Q[i][l] * P[l][j];
				T[i][j] = x;
			}
		for(i = 0;i < Dim;i++)
			for(j = k;j < Dim;j++)
				Q[i][j] = T[i][j];
	}
	IOfreemat(P,Dim);
	IOfreemat(T,Dim);
	IOfreevect(v,Dim);
}

void FCqrgivens(A,Q,epsA,eps,Dim)
matriz A,Q;
escalar epsA,eps;
unsigned int Dim;
{
	indice i,j;
	escalar s,c,h;

	OPidentidad(Q,Dim);
	for(j = 0;j < (Dim-1);j++)
	{
		for(i = (Dim-1);i > j;i--)
		{
			if(fabs(A[i][j]) > epsA)
			{
				h = sqrt(A[i][j] * A[i][j] + A[i-1][j] * A[i-1][j]);
				c = A[i-1][j] / h;
				s = A[i][j] / h;
				OPgivenspre(A,i-1,i,j,c,s,epsA,Dim);
				OPgivenspre(Q,i-1,i,0,c,s,eps,Dim);
			}
		}
	}
}						               	

void FClrelem(L,A,epsA,Dim)
matriz L,A;
escalar epsA;
unsigned int Dim;
{
	indice i,j;
	escalar x;

    OPidentidad(L,Dim);
    for(j = 0;j < (Dim-1);j++)
    {
    	for(i = (Dim-1);i > j;i--)
		{
			if((fabs(A[i-1][j]) >= epsA) || (fabs(A[i][j]) >= epsA))
			{
				if(fabs(A[i-1][j]) >= fabs(A[i][j]))
	            {
					x = - A[i][j] / A[i-1][j];
					OPelempre1(A,i-1,i,j,x,epsA,Dim);
            	    OPelempre1(L,i-1,i,j,x,epsA,Dim);
	            }
    	        else
        	    {
					x = - A[i-1][j] / A[i][j];
	                OPelempre2(A,i-1,i,j,x,epsA,Dim);
    	            OPelempre2(L,i-1,i,j,x,epsA,Dim);
        	    }
			}
		}
    }
}

void FChesstriorto(A,B,Q,Z,eps,epsA,epsB,Dim)
matriz A,B,Q,Z;
escalar eps,epsA,epsB;
unsigned int Dim;
{
	indice i,j;
	escalar c,s,h;

	for(j = 0;j < (Dim-2);j++)
		for(i = (Dim-1);i >= (j+2);i--)
		{
			if(fabs(A[i][j]) >= epsA)
			{
				h = sqrt(A[i-1][j] * A[i-1][j] + A[i][j] * A[i][j]);
				c = A[i-1][j] / h;
				s = A[i][j] / h;
                OPgivenspre(A,i-1,i,j,c,s,epsA,Dim);
                OPgivenspre(B,i-1,i,j,c,s,epsB,Dim);
                OPgivenspre(Q,i-1,i,0,c,s,eps,Dim);
			}
			if(fabs(B[i][i-1]) >= epsB)
			{
				h = sqrt(B[i][i] * B[i][i] + B[i][i-1] * B[i][i-1]);
				c = B[i][i] / h;
				s = B[i][i-1] / h;
				OPgivenspost(A,i-1,i,0,c,s,epsA,Dim);
                OPgivenspost(B,i-1,i,0,c,s,epsB,Dim);
                OPgivenspost(Z,i-1,i,0,c,s,eps,Dim);
			}
		}
}

void FChesstriorto2(A,B,epsA,epsB,Dim)
matriz A,B;
escalar epsA,epsB;
unsigned int Dim;
{
	indice i,j;
	escalar c,s,h;

	#ifdef _DEBUG_
	printf("\nfile:%s (%d)",__FILE__,__LINE__);
	#endif
	for(j = 0;j < (Dim-2);j++)
	{
		#ifdef _DEBUG_
		printf("\nfile:%s (%d)",__FILE__,__LINE__);
		#endif
		for(i = (Dim-1);i >= (j+2);i--)
		{
			#ifdef _DEBUG_
			printf("\nfile:%s (%d)",__FILE__,__LINE__);
			#endif
			if(fabs(A[i][j]) >= epsA)
			{
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				h = sqrt(A[i-1][j] * A[i-1][j] + A[i][j] * A[i][j]);
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				c = A[i-1][j] / h;
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				s = A[i][j] / h;
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				OPgivenspre(A,i-1,i,j,c,s,epsA,Dim);
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				OPgivenspre(B,i-1,i,j,c,s,epsB,Dim);
			}
			if(fabs(B[i][i-1]) >= epsB)
			{
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				h = sqrt(B[i][i] * B[i][i] + B[i][i-1] * B[i][i-1]);
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				c = B[i][i] / h;
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				s = B[i][i-1] / h;
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				OPgivenspost(A,i-1,i,0,c,s,epsA,Dim);
				#ifdef _DEBUG_
				printf("\nfile:%s (%d)",__FILE__,__LINE__);
				#endif
				OPgivenspost(B,i-1,i,0,c,s,epsB,Dim);
			}
		}
	}
}

void FChesstrielem(A,B,L,M,eps,epsA,epsB,Dim)
matriz A,B,L,M;
escalar eps,epsA,epsB;
unsigned int Dim;
{
	indice i,j;
    escalar x;

    OPidentidad(L,Dim);
    OPidentidad(M,Dim);
    for(j = 0;j < (Dim-2);j++)
    {
    	for(i = (Dim-1);i > (j+1);i--)
        {
        	if((fabs(A[i][j]) >= epsA) || (fabs(A[i-1][j]) >= epsA))
			{
				if(fabs(A[i-1][j]) >= fabs(A[i][j]))
                {
                	x = - A[i][j] / A[i-1][j];
                    OPelempre1(A,i-1,i,j,x,epsA,Dim);
                    OPelempre1(B,i-1,i,j,x,epsB,Dim);
                    OPelempre1(L,i-1,i,0,x,eps,Dim);
				}
                else
                {
                	x = - A[i-1][j] / A[i][j];
                    OPelempre2(A,i-1,i,j,x,epsA,Dim);
                    OPelempre2(B,i-1,i,j,x,epsB,Dim);
                    OPelempre2(L,i-1,i,0,x,eps,Dim);
                }
                if(fabs(B[i][i]) >= fabs(B[i][i-1]))
                {
                	x = - B[i][i-1] / B[i][i];
                    OPelempost1(A,i-1,i,0,x,epsA,Dim);
                    OPelempost1(B,i-1,i,0,x,epsB,Dim);
                    OPelempost1(M,i-1,i,0,x,eps,Dim);
                }
                else
                {
                	x = - B[i][i] / B[i][i-1];
                    OPelempost2(A,i-1,i,0,x,epsA,Dim);
                    OPelempost2(B,i-1,i,0,x,epsB,Dim);
                    OPelempost2(M,i-1,i,0,x,eps,Dim);
                }
            }
        }
    }
}

void FChesstrielem2(A,B,epsA,epsB,Dim)
matriz A,B;
escalar epsA,epsB;
unsigned int Dim;
{
	indice i,j;
    escalar x;

    for(j = 0;j < (Dim-2);j++)
    {
    	for(i = (Dim-1);i > (j+1);i--)
        {
        	if((fabs(A[i][j]) >= epsA) || (fabs(A[i-1][j]) >= epsA))
			{
				if(fabs(A[i-1][j]) >= fabs(A[i][j]))
                {
                	x = - A[i][j] / A[i-1][j];
                    OPelempre1(A,i-1,i,j,x,epsA,Dim);
                    OPelempre1(B,i-1,i,j,x,epsB,Dim);
				}
                else
                {
                	x = - A[i-1][j] / A[i][j];
                    OPelempre2(A,i-1,i,j,x,epsA,Dim);
                    OPelempre2(B,i-1,i,j,x,epsB,Dim);
                }
                if(fabs(B[i][i]) >= fabs(B[i][i-1]))
                {
                	x = - B[i][i-1] / B[i][i];
                    OPelempost1(A,i-1,i,0,x,epsA,Dim);
                    OPelempost1(B,i-1,i,0,x,epsB,Dim);
                }
                else
                {
                	x = - B[i][i] / B[i][i-1];
                    OPelempost2(A,i-1,i,0,x,epsA,Dim);
                    OPelempost2(B,i-1,i,0,x,epsB,Dim);
                }
            }
        }
    }
}

