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

void VAvap1x1Re(vector Lr, escalar A00, escalar B00, indice k)
{
	if(B00 == 0.0)
		Lr[k] = INFINITO;
	else
		Lr[k] = A00 / B00;
}

void VAvap2x2Re(vector Lr, escalar A00, escalar A01, escalar A10, escalar A11,
		escalar B00, escalar B01, escalar B11, indice k)
{
	escalar delta,p0,p1;

	if((B00 == 0.0) || (B11 == 0.0))
	{
		delta = A00 * A11 - A01 * A10;
		if(B00 == 0.0)
		{
			p0 = A00 * B11 - A10 * B01;
			Lr[k] = INFINITO;
			Lr[k+1] = delta / p0;
		}
		else
		{
			p1 = A11 * B00 - A10 * B01;
			Lr[k] = delta / p1;
			Lr[k+1] = INFINITO;
		}
	}
	else
    {
		p1 = (A10/B00)*(B01/B11) - A00/B00 - A11/B11;
		p0 = (A00*A11 - A01*A10) / (B00*B11);
		delta = p1*p1 - 4*p0;
		Lr[k] = .5*(-p1 + sqrt(delta));
		Lr[k+1] = .5*(-p1 - sqrt(delta));
	}
}

void VAvap2x2Im(vector Lr, vector Li, escalar A00, escalar A01, escalar A10, escalar A11,
		escalar B00, escalar B01, escalar B11, indice k)
{
	escalar delta,p0,p1;

	if((B00 == 0.0) || (B11 == 0.0))
	{
		delta = A00 * A11 - A01 * A10;
		if(B00 == 0.0)
		{
			p0 = A00 * B11 - A10 * B01;
			Lr[k] = INFINITO;
			Lr[k+1] = delta / p0;
			Li[k] = INFINITO;
            Li[k+1] = 0.0;
		}
		else
		{
			p1 = A11 * B00 - A10 * B01;
			Lr[k] = delta / p1;
			Lr[k+1] = INFINITO;
			Li[k] = 0.0;
            Li[k+1] = INFINITO;
		}
	}
	else
    {
		p1 = (A10/B00)*(B01/B11) - A00/B00 - A11/B11;
		p0 = (A00*A11 - A01*A10) / (B00*B11);
		delta = p1*p1 - 4*p0;
		if(delta > 0.0)
		{
			Lr[k] = .5*(-p1 + sqrt(delta));
			Lr[k+1] = .5*(-p1 - sqrt(delta));
			Li[k] = Li[k+1] = 0.0;
		}
		else
		{
			Lr[k] = Lr[k+1] = -.5*p1;
			Li[k] = .5*sqrt(-delta);
			Li[k+1] = -Li[k];
		}
	}
}

void VAvapRe(vector Lr, matriz A, matriz B, escalar eps, unsigned int Dim)
{
	indice k,l;

	for(k = 1;k < Dim;k++)
		if(fabs(A[k][k-1]) <= eps * (fabs(A[k][k]) + fabs(A[k-1][k-1])))
			A[k][k-1] = 0.0;

	k = Dim;
	l = k - 1;
	while(l > 0)
	{
		if(A[l][l-1] == 0.0)
		{
			if((k-l) == 1)
			{
				VAvap1x1Re(Lr,A[l][l],B[l][l],l);
			}
			else
			{
				VAvap2x2Re(Lr,A[l][l],A[l][l+1],A[l+1][l],A[l+1][l+1],B[l][l],B[l][l+1],B[l+1][l+1],l);
			}
			k = l;
		}
		l--;
	}
	if(k == 1)
		VAvap1x1Re(Lr,A[0][0],B[0][0],0);
	else
		VAvap2x2Re(Lr,A[0][0],A[0][1],A[1][0],A[1][1],B[0][0],B[0][1],B[1][1],0);
}

void VAvapIm(vector Lr, vector Li, matriz A, matriz B, escalar eps, unsigned int Dim)
{
	indice k,l;

	for(k = 1;k < Dim;k++)
		if(fabs(A[k][k-1]) <= eps * (fabs(A[k][k]) + fabs(A[k-1][k-1])))
			A[k][k-1] = 0.0;

	k = Dim;
	l = k - 1;
	while(l > 0)
	{
		if(A[l][l-1] == 0.0)
		{
			if((k-l) == 1)
			{
				VAvap1x1Re(Lr,A[l][l],B[l][l],l);
                Li[l] = 0.0;
			}
			else
			{
				VAvap2x2Im(Lr,Li,A[l][l],A[l][l+1],A[l+1][l],A[l+1][l+1],B[l][l],B[l][l+1],B[l+1][l+1],l);
			}
			k = l;
		}
		l--;
	}
	if(k == 1)
    {
		VAvap1x1Re(Lr,A[0][0],B[0][0],0);
		Li[0] = 0.0;
	}
	else
		VAvap2x2Im(Lr,Li,A[0][0],A[0][1],A[1][0],A[1][1],B[0][0],B[0][1],B[1][1],0);
}

