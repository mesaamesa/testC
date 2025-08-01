/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa												*/
/*																		*/
/*	Versi�n : 															*/
/*		   		v 1.0		1-XI-1994			10-II-1995				*/
/*																		*/
/************************************************************************/

#include <math.h>

#include "defines.h"
#include "io.h"

void BImirasturm(unsigned int *pS, vector P, unsigned int Dim)
{
	unsigned int s,r,pr,ps;

	s = 0;
	pr = ps = POSITIVO;
	for(r = 0;r < Dim;r++)
	{
		if(fabs(P[r]) < EPS)
			pr = CERO;
		else
			if(P[r] > 0.0)
				pr = POSITIVO;
			else
				pr = NEGATIVO;
		if(pr == CERO)
			if(ps == POSITIVO)
				pr = NEGATIVO;
			else
				pr = POSITIVO;
		else
			if(pr == ps)
				s++;
		ps = pr;
	}
	*pS = s;
}																				 

void BIsturm(vector P, matriz A, matriz B, escalar mu, unsigned int Dim)
{
    unsigned int r,i;
    int signo;
    escalar p;
	vector x;

	IOmemvect(&x,Dim);
	signo = +1;
	P[0] = A[Dim-1][Dim-1] - mu * B[Dim-1][Dim-1];
	x[Dim-2] = -P[0];
	for(r = 2;r <= Dim;r++)
	{
		p = A[Dim-r][Dim-1] - mu * B[Dim-r][Dim-1];
		for(i = 2;i <= r;i++)
		{
			p += (A[Dim-r][Dim-i] - mu * B[Dim-r][Dim-i]) * x[Dim-i];
		}
		signo = - signo;
		if(signo > 0)
			P[r-1] = p;
		else
			P[r-1] = -p;
        if(r != Dim)
			x[Dim-r-1] = -p;
	}
    IOfreevect(x,Dim);
}										 

void BIbuscainterval(unsigned int *pK, escalar *pa,escalar *pb, matriz A, matriz B,
 escalar mu, unsigned int Dim)
{
    bool fin;
    unsigned int k,kpos,kneg;
	escalar incmu,mupos,muneg,a0,b0;
	vector P;

	IOmemvect(&P,Dim);

    fin = FALSE; 
	incmu = 2 * mu * EPS;
	BIsturm(P,A,B,mu,Dim);
	BImirasturm(&k,P,Dim);
/*trazaprintf("\nmu = %1.17g k = %d",mu,k);IOputvect(P,"P",'s',Dim);getchar();*/
	while(!fin)
	{
/*trazagetchar();printf("\nincmu = %g",incmu);*/
		muneg = mu - incmu;
		mupos = mu + incmu;
		BIsturm(P,A,B,mupos,Dim);
		BImirasturm(&kpos,P,Dim);
/*trazaprintf("\nmupos = %1.17g kpos = %d",mupos,kpos);*/
		BIsturm(P,A,B,muneg,Dim);
		BImirasturm(&kneg,P,Dim);
/*trazaprintf("\nmuneg = %1.17g kneg = %d",muneg,kneg);*/
		if((k == kpos) && (k == kneg))
		{
			incmu *= 2;
			if(incmu > COTASUP)
				fin = TRUE;
		}
		else
		{
			if(kneg != k)
			{
				k = kneg;
				a0 = muneg;
				b0 = mu;
			}
			else
			{
				a0 = mu;
				b0 = mupos;
			}
			fin = TRUE;
		}
	}
	*pK = k;
	*pa = a0;
    *pb = b0;

	IOfreevect(P,Dim);
}												

void BIbisec(escalar *plambda, escalar a0, escalar b0, unsigned int k, matriz A, matriz B, unsigned int Dim)
{
    unsigned int sk;
	escalar a,b,c;
	vector P;

    IOmemvect(&P,Dim);

	a = a0;
	b = b0;
	while((b - a) > EPS*c)
	{
		c = (b + a) / 2;
		BIsturm(P,A,B,c,Dim);
		BImirasturm(&sk,P,Dim);
		if(sk < k)	 
			b = c;
		else
			a = c;
	}
	*plambda = c;

	IOfreevect(P,Dim);
}							    	


			