/**
    Autor : Antonio Mesa
    Version : v 1.0
        1-XI-1994
        10-II-1995
**/

#include <math.h>
#include "nuevo.h"
#include "iomatriz.h"
#include "fcmatriz.h"

escalar OPeps(div)
escalar div;
{
    escalar e,e2,x;
    e = 1.0;
    do
    {
        e2 = e;
        e /= div;
        x = 1.0 + e;
    }while(x > 1.0);
    return e2;
}

void OPsumamat(matriz C, matriz A, matriz B, unsigned int Dim)
{
	indice i,j;

	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
			C[i][j] = A[i][j] + B[i][j];
}

void OPsumasubmat(C,A,B,inf,sup)
matriz C,A,B;
indice inf,sup;
{
    indice i,j;

    for(i = inf;i < sup;i++)
        for(j = inf;j < sup;j++)
            C[i][j] = A[i][j] + B[i][j];
}

void OPrestamat(matriz C, matriz A, matriz B, unsigned int Dim)
{
	indice i,j;

	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
			C[i][j] = A[i][j] - B[i][j];
}

void OPrestasubmat(C,A,B,inf,sup)
matriz C,A,B;
indice inf,sup;
{
    indice i,j;

    for(i = inf;i < sup;i++)
        for(j = inf;j < sup;j++)
            C[i][j] = A[i][j] - B[i][j];
}

void OPprodmat(matriz C, matriz A, matriz B, unsigned int Dim)
{
	indice i,j,k;
	escalar ss;
	matriz T;

    switch(sel)
    {
    case PROD :
	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
		{
			ss = 0.0;
			for(k = 0;k < Dim;k++)
				ss += A[i][k] * B[k][j];
                C[i][j] = ss;
		}
        break;
    case PRETRANS :
	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
            {
                ss = 0.0;
                for(k = 0;k < Dim;k++)
                    ss += A[k][i] * B[k][j];
                C[i][j] = ss;
            }
        break;
    case POSTTRANS :
        for(i = 0;i < Dim;i++)
            for(j = 0;j < Dim;j++)
            {
                ss = 0.0;
                for(k = 0;k < Dim;k++)
                    ss += A[i][k] * B[j][k];
                C[i][j] = ss;
            }
        break;
    default : break;
    }
}

void OPprodmatmxn(C,A,filA,colA,B,filB,colB,sel)
matriz C,A;
unsigned int filA,colA;
matriz B;
unsigned int filB,colB;
int sel;
{
    indice i,j,k;
    escalar ss;

    switch(sel)
    {
    case PROD :
        for(i = 0;i < filA;i++)
            for(j = 0;j < colB;j++)
            {
                ss = 0.0;
                for(k = 0;k < filB;k++)
                    ss += A[i][k] * B[k][j];
                C[i][j] = ss;
            }
        break;
    case PRETRANS :
        for(i = 0;i < colA;i++)
            for(j = 0;j < colB;j++)
            {
                ss = 0.0;
                for(k = 0;k < filA;k++)
                    ss += A[k][i] * B[k][j];
                C[i][j] = ss;
            }
        break;
    case POSTTRANS :
        for(i = 0;i < filA;i++)
            for(j = 0;j < filB;j++)
            {
                ss = 0.0;
                for(k = 0;k < colA;k++)
                    ss += A[i][k] * B[j][k];
                C[i][j] = ss;
            }
        break;
    default : break;
    }
}

void OPprodsubmat(C,A,B,inf,sup)
matriz C,A,B;
indice inf,sup;
{
    indice i,j,k;
    escalar ss;

       for(i = inf;i < sup;i++)
        for(j = inf;j < sup;j++)
        {
            ss = 0.0;
            for(k = inf;k < sup;k++)
                ss += A[i][k] * B[k][j];
            C[i][j] = ss;
        }
}

void OPcpymat(matriz A, matriz B, unsigned int Dim)
{
	indice i,j;

    for(i = 0;i < Dim;i++)
    	for(j = 0;j < Dim;j++)
        	A[i][j] = B[i][j];
}

void OPsumavect(vector W, vector U, vector V, unsigned int Dim)
{
	indice i;

	for(i = 0;i < Dim;i++)
		W[i] = U[i] + V[i];
}

void OPrestavect(vector W, vector U, vector V, unsigned int Dim)
{
	indice i;

	for(i = 0;i < Dim;i++)
		W[i] = U[i] - V[i];
}

escalar OPprodint(vector U, vector V, unsigned int Dim)
{
	indice i;
	escalar pp;

	pp = 0.0;
	for(i = 0;i < Dim;i++)
			pp += U[i] * V[i];
	return pp;
}

escalar OPprodintmat(vector X, matriz A, vector Y, unsigned int Dim)
{
	indice i,j;
	escalar r;
	vector T;

	IOmemvect(&T,Dim);
	for(i = 0;i < Dim;i++)
	{
		r = 0.0;
		for(j = 0;j < Dim;j++)
			r += A[i][j] * Y[j];
		T[i] = r;
	}
	r = 0.0;
	for(i = 0;i < Dim;i++)
		r += X[i] * T[i];
	IOfreevect(T,Dim);
	return r;
}

void OPctemat(matriz C, matriz A, escalar K, unsigned int Dim)
{
	indice i,j;

	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
			C[i][j] = K * A[i][j];
}

void OPctesubmat(C,A,K,inf,sup)
matriz C,A;
escalar K;
indice inf,sup;
{
    indice i,j;

    for(i = inf;i < sup;i++)
        for(j = inf;j < sup;j++)
            C[i][j] = K * A[i][j];
}

void OPctevect(vector V, vector U, escalar K, unsigned int Dim)
{
	indice i;

	for(i = 0;i < Dim;i++)
		V[i] = K * U[i];
}

void OPcpyvect(vector V, vector U, unsigned int Dim)
{
	indice i;

    for(i = 0;i < Dim;i++)
    	V[i] = U[i];
}

void OPmatvect(vector V, matriz A, vector U, unsigned int Dim)
{
	indice i,j;
	escalar ss;
	vector W;

	IOmemvect(&W,Dim);
	for(i = 0;i < Dim;i++)
	{
		ss = 0.0;
		for(j = 0;j < Dim;j++)
			ss += A[i][j] * U[j];
		W[i] = ss;
	}
	for(i = 0;i < Dim;i++)
		V[i] = W[i];
	IOfreevect(W,Dim);
}

void OPmatdiag(matriz C, matriz A, vector D, unsigned int Dim)
{
	indice i,j;
	escalar dd;
	matriz T;

	IOmemmat(&T,Dim);
	for(j = 0;j < Dim;j++)
	{
		dd = D[j];
		for(i = 0;i < Dim;i++)
			T[i][j] = A[i][j] * dd;
	}
	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
			C[i][j] = T[i][j];
	IOfreemat(T,Dim);
}

void OPdiagmat(matriz C, matriz A, vector D, unsigned int Dim)
{
	indice i,j;
	escalar dd;
	matriz T;

	IOmemmat(&T,Dim);
	for(i = 0;i < Dim;i++)
	{
		dd = D[i];
		for(j = 0;j < Dim;j++)
			T[i][j] = dd * A[i][j];
	}
	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
			C[i][j] = T[i][j];
	IOfreemat(T,Dim);
}

void OPinvtsup(matriz Ui, matriz U, unsigned int Dim)
{
	indice i,j,k;
	escalar vv;

	for(i = 0;i < Dim;i++)
		Ui[i][i] = 1.0 / U[i][i];
	for(i = 1;i < Dim;i++)
		for(j = 0;j < (Dim-i);j++)
		{
			vv = U[j][i+j] * Ui[i+j][i+j];
			for(k = 1;k < i;k++)
				vv += U[j][i+j-k] * Ui[i+j-k][j+i];
			Ui[j][i+j] = -vv / U[j][j];
		}
	for(i = 1;i < Dim;i++)
		for(j = 0;j < i;j++)
			Ui[i][j] = 0.0;
}

void OPinvtinf(matriz Li, matriz L, unsigned int Dim)
{
	indice i,j,k;
	escalar mm;

	for(i = 0;i < Dim;i++)
		Li[i][i] = 1.0 / L[i][i];
	for(i = 1;i < Dim;i++)
		for(j = 0;j < (Dim-i);j++)
		{
			mm = L[i+j][j] * Li[j][j];
			for(k = 1;k < i;k++)
				mm += L[i+j][i+j-k] * Li[i+j-k][j];
			Li[i+j][j] = -mm / L[i+j][i+j];
		}
	for(i = 1;i < Dim;i++)
		for(j = i+1;j < Dim;j++)
			Li[i][j] = 0.0;
}

void OPinvdiag(vector Di, vector D, unsigned int Dim)
{
	indice i;

	for(i = 0;i < Dim;i++)
		Di[i] = 1 / D[i];
}

escalar OPnorma2(vector V, unsigned int Dim)
{
	indice i;
	escalar x;

	x = 0.0;
	for(i = 0;i < Dim;i++)
		x += V[i] * V[i];
	return sqrt(x);
}

escalar OPnormamax(vector V, unsigned int Dim)
{
	indice i;
    escalar x;

    x = 0.0;
    for(i = 0; i < Dim;i++)
    	if(fabs(V[i]) > x) x = fabs(V[i]);

    return x;
}

escalar OPnormamaxmat(matriz A, unsigned int Dim)
{
	indice i,j;
    escalar x,y;

    x = 0.0;
    for(i = 0;i < Dim;i++)
    {
    	y = 0.0;
    	for(j = 0;j < Dim;j++)
        {
        	y += fabs(A[i][j]);
        }
        if(y > x) x = y;
    }
    return x;
}

void OPtransmat(matriz A, unsigned int Dim)
{
	indice i,j;
	matriz T;

	IOmemmat(&T,Dim);
	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
			T[i][j] = A[j][i];
	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
			A[i][j] = T[i][j];
	IOfreemat(T,Dim);
}

void OPfiltramat(matriz A, escalar e, unsigned int Dim)
{
	indice i,j;

    for(i = 0;i < Dim;i++)
    	for(j = 0;j < Dim;j++)
        	if(fabs(A[i][j]) < e) A[i][j] = 0.0;
}

void OPsisteqsup(matriz U, vector X, vector Y, unsigned int Dim)
{
	indice i,j,k;
	escalar u;

	for(i = 0;i < Dim;i++)
	{
		k = Dim - 1 - i;
		u = 0.0;
		for(j = (Dim-i);j < Dim;j++)
			u += X[j] * U[k][j];
		X[k] = (Y[k] - u) / U[k][k];
	}
}

void OPsisteqinf(matriz L, vector X, vector Y, unsigned int Dim)
{
	indice i,j;
    escalar l;

    for(i = 0;i < Dim;i++)
    {
        l = 0.0;
        for(j = 0;j < i;j++)
        	l += X[j] * L[i][j];
        X[i] = (Y[i] - l) / L[i][i];
    }
}

void OPidentidad(matriz I, unsigned int Dim)
{
	indice i,j;

	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
		{
			if(i == j) I[i][i] = 1.0;
			else I[i][j] = 0.0;
		}
}

void OPceromat(matriz O, unsigned int Dim)
{
	indice i,j;
	for(i = 0;i < Dim;i++)
		for(j = 0;j < Dim;j++)
			O[i][j] = 0.0;
}

void OPunovect(vector V, unsigned int Dim)
{
	indice i;

    for(i = 0;i < Dim;i++)
    	V[i]  = 1.0;
}

void OPcerovect(vector V, unsigned int Dim)
{
	indice i;

    for(i = 0;i < Dim;i++)
    	V[i]  = 0.0;
}

void OPsort(vector X, unsigned int Dim)
{
	indice i,j,jmin;
    escalar xmin,t;

    for(i = 0;i < Dim;i++)
    {
       	xmin = X[i];
        jmin = i;
    	for(j = i;j < Dim;j++)
        {
        	if(fabs(X[j]) < xmin)
            {
            	jmin = j;
            	xmin = fabs(X[j]);
            }
        }
		t = X[jmin];
		X[jmin] = X[i];
		X[i] = t;
    }
}

void OPgivenspre(matriz A, indice i, indice j, escalar c, escalar s, unsigned int Dim)
{
	indice k;
    escalar p,q;

    for(k = colini;k < Dim;k++)
	{
		p = A[i][k];
		q = A[j][k];
        if(fabs(A[i][k] = p * c + q * s) <= epsA) A[i][k] = 0.0;
        if(fabs(A[j][k] = -p * s + q * c) <= epsA) A[j][k] = 0.0;
	}
}

void OPgivenspost(matriz A, indice i, indice j, escalar c, escalar s, unsigned int Dim)
{
	indice k;
    escalar p,q;

    for(k = filini;k < Dim;k++)
	{
		p = A[k][i];
		q = A[k][j];
        if(fabs(A[k][i] = p * c - q * s) <= epsA) A[k][i] = 0.0;
        if(fabs(A[k][j] = p * s + q * c) <= epsA) A[k][j] = 0.0;
	}
}

void OPelempre1(matriz A, indice i, indice j, escalar x, unsigned int Dim)
{
	indice k;
    escalar p,q;

    for(k = colini;k < Dim;k++)
    {
    	p = A[i][k];
        q = A[j][k];
        if(fabs(A[j][k] = x * p + q) <= epsA) A[j][k] = 0.0;
    }
}

void OPelempre2(matriz A, indice i, indice j, escalar x, unsigned int Dim)
{
	indice k;
    escalar p,q;

    for(k = colini;k < Dim;k++)
    {
    	p = A[i][k];
        q = A[j][k];
        A[i][k] = q;
        if(fabs(A[j][k] = x * q + p) <= epsA) A[j][k] = 0.0;
    }
}

void OPelempost1(matriz A, indice i, indice j, escalar x, unsigned int Dim)
{
	indice k;
    escalar p,q;

    for(k = filini;k < Dim;k++)
    {
    	p = A[k][i];
        q = A[k][j];
        if(fabs(A[k][i] = p + x * q) <= epsA) A[k][i] = 0.0;
    }
}

void OPelempost2(matriz A, indice i, indice j, escalar x, unsigned int Dim)
{
	indice k;
    escalar p,q;

    for(k = filini;k < Dim;k++)
    {
    	p = A[k][i];
        q = A[k][j];
        if(fabs(A[k][i] = q + x * p) <= epsA) A[k][i] = 0.0;
        A[k][j] = p;
    }
}

void OPcambiafil(matriz A, indice i, indice j, unsigned int Dim)
{
	indice k;
    escalar t;

    for(k = 0;k < Dim;k++)
    {
    	t = A[i][k];
        A[i][k] = A[j][k];
        A[j][k] = t;
    }
}

void OPcambiacol(matriz A, indice i, indice j, unsigned int Dim)
{
	indice k;
    escalar t;

    for(k = 0;k < Dim;k++)
    {
    	t = A[k][i];
        A[k][i] = A[k][j];
        A[k][j] = t;
    }
}

void OPsubdiaghess(matriz A, matriz B, unsigned int Dim)
{
	indice i,j;
    escalar x,y;

	for(i = 1;i < Dim;i++)
	{
        x = y / A[i][i-1];
        A[i][i-1] = y;
		for(j = i;j < Dim;j++)
		{
			A[i][j] = x * A[i][j];
			B[i][j] = x * B[i][j];
		}
	}
}
        	
        	