/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa												*/
/*																		*/
/*	Versi¢n : 															*/
/*		   		v 1.0		1-XI-1994			...						*/
/*																		*/
/************************************************************************/

#include "nuevo.h"
#include "iomatriz.h"
#include "opmatriz.h"
#include "fcmatriz.h"
#include <math.h>
#include <stdio.h>
#include <conio.h>

void house3(escalar *M3, escalar *v3)
{
	escalar y,w;
    y = sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
    w = y*(v3[0] - y);
    M3[0] = v3[0] / y;
    M3[1] = v3[1] / y;
    M3[2] = v3[2] / y;
    M3[3] = 1 + v3[1]*v3[1] / w;
    M3[4] = v3[1]*v3[2] / w;
    M3[5] = 1 + v3[2]*v3[2] / w;
}

void house2(escalar *M2, escalar *v2)
{
	escalar y;
    y = sqrt(v2[0]*v2[0] + v2[1]*v2[1]);
    M2[0] = v2[0] / y;
    M2[1] = v2[1] / y;
    M2[2] = -v2[0] / y;
}

void EIqzimplicito(matriz A, matriz B, unsigned int Dim)
{
	indice i,k,m,n;
    escalar r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,t0,t1,t2;
    escalar v3[3],v2[2],M3[6],M2[3];

    m = Dim-2;
    n = Dim-1;
    r0 = A[m][m] / B[m][m];
    r1 = A[0][0] / B[0][0];
    r2 = A[n][n] / B[n][n];
    r3 = A[m][n] / B[n][n];
    r4 = A[n][m] / B[m][m];
    r5 = B[m][n] / B[n][n];
    r6 = B[0][0] / A[1][0];
    r7 = A[0][1] / B[1][1];
    r8 = B[0][1] / B[1][1];
    r9 = A[1][1] / B[1][1];
    r10 = A[1][0] / B[0][0];
    r11 = A[2][1] / B[1][1];
	v3[0] = ((r0 - r1)*(r2 - r1) - r3*r4 + r4*r5*r1)*r6 + r7 - r1*r8;
    v3[1] = (r9 - r1) - r10*r8 - (r0 - r1) - (r2 - r1) + r4*r5;
    v3[2] = r11;

    for(k = 0;k < (Dim-2);k++)
    {
    	house3(M3,v3);
        for(i = 0;i < Dim;i++)
        {
        	t0 = A[k][i];
            t1 = A[k+1][i];
            t2 = A[k+2][i];
            A[k][i] = M3[0]*t0 + M3[1]*t1 + M3[2]*t2;
            A[k+1][i] = M3[1]*t0 + M3[3]*t1 + M3[4]*t2;
            A[k+2][i] = M3[2]*t0 + M3[4]*t1 + M3[5]*t2;
        }
        for(i = 0;i < Dim;i++)
        {
        	t0 = B[k][i];
            t1 = B[k+1][i];
            t2 = B[k+2][i];
            B[k][i] = M3[0]*t0 + M3[1]*t1 + M3[2]*t2;
            B[k+1][i] = M3[1]*t0 + M3[3]*t1 + M3[4]*t2;
            B[k+2][i] = M3[2]*t0 + M3[4]*t1 + M3[5]*t2;
        }
/*        for(i = 0;i < Dim;i++)
        {
        	t0 = Q[k][i];
            t1 = Q[k+1][i];
            t2 = Q[k+2][i];
            Q[k][i] = M3[0]*t0 + M3[1]*t1 + M3[2]*t2;
            Q[k+1][i] = M3[1]*t0 + M3[3]*t1 + M3[4]*t2;
            Q[k+2][i] = M3[2]*t0 + M3[4]*t1 + M3[5]*t2;
        }
*/
        v3[0] = B[k+2][k];
        v3[1] = B[k+2][k+1];
        v3[2] = B[k+2][k+2];
        house3(M3,v3);
        for(i = 0;i < Dim;i++)
        {
        	t0 = A[i][k];
            t1 = A[i][k+1];
            t2 = A[i][k+2];
            A[i][k] = M3[2]*t0 + M3[4]*t1 + M3[5]*t2;
            A[i][k+1] = M3[1]*t0 + M3[3]*t1 + M3[4]*t2;
            A[i][k+2] = M3[0]*t0 + M3[1]*t1 + M3[2]*t2;
        }
        for(i = 0;i < Dim;i++)
        {
        	t0 = B[i][k];
            t1 = B[i][k+1];
            t2 = B[i][k+2];
            B[i][k] = M3[2]*t0 + M3[4]*t1 + M3[5]*t2;
            B[i][k+1] = M3[1]*t0 + M3[3]*t1 + M3[4]*t2;
            B[i][k+2] = M3[0]*t0 + M3[1]*t1 + M3[2]*t2;
        }
/*        for(i = 0;i < Dim;i++)
        {
        	t0 = Z[i][k];
            t1 = Z[i][k+1];
            t2 = Z[i][k+2];
            Z[i][k] = M3[2]*t0 + M3[4]*t1 + M3[5]*t2;
            Z[i][k+1] = M3[1]*t0 + M3[3]*t1 + M3[4]*t2;
            Z[i][k+2] = M3[0]*t0 + M3[1]*t1 + M3[2]*t2;
        }
*/
        v2[0] = B[k+1][k];
        v2[1] = B[k+1][k+1];
        house2(M2,v2);
        for(i = 0;i < Dim;i++)
        {
        	t0 = A[i][k];
            t1 = A[i][k+1];
            A[i][k] = M2[1]*t0 + M2[2]*t1;
            A[i][k+1] = M2[0]*t0 + M2[1]*t1;
		} 
        for(i = 0;i < Dim;i++)
        {
        	t0 = B[i][k];
            t1 = B[i][k+1];
            B[i][k] = M2[1]*t0 + M2[2]*t1;
            B[i][k+1] = M2[0]*t0 + M2[1]*t1;
		} 
/*        for(i = 0;i < Dim;i++)
        {
        	t0 = Z[i][k];
            t1 = Z[i][k+1];
            Z[i][k] = M2[1]*t0 + M2[2]*t1;
            Z[i][k+1] = M2[0]*t0 + M2[1]*t1;
		}
*/
        v3[0] = A[k+1][k];
        v3[1] = A[k+2][k];
        if(k < (Dim-3)) v3[2] = A[k+3][k];
	}
    v2[0] = v3[0];
    v2[1] = v3[1];
    house2(M2,v2);
    for(i = (Dim-3);i < Dim;i++)
    {
    	t0 = A[Dim-2][i];
        t1 = A[Dim-1][i];
        A[Dim-2][i] = M2[0]*t0 + M2[1]*t1;
        A[Dim-1][i] = M2[1]*t0 + M2[2]*t1;
    }
    for(i = (Dim-2);i < Dim;i++)
    {
    	t0 = B[Dim-2][i];
        t1 = B[Dim-1][i];
        B[Dim-2][i] = M2[0]*t0 + M2[1]*t1;
        B[Dim-1][i] = M2[1]*t0 + M2[2]*t1;
    }
/*    for(i = 0;i < Dim;i++)
    {
    	t0 = Q[Dim-2][i];
        t1 = Q[Dim-1][i];
        Q[Dim-2][i] = M2[0]*t0 + M2[1]*t1;
        Q[Dim-1][i] = M2[1]*t0 + M2[2]*t1;
    }
*/
    v2[0] = B[Dim-1][Dim-2];
    v2[1] = B[Dim-1][Dim-1];
    house2(M2,v2);
    for(i = 0;i < Dim;i++)
    {
    	t0 = A[i][Dim-2];
        t1 = A[i][Dim-1];
        A[i][Dim-2] = M2[1]*t0 + M2[2]*t1;
        A[i][Dim-1] = M2[0]*t0 + M2[1]*t1;
    }
    for(i = 0;i < Dim;i++)
    {
    	t0 = B[i][Dim-2];
        t1 = B[i][Dim-1];
        B[i][Dim-2] = M2[1]*t0 + M2[2]*t1;
        B[i][Dim-1] = M2[0]*t0 + M2[1]*t1;
    }
/*    for(i = 0;i < Dim;i++)
    {
    	t0 = Z[i][Dim-2];
        t1 = Z[i][Dim-1];
        Z[i][Dim-2] = M2[1]*t0 + M2[2]*t1;
        Z[i][Dim-1] = M2[0]*t0 + M2[1]*t1;
    }
*/
}

void EIqzimplicqr(matriz A, matriz B, unsigned int Dim)
{
	indice i,k,m,n;
    escalar r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11;
	escalar x,y,z,c,s,h;

    m = Dim-2;
    n = Dim-1;
    r0 = A[m][m] / B[m][m];
    r1 = A[0][0] / B[0][0];
    r2 = A[n][n] / B[n][n];
    r3 = A[m][n] / B[n][n];
    r4 = A[n][m] / B[m][m];
    r5 = B[m][n] / B[n][n];
    r6 = B[0][0] / A[1][0];
    r7 = A[0][1] / B[1][1];
    r8 = B[0][1] / B[1][1];
    r9 = A[1][1] / B[1][1];
    r10 = A[1][0] / B[0][0];
    r11 = A[2][1] / B[1][1];
	x = ((r0 - r1)*(r2 - r1) - r3*r4 + r4*r5*r1)*r6 + r7 - r1*r8;
    y = (r9 - r1) - r10*r8 - (r0 - r1) - (r2 - r1) + r4*r5;
    z = r11;

    for(k = 0;k < (Dim-2);k++)
    {
        h = x*x + z*z;
        c = x / h;
        s = z / h;
        OPgivenspre(A,k,k+2,c,s,Dim);
        OPgivenspre(B,k,k+2,c,s,Dim);
        if(k == 0) x = A[k][k];
        else x = A[k][k-1];
        h = x*x + y*y;
        c = x / h;
        s = y / h;
        OPgivenspre(A,k,k+1,c,s,Dim);
        OPgivenspre(B,k,k+1,c,s,Dim);
/*trazaIOputmat(A,"QA",'s',Dim);getchar();IOputmat(B,"QB",'s',Dim);getchar();
*/
        h = B[k+2][k]*B[k+2][k] + B[k+2][k+2]*B[k+2][k+2];
        c = B[k+2][k+2] / h;
        s = B[k+2][k] / h;
        OPgivenspost(A,k,k+2,c,s,Dim);
        OPgivenspost(B,k,k+2,c,s,Dim);
        h = B[k+2][k+1]*B[k+2][k+1] + B[k+2][k+2]*B[k+2][k+2];
        c = B[k+2][k+2] / h;
        s = B[k+2][k+1] / h;
        OPgivenspost(A,k+1,k+2,c,s,Dim);
        OPgivenspost(B,k+1,k+2,c,s,Dim);
/*trazaIOputmat(A,"QAZ1",'s',Dim);getchar();IOputmat(B,"QBZ1",'s',Dim);getchar();
*/
        h = B[k+1][k]*B[k+1][k] + B[k+1][k+1]*B[k+1][k+1];
        c = B[k+1][k+1] / h;
        s = B[k+1][k] / h;
        OPgivenspost(A,k,k+1,c,s,Dim);
        OPgivenspost(B,k,k+1,c,s,Dim);
/*trazaIOputmat(A,"QAZ1Z2",'s',Dim);getchar();IOputmat(B,"QBZ1Z2",'s',Dim);getchar();
*/
        x = A[k+1][k];
        y = A[k+2][k];
        if(k < (Dim-3)) z = A[k+3][k];
	}
    h = x*x + y*y;
    c = x / h;
    s = y / h;
    OPgivenspre(A,Dim-2,Dim-1,c,s,Dim);
    OPgivenspre(B,Dim-2,Dim-1,c,s,Dim);
/*trazaIOputmat(A,"QA",'s',Dim);getchar();IOputmat(B,"QB",'s',Dim);getchar();
*/
    h = B[Dim-1][Dim-2]*B[Dim-1][Dim-2] + B[Dim-1][Dim-1]*B[Dim-1][Dim-1];
    c = B[Dim-1][Dim-1] / h;
    s = B[Dim-1][Dim-2] / h;
    OPgivenspost(A,Dim-2,Dim-1,c,s,Dim);
    OPgivenspost(B,Dim-2,Dim-1,c,s,Dim);
/*trazaIOputmat(A,"QAZ",'s',Dim);getchar();IOputmat(B,"QBZ",'s',Dim);getchar();
*/
}

