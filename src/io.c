/**
    Autor : Antonio Mesa
    Version : v 1.0
        1-XI-94
        10-II-1995
**/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
#include "error.h"

void IOmemmat(matriz *pA, unsigned int Dim)
{
	matriz A;
	indice i,j;

	A = (matriz)malloc(Dim*sizeof(vector));
	for(i = 0;i < Dim;i++)
    {
		A[i] = (vector)malloc(Dim*sizeof(escalar));
		for(j = 0;j < Dim;j++)
        	A[i][j] = 0.0;
	}
	*pA = A;
}

void IOmemmatmxn(pA,filA,colA)
matriz *pA;
unsigned int filA,colA;
{
    matriz A;
    indice i,j;

    A = (matriz)malloc(filA*sizeof(vector));
    for(i = 0;i < filA;i++)
    {
        A[i] = (vector)malloc(colA*sizeof(escalar));
        for(j = 0;j < colA;j++)
            A[i][j] = 0.0;
    }
    *pA = A;
}

void IOmemvect(vector *pV, unsigned int Dim)
{
	indice i;
	vector V;

	V = (vector)malloc(Dim*sizeof(escalar));
	for(i = 0;i < Dim;i++)
		V[i] = 0.0;
	*pV = V;
}

void IOfreemat(matriz A, unsigned int Dim)
{
	indice i;

	for(i = 0;i < Dim; i++)
		free(A[i]);
	free(A);
}

void IOfreematmxn(A,filA,colA)
matriz A;
unsigned int filA,colA;
{
    indice i;

    for(i = 0;i < filA; i++)
        free(A[i]);
    free(A);
}

void IOfreevect(vector V, unsigned int Dim)
{
	free(V);
}

void IOrdmat(matriz *pA, string Nombre, unsigned int *pDim)
{
	FILE *f;
	unsigned int ch,Dim;
	indice i,j,k;
	string buffer;
    matriz A;

	i = j = k = 0;
	if((f = fopen(Nombre,"r")) != NULL)
	{
    	while(((ch = getc(f)) != ',') && (!feof(f)))
		{
        	buffer[k] = ch;
            if(k < MAXCHARS) k++;
        }
        buffer[k] = '\0';
        Dim = atoi(buffer);
        IOmemmat(&A,Dim);
		while((i < Dim) && !feof(f))
		{
			k = 0;
			while(((ch = getc(f)) != ',') && (!feof(f)))
			{
				buffer[k] = ch;
				if(k < MAXCHARS) k++;
			}
			buffer[k] = '\0';
			A[i][j] = atof(buffer);
			j++;
			if(j >= Dim) {j = 0; i++;}
		}
	}
	fclose(f);
    *pA = A;
    *pDim = Dim;
	if(i != Dim) ERfatal(ERRFICHRD);
}

void IOwrmat(matriz A, string Nombre, unsigned int Dim)
{
	FILE *f;
	indice i,j;
	string buffer;

	i = j = 0;
	if((f = fopen(Nombre,"w")) != NULL)
	{
    	//TODO: itoa(Dim,buffer,10);
        fputs(buffer,f);
        putc(',',f);
        putc('\n',f);
		while((i < Dim) && !feof(f))
		{
			gcvt(A[i][j],NDEC,buffer);
			fputs(buffer, f);
			putc(',', f);
			j++;
			if(j == Dim)
			{
				i++;
				j = 0;
				putc('\n', f);
			}
		}
	}
	fclose(f);
	if(i != Dim) ERfatal(ERRFICHWR);
}

void IOrdvect(vector *pV, string Nombre, unsigned int *pDim)
{
	FILE *f;
	unsigned int ch,Dim;
	indice i,k;
	string buffer;
    vector V;

	i = k = 0;
	if((f = fopen(Nombre,"r")) != NULL)
	{
    	while(((ch = getc(f)) != ',') && (!feof(f)))
        {
        	buffer[k] = ch;
            if(k < MAXCHARS) k++;
        }
        buffer[k] = '\0';
        Dim = atoi(buffer);
        IOmemvect(&V,Dim);
		while((i < Dim) && !feof(f))
		{
			k = 0;
			while(((ch = getc(f)) != ',') && (!feof(f)))
			{
				buffer[k] = ch;
				if(k < MAXCHARS) k++;
			}
			buffer[k] = '\0';
			V[i] = atof(buffer);
			i++;
		}
	}
	fclose(f);
    *pV = V;
    *pDim = Dim;
	if(i != Dim) ERfatal(ERRFICHRD);
}

void IOwrvect(vector V, string Nombre, unsigned int Dim)
{
	FILE *f;
	indice i;
	string buffer;

	i = 0;
	if((f = fopen(Nombre,"w")) != NULL)
	{
		while((i < Dim) && !feof(f))
		{
			gcvt(V[i],NDEC,buffer);
			fputs(buffer, f);
			putc(' ', f);
			i++;
		}
	}
	fclose(f);
	if(i != Dim) ERfatal(ERRFICHWR);
}

void IOwrvectCo(Vr,Vi,Nombre,Dim)
vector Vr,Vi;
string Nombre;
unsigned int Dim;
{
    FILE *f;
    indice i;
    string buffer;

    i = 0;
    if((f = fopen(Nombre,"w")) != NULL)
    {
        while((i < Dim) && !feof(f))
        {
            gcvt(Vr[i],NDEC,buffer);
            fputs(buffer,f);
            if(Vi[i] >= 0.0)
            {
                fputs(" + i",f);
                gcvt(Vi[i],NDEC,buffer);
            }
            else
            {
                fputs(" - i",f);
                gcvt(-Vi[i],NDEC,buffer);
            }
            fputs(buffer,f);
            fputc(' ',f);
            fputc('\n',f); 
            i++;
        }
    }
    else ERfatal(ERRFICHWR);
    fclose(f);
    if(i != Dim) ERfatal(ERRFTO);
}

void IOwrvect2(vector V, string Nombre, unsigned int Dim)
{
	FILE *f;
	indice i;
	string buffer;

	i = 0;
	if((f = fopen(Nombre,"w")) != NULL)
	{
    	//TODO: itoa(Dim,buffer,10);
        fputs(buffer,f);
        putc(',',f);
        putc('\n',f);
		while((i < Dim) && !feof(f))
		{
			gcvt(V[i],NDEC,buffer);
			fputs(buffer, f);
			putc(',', f);
			i++;
			putc('\n', f);
		}
	}
	fclose(f);
	if(i != Dim) ERfatal(ERRFICHWR);
}

void IOputmat(matriz A, string Nombre, char Formato, unsigned int Dim)
{
	indice i,j;

	for(j = 0;j < Dim; j++)
		for(i = 0;i < Dim; i++)
			if(Formato == 'l')
				printf("\n%s (%d)(%d) = %1.17g",Nombre,i,j,A[i][j]);
			else
				printf("\n%s (%d)(%d) = %1.4g",Nombre,i,j,A[i][j]);
}

void IOputvect(vector V, string Nombre, char Formato, unsigned int Dim)
{
	indice i;

	for(i = 0;i < Dim; i++)
		if(Formato == 'l')
			printf("\n%s (%d) = %1.17g",Nombre,i,V[i]);
			else
				printf("\n%s (%d) = %1.4g",Nombre,i,V[i]);
}

void IOwrpar(matriz X, vector L, string VEP, string VAP, unsigned int Dim)
{
	FILE *fvap,*fvep;
	indice i,j;
	string buffer,fichVAP,fichVEP;

	for(i = 0;i < Dim;i++)
	{
		strcpy(fichVAP,VAP);
        strcpy(fichVEP,VEP);
		//TODO: itoa(i,buffer,10);
		strcat(fichVAP,buffer);
        strcat(fichVEP,buffer);
		if(((fvap = fopen(fichVAP,"w")) != NULL)
        	&& ((fvep = fopen(fichVEP,"w")) != NULL))
		{
			gcvt(L[i],NDEC,buffer);
            fputs(buffer,fvap);
            j = 0;
			while((j < Dim) && !feof(fvep))
			{
				gcvt(X[j][i],NDEC,buffer);
				fputs(buffer,fvep);
				putc(' ', fvep);
				j++;
			}
		}
		fclose(fvap);
        fclose(fvep);
		if(j != Dim) ERfatal(ERRFICHWR);
	}
}
