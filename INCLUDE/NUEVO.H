/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa  												*/
/*																		*/
/*	Prefijo mnemot�cnico del m�dulo : 									*/
/*	Modelo de compilaci�n :	DOS											*/
/* 	Ultima versi�n : v 1.0										        */
/*	Descripci�n : definici�n de nuevos tipos de variables y de          */
/* 		constantes globales												*/
/*																		*/
/************************************************************************/

#ifndef NUEVO
#define NUEVO

#define FALSE 0
#define TRUE 1
#define MAXCHARS 20
#define TAMANOERR 4
#define KITERMAX 3
#define NDEC 12
#define ERROR 1
#define COTA 1.0e-5
#define EPSDIV 2
#define INFINITO 1.0e+300
#define EPS 1.0e-3

enum {NOERROR,ERRMEMO,ERRFTO,ERRFICHRD,ERRFICHWR,ERRITER};
enum {PROD,PRETRANS,POSTTRANS};
enum {CERO,POSITIVO,NEGATIVO};

typedef double escalar;
typedef double * vector;
typedef double ** matriz;
typedef unsigned int indice;
typedef unsigned char bool;
typedef char string[MAXCHARS];

#endif
