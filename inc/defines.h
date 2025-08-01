/**
    Autor : Antonio Mesa
    Prefijo mnemotecnico del modulo :
    Modelo de compilacion :    DOS
    Ultima version : v 1.0
    Descripcion : definicion de nuevos tipos de variables y de constantes globales.
 **/

#ifndef NUEVO
#define NUEVO

#define FALSE 0
#define TRUE 1
#define MAXCHARS 20
#define TAMANOERR 4
#define KITERMAX 1
#define NDEC 12
#define ERROR 1
#define COTA 1.0e-5
#define EPSDIV 2
#define INFINITO 1.0e+300
#define EPS 1.0e-15
enum {NOERROR,ERRMEMO,ERRFTO,ERRFICHRD,ERRFICHWR,ERRITER};
enum {PROD,PRETRANS,POSTTRANS};
/*constantes bisec*/
enum {CERO,POSITIVO,NEGATIVO};
#define COTASUP 10.0e10

typedef double escalar;
typedef double * vector;
typedef double ** matriz;
typedef unsigned int indice;
typedef unsigned char bool;
typedef char string[MAXCHARS];

#endif
