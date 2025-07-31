/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa  												*/
/*																		*/
/*	Prefijo mnemot‚cnico del m¢dulo : QZ    							*/
/*	Modelo de compilaci¢n :	WINDOWS										*/
/* 	Ultima versi¢n : v 1.0										        */
/*	Descripci¢n : Calcula los vap's aplicando el método QZ explicito	*/
/*				y los vep's mediante iteración inversa.					*/
/*																		*/
/************************************************************************/

#ifndef QZEXPLIC
#define QZEXPLIC

void QZexplicito(matriz, matriz, unsigned int);
/*Pre : */
/*Post : */

void QZvepslu(matriz, matriz, matriz, vector, unsigned int);
/*Pre : */
/*Post : */

void QZveplu1(vector, matriz, matriz, escalar, unsigned int);
/*Pre : */
/*Post : */

#endif