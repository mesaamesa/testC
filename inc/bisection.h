/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa  												*/
/*																		*/
/*	Prefijo mnemot‚cnico del m¢dulo : BI    							*/
/*	Modelo de compilaci¢n :	WINDOWS										*/
/* 	Ultima versi¢n : v 1.0										        */
/*	Descripci¢n : Implementa el método de la bisección para buscar  	*/
/*				un vap, utilizando la serie de Sturm.					*/
/*																		*/
/************************************************************************/

#ifndef BISEC
#define BISEC

void BImirasturm(unsigned int *, vector, unsigned int);
/*Pre : */
/*Post : */

void BIsturm(vector, matriz, matriz, escalar, unsigned int);
/*Pre : */
/*Post : */

void BIbuscainterval(unsigned int *, escalar *, escalar *, matriz, matriz, escalar, unsigned int);
/*Pre : */
/*Post : */

void BIbisec(escalar *, escalar, escalar, unsigned int, matriz, matriz, unsigned int);
/*Pre : */
/*Post : */

#endif
