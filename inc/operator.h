/**
    Autor : Antonio Mesa
    Prefijo mnemotecnico del modulo : OP
    Modelo de compilacion :    WINDOWS
    Ultima version : v 1.0
    Descripcion : Realiza las operaciones entre matrices, vectores y escalares.
**/

#ifndef OPMATRIZ
#define OPMATRIZ

escalar OPeps();
/* Pre : */
/* Post : */

void OPsumamat(matriz, matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void OPsumasubmat();
/* Pre : */
/* Post : */

void OPrestamat(matriz, matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void OPrestasubmat();
/* Pre : */
/* Post : */

void OPprodmat(matriz, matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void OPprodmatmxn();
/* Pre : */
/* Post : */

void OPprodsubmat();
/* Pre : */
/* Post : */

void OPcpymat(matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void OPsumavect(vector, vector, vector, unsigned int);
/* Pre : */
/* Post : */

void OPrestavect(vector, vector, vector, unsigned int);
/* Pre : */
/* Post : */

escalar OPprodint(vector, vector, unsigned int);
/* Pre : */
/* Post : */

escalar OPprodintmat(vector, matriz, vector, unsigned int);
/* Pre : */
/* Post : */

escalar OPnormamaxmat(matriz, unsigned int);
/* Pre : */
/* Post : */

void OPctemat(matriz, matriz, escalar, unsigned int);
/* Pre : */
/* Post : */

void OPctevect(vector, vector, escalar, unsigned int);
/* Pre : */
/* Post : */

void OPctevect(   );
/* Pre : */
/* Post : */

void OPcpyvect(vector, vector, unsigned int);
/* Pre : */
/* Post : */

void OPmatdiag(matriz, matriz, vector, unsigned int);
/* Pre : */
/* Post : */

void OPdiagmat(matriz, matriz, vector, unsigned int);
/* Pre : */
/* Post : */

void OPmatvect(vector, matriz, vector, unsigned int);
/* Pre : */
/* Post : */

void OPinvtsup(matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void OPinvtinf(matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void OPinvdiag(vector, vector, unsigned int);
/* Pre : */
/* Post : */

double OPnorma2(vector, unsigned int);
/* Pre : */
/* Post : */

double OPnormamax(vector, unsigned int);
/* Pre : */
/* Post : */

void OPfiltramat(matriz, escalar, unsigned int);
/* Pre : */
/* Post : */

void OPtransmat(matriz, unsigned int);
/* Pre : */
/* Post : */

void OPsisteqsup(matriz, vector, vector, unsigned int);
/* Pre : */
/* Post : */

void OPsisteqinf(matriz, vector, vector, unsigned int);
/* Pre : */
/* Post : */

void OPidentidad(matriz, unsigned int);
/* Pre : */
/* Post : */

void OPceromat(matriz, unsigned int);
/* Pre : */
/* Post : */

void OPunovect(vector, unsigned int);
/* Pre : */
/* Post : */

void OPcerovect(vector, unsigned int);
/* Pre : */
/* Post : */

void OPsort(vector, unsigned int);
/* Pre : */
/* Post : */

void OPgivenspre(matriz, indice, indice, escalar, escalar, unsigned int);
/* Pre : */
/* Post : */

void OPgivenspost(matriz, indice, indice, escalar, escalar, unsigned int);
/* Pre : */
/* Post : */

void OPelempre1(matriz, indice, indice, escalar, unsigned int);
/* Pre : */
/* Post : */

void OPelempre2(matriz, indice, indice, escalar, unsigned int);
/* Pre : */
/* Post : */

void OPelempost1(matriz, indice, indice, escalar, unsigned int);
/* Pre : */
/* Post : */

void OPelempost2(matriz, indice, indice, escalar, unsigned int);
/* Pre : */
/* Post : */

void OPcambiafil(matriz, indice, indice, unsigned int);
/* Pre : */
/* Post : */

void OPcambiacol(matriz, indice, indice, unsigned int);
/* Pre : */
/* Post : */

void OPsubdiaghess(matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

#endif
