/**
    Autor : Antonio Mesa
    Prefijo mnemotecnico del modulo : QZ
    Modelo de compilacion :    WINDOWS
    Ultima version : v 1.0
    Descripcion : Calcula los vap's aplicando el metodo QZ explicito y los vep's mediante iteracion inversa.
**/

#ifndef QZEXPLIC
#define QZEXPLIC

void QZexplicito(matriz, matriz, unsigned int);
/*Pre : */
/*Post : */

void QZexplic2();
/*Pre : */
/*Post : */

void QZvepslu(matriz, matriz, matriz, vector, unsigned int);
/*Pre : */
/*Post : */

void QZveplu1(vector, matriz, matriz, escalar, unsigned int);
/*Pre : */
/*Post : */

#endif