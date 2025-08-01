/**
    Autor : Antonio Mesa
    Prefijo mnemotecnico del modulo : FC
    Modelo de compilacion :    WINDOWS
    Ultima version : v 1.0
    Descripcion : Convierte una matriz cuadrada en una matrices de una estructura determinada.
**/

#ifndef FCMATRIZ
#define FCMATRIZ

void FCldu(matriz, matriz, vector, matriz, unsigned int);
/* Pre : */
/* Post : */

void FClu(matriz, matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void FCchol();
/* Pre : */
/* Post : */

void FCqrgram(matriz, matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void FCqrhouse(matriz, matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void FCqrgivens();
/* Pre : */
/* Post : */

void FClrelem(matriz, matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void FChesstriorto(matriz, matriz, matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void FChesstriorto2(matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void FChesstrielem(matriz, matriz, matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

void FChesstrielem2(matriz, matriz, unsigned int);
/* Pre : */
/* Post : */

#endif
