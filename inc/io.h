/**
    Autor : Antonio Mesa
    Prefijo mnemotecnico del modulo : IO
    Modelo de compilacion : WINDOWS
    Ultima version : v 1.0
    Descripcion : Comunica las matrices y los resultados con pantalla, ficheros, etc.
**/

#ifndef IOMATRIZ
#define IOMATRIZ

void IOrdmat(matriz *, string, unsigned int *);
/* Pre : */
/* Post : */

void IOwrmat(matriz, string, unsigned int);
/* Pre : */
/* Post : */

void IOputmat(matriz, string, char, unsigned int);
/* Pre : */
/* Post : */

void IOmemmat(matriz *, unsigned int);
/* Pre : */
/* Post : */

void IOmemmatmxn();
/* Pre : */
/* Post : */

void IOfreemat(matriz, unsigned int);
/* Pre : */
/* Post : */

void IOfreematmxn();
/* Pre : */
/* Post : */

void IOrdvect(vector *, string, unsigned int *);
/* Pre : */
/* Post : */

void IOwrvectRe();
/* Pre : */
/* Post : */

void IOwrvectCo();
/* Pre : */
/* Post : */

void IOwrvect(vector, string, unsigned int);
/* Pre : */
/* Post : */

void IOwrvect2(vector, string, unsigned int);
/* Pre : */
/* Post : */

void IOputvect(vector, string, char, unsigned int);
/* Pre : */
/* Post : */

void IOmemvect(vector *, unsigned int);
/* Pre : */
/* Post : */

void IOfreevect(vector, unsigned int);
/* Pre : */
/* Post : */

void IOwrpar(matriz, vector, string, string, unsigned int);
/* Pre : */
/* Post : */

#endif
