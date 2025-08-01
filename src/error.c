/**
    Autor : Antonio Mesa
    Version : v 1.0
        1-XI-1994
        10-II-1995
**/

#include <stdio.h>
#include <stdlib.h>

#include "defines.h"

void ERfatal(unsigned int Coderr)
{
	switch (Coderr)
	{
		case ERRMEMO : puts("\n	!!! ERROR FATAL: MEMORIA INSUFICIENTE !!!"); break;
		case ERRFICHRD : puts("\n !!! ERROR FATAL: READ FICHERO INEXISTENTE !!!"); break;
		case ERRFICHWR : puts("\n !!! ERROR FATAL: WRITE FICHERO INEXISTENTE !!!"); break;
		default : puts("\n !!! ERROR FATAL: CAUSA IMPREVISTA !!!");
	}
	exit(0);
}
