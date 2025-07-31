/**
    Autor : Antonio Mesa
    Version : v 1.0
        1-XI-1994
        10-II-1995
**/

#include "nuevo.h"
#include <stdio.h>
#include <stdlib.h>
#include <process.h>

void ERfatal(unsigned int Coderr)
{
	switch (Coderr)
	{
		case ERRMEMO : puts("\n	!!! ERROR FATAL: MEMORIA INSUFICIENTE !!!"); break;
		case ERRFICH : puts("\n !!! ERROR FATAL: FICHERO INEXISTENTE !!!"); break;
		default : puts("\n !!! ERROR FATAL: CAUSA IMPREVISTA !!!");
	}
	exit(0);
}
