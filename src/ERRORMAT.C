/************************************************************************/
/*																		*/
/*	Autor : Antonio Mesa												*/
/*																		*/
/*	Versi¢n : 															*/
/*		   		v 1.0		1-XI-1994			10-II-1995				*/
/*																		*/
/************************************************************************/

#include "nuevo.h"
#include <stdio.h>
#include <stdlib.h>
#include <process.h>

void ERfatal(Coderr)
unsigned int Coderr;
{
	switch (Coderr)
	{
		case ERRMEMO:
			printf("Error Memo");
		break;
		default:
			printf("Error n° %d",Coderr);
		break;
	}
}