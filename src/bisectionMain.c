#include "nuevo.h"
#include "drvbisec.h"
#include <stdlib.h>
#include <math.h>

int main(int argcnt, char **argval)
{
	DRdriverbisec(argval[1],argval[2],atof(argval[3]));
}