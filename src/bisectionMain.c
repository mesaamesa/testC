#include <stdlib.h>

#include "drvBisection.h"

int main(int argcnt, char **argval)
{
	DRdriverbisec(argval[1],argval[2],atof(argval[3]));
}