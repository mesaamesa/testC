#include "nuevo.h"
#include "drveig.h"
#include <stdlib.h>
#include <math.h>

int main(int argcnt, char **argval)
{
    DRdrivereig(argval[1],argval[2],(bool)atoi(argval[3]),atof(argval[4]));
}