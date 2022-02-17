
#include "armcalypso.h"


int ISCALYPSOVERSION = 0;


void SwitchToCALYPSO()
{
	ISCALYPSOVERSION = 1;
}


int GetCALYPSOVersion()
{
	return ISCALYPSOVERSION;
}
