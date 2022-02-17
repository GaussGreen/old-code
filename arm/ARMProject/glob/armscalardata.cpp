#include "armscalardata.h"


int ARM_ScalarData::operator == (ARM_ScalarData& srcVal)
{
	return	( itsValue == srcVal.itsValue );
}


int ARM_ScalarData::operator < (ARM_ScalarData& srcVal)
{
	return	( itsValue < srcVal.itsValue );
}


void ARM_ScalarData::View(char* id, FILE* ficOut)
{
    FILE* fOut;

    char fOutName[200];
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
 
       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

    fprintf(fOut, "\n ARM_ScalarData");

    ARM_AbstractMarketClass::View(id, fOut);

	fprintf(fOut, "\n Value : %10.5lf\n\n", GetValue());

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}
