#include "ICMKernel\random\icm_RandomGenerator.h"


void ICM_RandomGenerator::Init()
{
}

ICM_RandomGenerator::ICM_RandomGenerator() 
{ 
	Init() ;
}

ICM_RandomGenerator::ICM_RandomGenerator(const ICM_RandomGenerator& ref):ARM_Object(ref)
{ 
}


ICM_RandomGenerator::~ICM_RandomGenerator()
{}

void ICM_RandomGenerator::View(char* id, FILE* ficOut)
{

	FILE* fOut;
	char  fOutName[200];

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

	int size =0;

	fprintf(fOut, "\t\t\t ----------------- Basis Random Generator ----------------- \n\n");

	
	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}
