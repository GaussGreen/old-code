#include "ICMKernel\\random\icm_RandomRanDef.h"
#include <nag.h>
#include <nagg05.h>
#include <stdlib.h>
#include <time.h>
#include "ICMKernel\util\icm_macro.h"

ICM_RandomRanDef::ICM_RandomRanDef() 
{ 
	Init(1,false) ;
}

ICM_RandomRanDef::ICM_RandomRanDef(int InitialSeed, bool IsRandomSeed) 
{
	Init(InitialSeed,IsRandomSeed) ;
}

ICM_RandomRanDef::ICM_RandomRanDef(const ICM_RandomRanDef& ref):ICM_RandomGenerator(ref)
{
	itsInitialSeed = ref.itsInitialSeed; 
	itsIsRandomSeed = ref.itsIsRandomSeed; 
	
}

ARM_Object * ICM_RandomRanDef::Clone(void)
{
	return new ICM_RandomRanDef(*this);
}

ICM_RandomRanDef& ICM_RandomRanDef::operator=(const ICM_RandomRanDef& ref)
{
	if (this!=&ref)
	{
		this->~ICM_RandomRanDef(); 
		new(this)ICM_RandomRanDef(ref);
	}
	return *this; 
}


ICM_RandomRanDef::~ICM_RandomRanDef() {}

/*void ICM_RandomRanDef::setInitialSeed(const int i)  
{  
	if(!(itsInitialSeed >0) )
		ICMTHROW(ERR_INVALID_ARGUMENT,"seed must be positive") ; 
	itsInitialSeed = i;
	srand(itsInitialSeed);
}

void ICM_RandomRanDef::setIsRandomSeed(const bool b)  
{ 
	
	itsIsRandomSeed = b;
	if(itsIsRandomSeed) 
		itsInitialSeed = time(0);
	srand(itsInitialSeed);
}
*/

void ICM_RandomRanDef::Init(int InitialSeed, bool IsRandomSeed)
{
		if(!(InitialSeed >0) )
			ICMTHROW(ERR_INVALID_ARGUMENT,"seed must be positive") ; 
		itsInitialSeed = InitialSeed;
		itsIsRandomSeed = IsRandomSeed;
		
		if (itsIsRandomSeed) 
			itsInitialSeed = time(0);

		srand(itsInitialSeed);

}




void ICM_RandomRanDef::reset(){
		//if (itsIsRandomSeed) 
		//	itsInitialSeed = time(0);

		srand(itsInitialSeed);
}


void ICM_RandomRanDef::setParameters(const std::string& sParamName, double dParamValue)
{
	if (sParamName == "InitialSeed")
		itsInitialSeed = dParamValue; 
}


void ICM_RandomRanDef::View(char* id, FILE* ficOut)
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

	fprintf(fOut, "\t\t\t ----------------- Random Generator RanDef ----------------- \n\n");

	fprintf(fOut, " Initial Seed : %d \n",itsInitialSeed);
	fprintf(fOut, " Random Seed : %d \n",itsIsRandomSeed);
//	fprintf(fOut, " Random nb : %f \n",itsRandom);
	fprintf(fOut, "\n");

	ICM_RandomGenerator::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}