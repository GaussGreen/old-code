#include "ICMKernel\\random\icm_RandomRan1.h"
#include <nag.h>
#include <nagg05.h>
#include "ICMKernel\util\icm_macro.h"

#define NTAB 32

ICM_RandomRan1::ICM_RandomRan1() 
{ 
	Init(-2) ;
}

ICM_RandomRan1::ICM_RandomRan1(long InitialSeed) 
{
	Init(InitialSeed) ;
}

ICM_RandomRan1::ICM_RandomRan1(const ICM_RandomRan1& ref):ICM_RandomGenerator(ref)
{
	itsInitialSeed = ref.itsInitialSeed; 
	itsCurrentSeed = ref.itsCurrentSeed; 
	its_iy=ref.its_iy;
	for (int i=0; i<NTAB; i++) 
	its_iv[i] = ref.its_iv[i];
}

ARM_Object * ICM_RandomRan1::Clone(void)
{
	return new ICM_RandomRan1(*this);
}

ICM_RandomRan1& ICM_RandomRan1::operator=(const ICM_RandomRan1& ref)
{
	if (this!=&ref)
	{
		this->~ICM_RandomRan1(); 
		new(this)ICM_RandomRan1(ref);
	}
	return *this; 
}


ICM_RandomRan1::~ICM_RandomRan1() {}


void ICM_RandomRan1::Init(long InitialSeed)
{
		if(InitialSeed > 0)
			ICMTHROW(ERR_INVALID_ARGUMENT,"seed must be <0 ") ; 
		itsInitialSeed = InitialSeed;
		itsCurrentSeed = InitialSeed;
		its_iy=0;
}




void ICM_RandomRan1::reset(){
	itsCurrentSeed = itsInitialSeed; 	
}


void ICM_RandomRan1::setParameters(const std::string& sParamName, double dParamValue)
{
	if (sParamName == "InitialSeed")
		itsInitialSeed = dParamValue; 
	if (sParamName == "CurrentSeed")
		itsCurrentSeed = dParamValue; 

}


void ICM_RandomRan1::View(char* id, FILE* ficOut)
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

	fprintf(fOut, "\t\t\t ----------------- Random Generator Ran1 ----------------- \n\n");

	fprintf(fOut, " Initial Seed : %d \n",itsInitialSeed);
	fprintf(fOut, " Current Seed : %d \n",itsCurrentSeed);
//	fprintf(fOut, " Random nb : %f \n",itsRandom);
	fprintf(fOut, "\n");

	ICM_RandomGenerator::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}
/*

*/