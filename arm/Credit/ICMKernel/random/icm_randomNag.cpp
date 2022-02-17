#include "ICMKernel\\random\icm_RandomNag.h"

#include "ICMKernel\util\icm_macro.h"


ICM_RandomNag::ICM_RandomNag() 
{ 
	Init(0,false) ;
}

ICM_RandomNag::ICM_RandomNag(long InitialSeed, bool IsRandomSeed) 
{ 
	Init(InitialSeed,IsRandomSeed) ;
}

ICM_RandomNag::ICM_RandomNag(const ICM_RandomNag& ref):ICM_RandomGenerator(ref)
{
	itsInitialSeed = ref.itsInitialSeed; 
	itsIsRandomSeed = ref.itsIsRandomSeed; 

}

ARM_Object * ICM_RandomNag::Clone(void)
{
	return new ICM_RandomNag(*this);
}

ICM_RandomNag& ICM_RandomNag::operator=(const ICM_RandomNag& ref)
{
	if (this!=&ref)
	{
		this->~ICM_RandomNag(); 
		new(this)ICM_RandomNag(ref);
	}
	return *this; 
}


ICM_RandomNag::~ICM_RandomNag(){}


void ICM_RandomNag::Init(long InitialSeed, bool IsRandomSeed)
{
		itsInitialSeed = InitialSeed;
		itsIsRandomSeed = IsRandomSeed;
		
		if (itsIsRandomSeed) 
			g05ccc(); 	
		else 
			g05cbc(itsInitialSeed); 	
}


void ICM_RandomNag::reset()
{
		if (itsIsRandomSeed) 
			g05ccc(); 	
		else 
			g05cbc(itsInitialSeed); 	
}


void ICM_RandomNag::setParameters(const std::string& sParamName, double dParamValue)
{
	if (sParamName == "InitialSeed")
		itsInitialSeed = dParamValue; 
}


void ICM_RandomNag::View(char* id, FILE* ficOut)
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

	fprintf(fOut, "\t\t\t ----------------- Random Generator Ran Nag ----------------- \n\n");

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