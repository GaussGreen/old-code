#include "ICMKernel\\random\icm_randomLaws.h"
#include "ICMKernel\util\icm_macro.h"
#include "ICMKernel\util\icm_utils.h"

ICM_RandomLaws::ICM_RandomLaws() 
{ 
	Init();
}

ICM_RandomLaws::ICM_RandomLaws(const ICM_RandomGenerator& pRandomUnif)
{
	Init();
	its_pRandomUnif = dyn_clone(&pRandomUnif); // pvUnifRand cannot be null
}

ICM_RandomLaws::ICM_RandomLaws(const ICM_RandomLaws& ref):ICM_RandomGenerator(ref)
{
	its_pRandomUnif = dyn_clone((ref.its_pRandomUnif));
}

ICM_RandomLaws::~ICM_RandomLaws() 
{
	if (its_pRandomUnif) 
		delete its_pRandomUnif;
	its_pRandomUnif = NULL;
}

void ICM_RandomLaws::Init()
{
	its_pRandomUnif = NULL;
}
	
ARM_Object *  ICM_RandomLaws::Clone(void)
{
	return new ICM_RandomLaws(*this);
}

const ICM_RandomGenerator* ICM_RandomLaws::getUniformGenerator() const {
	return its_pRandomUnif;
}

void ICM_RandomLaws::SetUniformGenerator(const ICM_RandomGenerator& pRandomUnif)
{
	if(its_pRandomUnif) delete its_pRandomUnif;
	its_pRandomUnif = dyn_clone(&pRandomUnif);
}

void ICM_RandomLaws::View(char* id, FILE* ficOut)
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
	ICM_RandomGenerator::View(id, fOut);
	fprintf(fOut, "\t\t\t ----------------- Random Generator Laws ----------------- \n\n");
	fprintf(fOut, "\n");
	
	if (its_pRandomUnif) {
		its_pRandomUnif->View("",fOut);
		
			
	} else {
		fprintf(fOut, "No uniform random set\n");
	}

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}
