
#include "ICMKernel\inst\icm_gen.h"

void ICM_GenCF::Init()
{
	SetName(ICM_GENCF);

	itsMatrix = NULL;
}


ICM_GenCF::ICM_GenCF (ICM_Matrix<ARM_Vector>* matrice)
{
	Init();

	Set(matrice);
	
}

void ICM_GenCF::Set(ICM_Matrix<ARM_Vector>* matrice)
{
	SetName(ICM_GENCF);

	if (itsMatrix)
		delete itsMatrix;

	itsMatrix = (ICM_Matrix<ARM_Vector>*) matrice->Clone();
}


void ICM_GenCF::BitwiseCopy(const ARM_Object* src)
{
    ICM_GenCF* gen = (ICM_GenCF*) src;

	if (gen->itsMatrix)
	{
		if (itsMatrix)
			delete itsMatrix;
		itsMatrix = (ICM_Matrix<ARM_Vector>*) gen->itsMatrix->Clone();
	}
}

void ICM_GenCF::Copy(const ARM_Object* src)
{
     ARM_Security::Copy(src);
 
     BitwiseCopy(src);
}


ARM_Object* ICM_GenCF::Clone(void)
{
     ICM_GenCF* theClone = new ICM_GenCF();

     theClone->Copy(this);
 
     return(theClone);
}

// *************************************************************
// View Matrix 
// *************************************************************
void ICM_GenCF::View(char* id, FILE* ficOut)
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

	// Affichage de la matrice
	fprintf(fOut, "\t\t\t ----------------- Pricing Matrix ----------------- \n");

	if (GetMatrix())
		GetMatrix()->View(id, fOut);

	if ( ficOut == NULL )
	{
	fclose(fOut);
	}
}
