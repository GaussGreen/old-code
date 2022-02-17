#include "ICMKernel\util\icm_RangeFactor.h"

// ----------------------------------------------------------------------
// Constructor
// ----------------------------------------------------------------------
void ICM_RangeFactor::Init()
{
	itsnbRow = 0;
	itsIndice = 0;
	itsVInf.clear();
	itsVMax.clear();
	itsVValueMin.clear();
	itsVValueMax.clear();
}

void ICM_RangeFactor::Set(const vector<double>& inf,
						  const vector<double>& max,
						  const vector<double>& valueMin,
						  const vector<double>& valueMax)
{
	int sizeInf = inf.size();
	int sizeMax = max.size();
	int sizeValueMin = valueMin.size();
	int sizeValueMax = valueMax.size();
	if ( (sizeInf != sizeMax)		||
		 (sizeInf != sizeValueMin)	||
		 (sizeInf != sizeValueMax)
	   )
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_RangeFactor Ctor: wrong size of vector"); 
	}

	itsnbRow = sizeInf;

	itsVInf.resize(itsnbRow);			itsVInf = inf;
	itsVMax.resize(itsnbRow);			itsVMax = max;
	itsVValueMin.resize(itsnbRow);		itsVValueMin = valueMin;
	itsVValueMax.resize(itsnbRow);		itsVValueMax = valueMax;

	computeMid();
}
// ----------------------------------------------------------------------
//	View Method
// ----------------------------------------------------------------------
void ICM_RangeFactor::View(char* id, FILE* ficOut)
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

	fprintf(fOut, "\n ======> Range Factor :\n\n");

	fprintf(fOut, "\n Nombre de ligne: %i \n\n",itsnbRow);
		
	fprintf(fOut, " Inf\t");
	fprintf(fOut, " Max\t");
	fprintf(fOut, " ValueMin\t");
	fprintf(fOut, " ValueMax\t");
	fprintf(fOut, " ValueMid\t");
	
	fprintf(fOut, "\n\n");
			
	for (int i = 0; i<itsnbRow; i++)
	{	
		fprintf(fOut, " %.2f\t",itsVInf[i]);
		fprintf(fOut, " %.2f\t",itsVMax[i]);
		fprintf(fOut, " %.3f\t",itsVValueMin[i]);
		fprintf(fOut, " %.2f\t",itsVValueMax[i]);
		fprintf(fOut, " %.2f\t\n",itsVValueMid[i]);
	}

	fprintf(fOut, "\n");
	
	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}
// -----------------------------------------------------------------------
// FindValue Method : prend en input un x et selon la valeur de l'intervalle
// dans lequel on se trouve renvoie la value Mid qui correspond
// -----------------------------------------------------------------------
int ICM_RangeFactor::FindIndex(double x)
{
	// On recupere l'indice de l'intervalle dans lequel on se trouve
	unsigned int indice = itsIndice;

	// Cas ou on sort de l'intervalle
		// Par le haut
	if (x>itsVMax[indice])
	{
		while ((indice < itsnbRow) && (x>itsVMax[indice]))
			indice ++;

		if (indice == itsnbRow) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_RangeFactor FindIndex: value out of range ("<<indice<<")");	

		SetIndice(indice - 1);
		return (indice - 1);
	}
	else if (x < itsVInf[indice])
	{
		while ((indice >= 0) && (x < itsVInf[indice]))
			indice --;

		if (indice < 0)  
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_RangeFactor FindIndex: value out of range ("<<indice<<")");	

		SetIndice(indice + 1);
		return (indice + 1);
	}
	else 
		return (indice);	
}

// -----------------------------------------------------------------------
// Compute Mid
// -----------------------------------------------------------------------
void ICM_RangeFactor::computeMid(void)
{
	itsVValueMid.resize(itsnbRow);
	for (int i = 0; i<itsnbRow; i++)
		itsVValueMid[i] = (itsVValueMin[i] + itsVValueMax[i])/2.0;
}