#include "ICMKernel\\random\icm_RandomKISS.h"
#include "ICMKernel\util\icm_macro.h"

#define ARM_ULONG_MAX 4294967295;

ICM_RandomKISS::ICM_RandomKISS() 
{ 
	Init(123456789,362436,521288629,7654321);
	its_a=698769069L;
}

ICM_RandomKISS::ICM_RandomKISS(unsigned long x, unsigned long y, unsigned long z, unsigned long c)
{
	Init(x,y,z,c);
}

ICM_RandomKISS::ICM_RandomKISS(const ICM_RandomKISS& ref):ICM_RandomGenerator(ref)
{

}

ARM_Object * ICM_RandomKISS::Clone(void)
{
	return new ICM_RandomKISS(*this);
}

ICM_RandomKISS& ICM_RandomKISS::operator=(const ICM_RandomKISS& ref)
{
	if (this!=&ref)
	{
		this->~ICM_RandomKISS(); 
		new(this)ICM_RandomKISS(ref);
	}
	return *this; 
}


ICM_RandomKISS::~ICM_RandomKISS() {}


void ICM_RandomKISS::Init(unsigned long x, unsigned long y, unsigned long z, unsigned long c)
{
		//0<=x<2^32, 0<y<2^32, 0<=z<2^32, 0<=c<698769069.
		if((x < 0) || (x> pow(2,32)))
			ICMTHROW(ERR_INVALID_ARGUMENT," must be  0<=x<2^32 ") ; 
		if((y < 1) || (y> pow(2,32)))
			ICMTHROW(ERR_INVALID_ARGUMENT," must be  0< y<2^32 ") ; 
		if((z < 0) || (z> pow(2,32)))
			ICMTHROW(ERR_INVALID_ARGUMENT," must be  0<=z<2^32 ") ; 
		if((c < 0 )|| (c > 698769068))
			ICMTHROW(ERR_INVALID_ARGUMENT," must be  0<=c<698769069 ") ; 
		its_x = x;
		its_y = y;
		its_z = z;
		its_c = c;
		its_a=698769069L;
}

void ICM_RandomKISS::GenerateRandoms(ARM_Vector& RandomVector)
{	
	GenerateRandomKISS(RandomVector);
}

void ICM_RandomKISS::GenerateRandomKISS(ARM_Vector& RandomVector)
{
	unsigned __int64 t=0;
	double * dVect = RandomVector.GetElt();
	int Dim = RandomVector.size();
	for (int i = 0 ; i<Dim; i++ ) {		
		its_x=69069*its_x+12345;
		its_y= pow(its_y,(its_y<<13)); 
		its_y= pow(its_y,(its_y>>17));
		its_y= pow(its_y,(its_y<<5));
		t=its_a*its_z+its_c;
		its_c=(t>>32);
		dVect[i]= (double)(its_x+its_y+(its_z=t))/ (double)ARM_ULONG_MAX;
	}
}


void ICM_RandomKISS::reset(){
	Init(123456789,362436,521288629,7654321);
}


void ICM_RandomKISS::setParameters(const std::string& sParamName, double dParamValue)
{
		if (sParamName == "X") {
			if((dParamValue < 0) || (dParamValue> pow(2,32)))
				ICMTHROW(ERR_INVALID_ARGUMENT," must be  0<=x<2^32 ") ; 
			its_x = dParamValue;
		}
		if (sParamName == "Y"){
			if((dParamValue < 1) || (dParamValue> pow(2,32)))
				ICMTHROW(ERR_INVALID_ARGUMENT," must be  0< y<2^32 ") ; 
			its_y = dParamValue; 
		}
		if (sParamName == "Z"){
			if((dParamValue < 0) || (dParamValue> pow(2,32)))
				ICMTHROW(ERR_INVALID_ARGUMENT," must be  0<=z<2^32 ") ; 
			its_z = dParamValue; 
		}
		if (sParamName == "C"){
			if(dParamValue < 0 || dParamValue > 698769068)
				ICMTHROW(ERR_INVALID_ARGUMENT," must be  0<=c<698769069 ") ; 
			its_c = dParamValue; 
		}

}


void ICM_RandomKISS::View(char* id, FILE* ficOut)
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

	fprintf(fOut, " Initial Seed X: %d \n",its_x);
	fprintf(fOut, " Initial Seed Y: %d \n",its_y);
	fprintf(fOut, " Initial Seed Z: %d \n",its_z);
	fprintf(fOut, " Initial Seed C: %d \n",its_c);

	fprintf(fOut, "\n");

	ICM_RandomGenerator::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}
/*
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB 
#undef NDIV
#undef EPS
#undef RNMX
*/