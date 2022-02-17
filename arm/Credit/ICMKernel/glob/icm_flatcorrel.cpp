#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\glob\icm_flatcorrel.h"
#include "ICMKernel\util\icm_utils.h"



ICM_FlatCorrel::ICM_FlatCorrel() {
	Init();
}
//	----------------------------------------------------------------------------------------------------
ICM_FlatCorrel::ICM_FlatCorrel(const ARM_Date& AsOf,
							   const std::string& structName,
							   const ARM_IRIndex* ix1, 
							   const ARM_IRIndex* ix2,
							   double correlValue)
{
	Init();
	Set(AsOf,structName,ix1,ix2,correlValue);
}
// --------------------------------------------------------------------------------
void ICM_FlatCorrel::Init(void) {
	SetName(ICM_FLAT_CORRELATION);
}

// -----------------------------------------------------------------------------------
// FIXMEFRED: mig.vc8 (28/05/2007 10:24:05):missing return type
void ICM_FlatCorrel::Set(const ARM_Date& AsOf,
					const std::string& structName,
					const ARM_IRIndex* ix1, 
					const ARM_IRIndex* ix2,
					double correlValue) {
	SetAsOf(AsOf);
	SetIndex1(ix1);
	SetIndex2(ix2);
	//SetStructName
	if (gt(correlValue,1.) || lt(correlValue,-1.) )
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_FlatCorrel: Can't create with correl="<<correlValue) ; 
	itsCorrel=correlValue; 
}
//	-----------------------------------------------------------------------------
ICM_FlatCorrel::ICM_FlatCorrel(const ICM_FlatCorrel&ref) : ICM_Correlation(ref)
{
	itsCorrel=ref.itsCorrel ;
}
//	-----------------------------------------------------------------------------
ICM_FlatCorrel& ICM_FlatCorrel::operator=(const ICM_FlatCorrel& ref)
{
	if (&ref!=this) 
	{
		this->~ICM_FlatCorrel(); 
		new(this)ICM_FlatCorrel(ref); 
	}
	return *this; 
}
//	-----------------------------------------------------------------------------
ICM_FlatCorrel::~ICM_FlatCorrel() 
{} 
//	-----------------------------------------------------------------------------
// virtual 
ARM_Object* ICM_FlatCorrel::Clone()
{
	ICM_FlatCorrel* theClone = new ICM_FlatCorrel();
	theClone->Copy(this);
	return(theClone);
}
//	-----------------------------------------------------------------------------
// virtual 
ARM_Object* ICM_FlatCorrel::Clone() const 
{
	return (unconst(this))->Clone(); 
}
// -------------
//	Copy Method 
// -------------
void ICM_FlatCorrel::Copy(const ARM_Object* src)
{
	ICM_Correlation::Copy(src);
	BitwiseCopy(src);
}

// ----------------------------
//	Copy of members data
// ----------------------------
void ICM_FlatCorrel::BitwiseCopy(const ARM_Object* src)
{
	const ICM_FlatCorrel* Correlation = (const ICM_FlatCorrel *) src;
	itsCorrel = Correlation->itsCorrel;
}

//	-----------------------------------------------------------------------------
// virtual 
ICM_Correlation* ICM_FlatCorrel::GenerateShiftCorrel(const std::string& label,
											 qSENSITIVITY_TYPE typesensi,
											 double epsilon )
{
	
	return ICM_FlatCorrel::GenerateShiftBetas("",typesensi,epsilon); 
}
//	-----------------------------------------------------------------------------
// virtual 
ICM_Correlation* ICM_FlatCorrel::GenerateShiftBetas(const std::string& label,
											qSENSITIVITY_TYPE typesensi,
											double epsilon  )
{
	ICM_FlatCorrel* item=0 ;
	switch(typesensi)
	{
	case ICMBETA_TYPE: 
		{
			double beta = sqrt(fabs(itsCorrel)) + epsilon; 
			double correl = beta*beta ;
			if (leq(correl,0.)) correl*=-1; 
			if (gt(correl,1.)) correl=1; 
			else if (lt(correl,1.)) correl=-1; 
			item = new ICM_FlatCorrel(GetAsOfDate(),GetStructName(),
				GetIndex1(),GetIndex2(),correl); 
		} 
		break ;
	case ICMCORRELATION_TYPE:
		{
			double correl = itsCorrel + epsilon; 
			if (gt(correl,1.)) correl=1; 
			else if (lt(correl,1.)) correl=-1; 
			item = new ICM_FlatCorrel(GetAsOfDate(),GetStructName(),
				GetIndex1(),GetIndex2(),correl); 
		}
		break ;
	default:
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_FlatCorrel::GenerateShiftBetas: not implemented "<<
			ICM_EnumsCnv::toString(typesensi)  ); 
	} ; 
	return item; 
}
//	-----------------------------------------------------------------------------
void 
ICM_FlatCorrel::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];
	int i = 0;

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

	fprintf(fOut, "\n ======> Flat Correlation Structure :\n\n");
	fprintf(fOut, "\t CorrelValue= %f\n",itsCorrel); 
	ICM_Correlation::View(id,fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}
