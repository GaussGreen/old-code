#include "ICMKernel\crv\ICM_interpoldefcrv.h"
#include "ARMKernel/ccy/currency.h"

//	-----------------------------------------------------------------------------------
ICM_InterpolDefCrv::ICM_InterpolDefCrv(const ICM_InterpolDefCrv&ref) : ICM_DefaultCurve(ref) 
{
	//itsInterpolationMethod=ref.itsInterpolationMethod; 
}
//	-----------------------------------------------------------------------------------
ICM_InterpolDefCrv::ICM_InterpolDefCrv(const ARM_Date& asOf,
						 const ARM_Vector& Dates,
						 const ARM_Vector& Inputs,
						 double Recovery,
						 ARM_ZeroCurve* zc, 
						 int intRule,
						 int adjStartDate,
						 qCDS_ADJ adj,
						 const std::string& ccy ,
						 const std::string& label ) 
						 // qINTERPOL_TYPE InterpolationMethod )
{
		
	Init();
	Set(asOf,Dates,Inputs, Recovery,zc,intRule,adjStartDate,adj,ccy,label);		
}

//	-----------------------------------------------------------------------------------
void 
ICM_InterpolDefCrv::Init(void)
{
	// SetName(ICM_INTERPOL_DEFCRV);
	// itsInterpolationMethod = qINTERPOL_LINEAR; 
}
//	-----------------------------------------------------------------------------------
// virtual 
/** void 
ICM_InterpolDefCrv::Copy(const ARM_Object* srcZc)
{
	ICM_DefaultCurve::Copy(srcZc);
	// BitwiseCopy(srcZc);
}**/ 
//	-----------------------------------------------------------------------------------
// virtual 
ARM_Object* 
ICM_InterpolDefCrv::Clone(void)
{
	return new ICM_InterpolDefCrv(*this); 
}

//	-----------------------------------------------------------------------------------
// virtual 
void 
ICM_InterpolDefCrv::Calibrate() 
{
	ICM_DefaultCurve::Calibrate(); 
}
//	-----------------------------------------------------------------------------------
// virtual 
// void 
// ICM_InterpolDefCrv::Calibrate_Stress_Test_Guess_Brent()
// {
// 	ICM_DefaultCurve::Calibrate_Stress_Test_Guess_Brent() ; 
// }


//	-----------------------------------------------------------------------------------
void 
ICM_InterpolDefCrv::Set (const ARM_Date& asOf,
						 const ARM_Vector& Dates,
						 const ARM_Vector& Inputs,
						 double Recovery,
						 ARM_ZeroCurve* zc,  	
						 int intRule,
						 int adjStartDate,
						 qCDS_ADJ adj,
						 const std::string& ccy ,
						 const std::string&  label )
						 // qINTERPOL_TYPE InterpolationMethod )
{
	ICM_DefaultCurve::Set(asOf,&Dates,&Inputs, Recovery,zc,intRule,adjStartDate,adj,true,ccy,label,NULL,qDEFCURVE_DICHO,"STD", ARM_Currency(ccy.c_str()).GetCreditStartDateLag());
}

//	-----------------------------------------------------------------------------------
void 
ICM_InterpolDefCrv::CptTermsSurvivalProba()
{
	// do nothing.
}

//	-----------------------------------------------------------------------------------
double 
ICM_InterpolDefCrv::SurvivalFunction(const double& yearterm) const 
{
	if (GetIsNameInDefault()) return 0.;

	bool equal = false;
	int i = 0;
	// Inf_Equal(yearterm,i,equal);
	i=position(yearterm); 
	if (eq(itsInterpolYF.Elt(i),yearterm)) 
	// if (equal) 
		return GetSurvivalProba(i);

	double SP_prev = GetSurvivalProba(i);
	double SP_after = GetSurvivalProba(i+1);
	double YF_prev = GetYearTerm(i); 
	double YF_after = GetYearTerm(i+1); 

	double Sproba = SP_prev+(SP_after-SP_prev)*(yearterm-YF_prev)/(YF_after-YF_prev);
				
	return Sproba;
}

//	-----------------------------------------------------------------------------------
void 
ICM_InterpolDefCrv::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[40];
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

   fprintf(fOut, "\n\n--------------------------------------------------------------------------------\n");
   fprintf(fOut, "---------------- Interpolated Default Curve ----------------------------------------\n");
   fprintf(fOut, "--------------------------------------------------------------------------------\n\n");

	ICM_DefaultCurve::View(id, fOut);

   	if ( ficOut == NULL )
    {
       fclose(fOut);
    }

}

// virtual 
/** ICM_DefaultCurve* 
ICM_InterpolDefCrv::xGenerateShiftCurve(double epsilon,qSENSITIVITY_TYPE mode ) const 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_InterpolDefCrv::GenerateShiftCurve: not implemented"); 
}**/ 
// virtual 
ICM_DefaultCurve* 
ICM_InterpolDefCrv::GenerateShiftCurve(const std::vector<std::string> & Terms, 
										  const ARM_Vector& epsilon ,
										  qSENSITIVITY_TYPE  ) const 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_InterpolDefCrv::GenerateShiftCurve: not implemented"); 
}
