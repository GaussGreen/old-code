
#ifndef _ICM_INTERPOL_DEFAULT_CURVE_H
#define _ICM_INTERPOL_DEFAULT_CURVE_H

#include "ICMKernel\crv\icm_defaultcurve.h"

/*********************************************************************************/
/*! \class  ICM_InterpolDefCrv ICM_InterpolDefCrv.h "ICM_InterpolDefCrv.h"
 *  \author Valentin R.
 *	\version 1.0
 *	\date   12 april 2005
 *	\file   ICM_InterpolDefCrv.h
 *	\brief  Default Curve */
/***********************************************************************************/

class ICM_InterpolDefCrv : public ICM_DefaultCurve
{

public:
	ICM_InterpolDefCrv(const ICM_InterpolDefCrv&ref); 
	ICM_InterpolDefCrv(const ARM_Date& asOf,
					 const ARM_Vector& Dates,
					 const ARM_Vector& Inputs,
					 double Recovery,
					 ARM_ZeroCurve* zc,  	
					 int intRule,
					 int adjStartDate,
					 qCDS_ADJ adj,
					 const std::string& ccy  ,
					 const std::string& label /*= NULL */ ); 
					 // qINTERPOL_TYPE InterpolationMethod = qINTERPOL_LINEAR);


	void Set (const ARM_Date& asOf,
					 const ARM_Vector& Dates,
					 const ARM_Vector& Inputs,
					 double  Recovery,
					 ARM_ZeroCurve* zc,  	
					 int intRule,
					 int adjStartDate,
					 qCDS_ADJ adj,
					 const std::string& ccy  ,
					 const std::string& label /* = NULL */ ); 
					 // qINTERPOL_TYPE InterpolationMethod = qINTERPOL_LINEAR);

    ICM_InterpolDefCrv(void) { Init();}


	// void BitwiseCopy(const ARM_Object* src)
	// {}

    // virtual void Copy(const ARM_Object* srcZc) ;
	/** {
	ICM_DefaultCurve::Copy(srcZc);

	BitwiseCopy(srcZc);
	}**/ 

    virtual ARM_Object* Clone(void) ;
    /** {
        ICM_InterpolDefCrv* theClone = new ICM_InterpolDefCrv();

         theClone->Copy(this);

        return(theClone);
    } */ 


	virtual void View(char* id = NULL, FILE* ficOut = NULL );

	void Init();
	
	virtual void CptTermsSurvivalProba();

private:
	virtual double SurvivalFunction(const double& yearterm) const ;
	virtual void Calibrate() ; 
	/*{
		ICM_DefaultCurve::Calibrate(); 
	}*/
	// virtual void Calibrate_Stress_Test_Guess_Brent();
	/*{
		ICM_DefaultCurve::Calibrate_Stress_Test_Guess_Brent() ; 
	}*/
protected:
	virtual ICM_DefaultCurve* GenerateShiftCurve(const std::vector<std::string> & Terms, 
											  const ARM_Vector& epsilon ,
											  qSENSITIVITY_TYPE  ) const  ; 
/** JLA Removed 
private:
	virtual ICM_DefaultCurve* xGenerateShiftCurve(double epsilon  ,
											  qSENSITIVITY_TYPE mode ) const  ; 
											  **/ 
	
private:
	virtual double DefProbInverse(const double& PriceIn) const
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_InterpolDefCrv::DefProbInverse:not implemented"); 
	}
	virtual ICM_DefaultCurve* GenDefCurve(const ARM_Vector& dates,
										  const ARM_Vector& streads,
										  const double& recovery,
										  const string& label,
										  const bool& isdates,
										  ARM_ZeroCurve* ircurve) const
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_InterpolDefCrv::GenDefCurve:not implemented"); 
	}
private:
	ICM_InterpolDefCrv& operator=(const ICM_InterpolDefCrv&); //N 
};


#endif /*---- End of file ----*/
