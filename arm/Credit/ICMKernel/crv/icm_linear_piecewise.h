
#ifndef _ICM_LINEAR_PIECEWISE_CURVE_H
#define _ICM_LINEAR_PIECEWISE_CURVE_H


#include "ICMKernel\crv\icm_defaultcurve.h"


/*********************************************************************************/
/*! \class  ICM_Linear_Piecewise ICM_Linear_Piecewise.h "ICM_Linear_Piecewise.h"
 *  \author 
 *	\version 1.0
 *	\date   4 may 2004
 *	\file   ICM_Linear_Piecewise.h
 *	\brief  Default Curve */
/***********************************************************************************/


class ICM_Linear_Piecewise : public ICM_DefaultCurve
{
public:
		ICM_Linear_Piecewise(const ICM_Linear_Piecewise&ref); 
		ICM_Linear_Piecewise(const ARM_Date& asOf,
						 const vector<string>& terms,
						 ARM_Vector* rates,
						 double Recovery,
						 ARM_ZeroCurve* zc,
						int intRule,
						int adjStartDate,
						 qCDS_ADJ adj /*= qCredit_Default */,  	
						 const std::string& ccy /*= ARM_DEFAULT_CURRENCY*/,
						 const std::string& label /* = NULL */ ,
						 int Lag,
						 bool issummitcurve = true
						 );
		void Set (const ARM_Date& asOf,
					   const vector<string>& terms,
                       ARM_Vector* rates,
					   double Recovery,
					   ARM_ZeroCurve* zc,
						int intRule,
						int adjStartDate,
					   qCDS_ADJ adj /*= qCredit_Default*/,  	
                       const std::string& ccy /*= ARM_DEFAULT_CURRENCY*/,
					   const std::string& label /* = NULL */ ,
					   int Lag,
					   bool issummitcurve = true
					   );
        ICM_Linear_Piecewise(void) ; // { Init();}

		void Init();

private:
	virtual double SurvivalFunction(const double& yearterm) const;
public:
		virtual void CptTermsSurvivalProba(void);
		virtual void ResetLambda(const int& indice,const double& lambda);
		virtual void View(char* id = NULL, FILE* ficOut = NULL );
		virtual ICM_DefaultCurve* GenerateShiftCurve(const vector<string>& terms, 
												  const ARM_Vector& epsilon ,
												  qSENSITIVITY_TYPE mode ) const ;
/** JLA REMOVED
private:
	virtual ICM_DefaultCurve* xGenerateShiftCurve(double epsilon  ,
												  qSENSITIVITY_TYPE mode ) const ;
												  **/ 
public:

		virtual double DefProbInverse(const double& PriceIn) const ; 

		// void BitwiseCopy(const ARM_Object* src)
		// {}

        // virtual void Copy(const ARM_Object* srcZc) ;
		// {
		// ICM_DefaultCurve::Copy(srcZc);
		// 
		// BitwiseCopy(srcZc);
		// }
 
        virtual ARM_Object* Clone(void) ;
        // {
        //     ICM_Linear_Piecewise* theClone = new ICM_Linear_Piecewise();
		// 
        //      theClone->Copy(this);
		// 
        //     return(theClone);
        // }
 
	virtual void Calibrate() ; 
	// {
	// 	ICM_DefaultCurve::Calibrate(); 
	// }
	// virtual void Calibrate_Stress_Test_Guess_Brent() ;
	// {
	// 	ICM_DefaultCurve::Calibrate(); 
	// }
private:
	virtual ICM_DefaultCurve* GenDefCurve(const ARM_Vector& dates,
									  const ARM_Vector& streads,
									  const double& recovery,
									  const string& label,
									  const bool& isdates,
									  ARM_ZeroCurve* ircurve) const 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Linear_Piecewise::GenDefCurve:not implemented"); 
	}
private:
	ICM_Linear_Piecewise& operator=(const ICM_Linear_Piecewise&); //NA
};


#endif /*---- End of file ----*/

