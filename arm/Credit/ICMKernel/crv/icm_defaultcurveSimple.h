#ifndef _ICM_DEFCURVE_SIMPLE_H
#define _ICM_DEFCURVE_SIMPLE_H

#include "ICMKernel\crv\icm_defaultcurve.h"

/*********************************************************************************/
/*! \class  ICM_DefCurveSimple icm_defaultcurveSimple.h "icm_defaultcurveSimple.h"
 *  \author 
 *	\version 1.0
 *	\date   4 may 2004
 *	\file   icm_defaultcurveSimple.h
 *	\brief  Default Curve simplified */
/***********************************************************************************/

class ICM_DefCurveSimple : public ICM_DefaultCurve
{
private:

	double itsSpread;
	double itsRecovery;

public:
	ICM_DefCurveSimple()
	{Init();}

	ICM_DefCurveSimple(const double spread,
					 const double recovery,
					 string label)
	{
		Init();

		itsSpread=spread;
		itsRecovery=recovery;
		SetLabel(label);
	}

private:
	virtual double SurvivalFunction(const double& yearterm) const 
	{
		double lambda = itsSpread*1.e-4/(1.-itsRecovery);
		return exp(-lambda*yearterm);
	}
public:
	void Init()
	{
		SetName(ICM_DEFCURVE_QUICK);

		itsSpread=-1.;
		itsRecovery=-1;
	}

	void BitwiseCopy(const ARM_Object* src)
	{
		ICM_DefCurveSimple* dc = (ICM_DefCurveSimple*)src;
		itsSpread=dc->itsSpread;
		itsRecovery=dc->itsRecovery;
	}

    virtual void Copy(const ARM_Object* srcZc)
	{
	ICM_DefaultCurve::Copy(srcZc);
	BitwiseCopy(srcZc);
	}

    virtual ARM_Object* Clone(void)
    {
        ICM_DefCurveSimple* theClone = new ICM_DefCurveSimple();
        theClone->Copy(this);
        return(theClone);
    }

	virtual void Calibrate() { ICM_DefaultCurve::Calibrate(); }
	// virtual void Calibrate_Stress_Test_Guess_Brent() { ICM_DefaultCurve::Calibrate_Stress_Test_Guess_Brent(); }

private:
	virtual ICM_DefaultCurve* GenerateShiftCurve(const std::vector<std::string> & Terms, 
									  const ARM_Vector& epsilon ,
									  qSENSITIVITY_TYPE  ) const
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefCurveSimple::GenerateShiftCurve:not implemented"); 
	}
/** JLA Removed
private:virtual ICM_DefaultCurve* xGenerateShiftCurve(double epsilon  ,
											  qSENSITIVITY_TYPE mode ) const
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefCurveSimple::GenerateShiftCurve:not implemented"); 
	}
	**/ 
public:
	virtual ICM_DefaultCurve* GenDefCurve(const ARM_Vector& dates,
										  const ARM_Vector& streads,
										  const double& recovery,
										  const string& label,
										  const bool& isdates,
										  ARM_ZeroCurve* ircurve) const 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefCurveSimple::GenDefCurve:not implemented"); 
	}
	virtual double DefProbInverse(const double& PriceIn) const 
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_DefCurveSimple::DefProbInverse:not implemented"); 
	}

private:
	ICM_DefCurveSimple(const ICM_DefCurveSimple&); //NA 
	ICM_DefCurveSimple& operator=(const ICM_DefCurveSimple&); //NA 
};


#endif /*---- End of file ----*/
