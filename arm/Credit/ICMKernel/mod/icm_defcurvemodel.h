
#ifndef _DEFAULT_CURVE_MODEL_H
#define _DEFAULT_CURVE_MODEL_H

#include "ARMKernel\mod\model.h"
#include "ARMKernel\crv\volcurv.h"
#include "ICMKernel\crv\icm_defaultcurve.h"

/*********************************************************************************/
/*! \class  ICM_DefaultCurveModel icm_defcurvemodel.h "icm_defcurvemodel.h"
 *  \author Damien Pouponneau 
 *	\version 1.0
 *	\date   May 2004
 *	\brief  Default Model */
/***********************************************************************************/

class ICM_DefaultCurveModel : public ARM_Model  
{        
private :
	ARM_ZeroCurve*		itsZeroCurve;		//Zero Curve
	const ICM_DefaultCurve*	itsDefaultCurve;	//Default Curve
	ARM_VolCurve*		itsVolCurve;		//Volatility Curve for options
	bool				itsFlgClone;		//Clone Object or Not
public:

    ICM_DefaultCurveModel(void)
    {
        Init();
    }
	ICM_DefaultCurveModel(const ICM_DefaultCurveModel&ref);  
	ICM_DefaultCurveModel(const ICM_DefaultCurve* DefProbCurve,
						  ARM_ZeroCurve* ZeroCurve,
						  ARM_VolCurve*	VolCurve ,
						  bool FlgClone ) ;

	ICM_DefaultCurveModel(const ICM_DefaultCurve* DefProbCurve,
						  ARM_ZeroCurve* ZeroCurve,
						  const ARM_VolCurve* VolCurve = NULL ) ;
	void Set(const ICM_DefaultCurve* DefProbCurve, ARM_ZeroCurve* ZeroCurve, ARM_VolCurve* VolCurve = NULL,bool FlgClone = true);
    virtual ~ICM_DefaultCurveModel(void) ;
	// virtual void Copy(const ARM_Object* src) ;
	virtual ARM_Object* Clone(void) ;

	const ICM_DefaultCurve* GetDefaultCurve() const { return itsDefaultCurve;};
	void SetDefaultCurve(const ICM_DefaultCurve* notcrv) ;

	virtual ARM_ZeroCurve* GetZeroCurve() const { return itsZeroCurve;};
	virtual void SetZeroCurve(ARM_ZeroCurve* crv) ;

	ARM_VolCurve* GetVolCurve() const { return itsVolCurve;};
	void SetVolCurve(ARM_VolCurve* crv) ;
	virtual void View(char* id, FILE* ficOut);
	bool GetFlgClone() {return itsFlgClone;};
protected:
	
	// void SetFlgClone(bool value) {itsFlgClone = value;};
public:
	ICM_DefaultCurveModel* GenerateShiftModel(qSENSITIVITY_TYPE typesensi,
		const std::string& , 
		double epsilon );
private:
	// void BitwiseCopy(const ARM_Object* src) ;
	void Init(void) ; // do not handle memory cleanup: private
private:
	
	ICM_DefaultCurveModel& operator=(const ICM_DefaultCurveModel&); // NA 
};

#endif
