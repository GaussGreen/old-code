
#ifndef _ICM_Zero_Curve_Fwd_
#define _ICM_Zero_Curve_Fwd_

#include "ICMKernel/util/icm_macro.h"
#include "ARMKernel/crv/zerocurv.h"


//	-------------------------------------------------------------------------------------------------------------
class	ICM_ZeroCurveForward : public ARM_ZeroCurve
{
public:
	ICM_ZeroCurveForward(const ARM_ZeroCurve&ref,const ARM_Date&newAsOf); 
	ICM_ZeroCurveForward(const ICM_ZeroCurveForward&ref); 
	ICM_ZeroCurveForward& operator=(const ICM_ZeroCurveForward&ref); 
	virtual ~ICM_ZeroCurveForward() ; 
	
	//	--	Overloaded. 
	virtual double DiscountFunction(double yearTerm) ; 
	virtual ARM_CLASS_NAME GetRootName(void) ;

	//	--	 Not Implemented 
	virtual double D1DiscountFunction(double yearTerm) { ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented"); }
	virtual double D2DiscountFunction(double yearTerm) { ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented"); }
	virtual ARM_Object& operator = (const ARM_Object& obj) { ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented"); }
	virtual void View(char* id = NULL, FILE* ficOut = NULL) { ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented"); }
	virtual void Print(void) { ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented"); }
	
	//	--	ARM_Object 
	virtual void Copy(const ARM_Object* src) ;
	virtual ARM_Object* Clone(void) ; 
private:
	const ARM_ZeroCurve& itsCurve; 
	double itsCorrectionFactor ; 
} ;


#endif // _ICM_Zero_Curve_Fwd