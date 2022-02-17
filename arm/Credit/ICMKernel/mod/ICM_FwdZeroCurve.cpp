
#include "icm_fwdzerocurve.h"



//	-------------------------------------------------------------------------------------------------------------
ICM_ZeroCurveForward::ICM_ZeroCurveForward(const ARM_ZeroCurve&ref,const ARM_Date&newAsOf)
: ARM_ZeroCurve(unconst(newAsOf),ref.GetCurrencyUnit()), itsCurve(ref)
{

	if (newAsOf<ref.GetAsOfDate())
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant create fwd curve at a date < asofdate");
	itsCorrectionFactor=unconst(ref).DiscountPrice(unconst(newAsOf)); 
}
//	-------------------------------------------------------------------------------------------------------------
ICM_ZeroCurveForward::ICM_ZeroCurveForward(const ICM_ZeroCurveForward&ref)
: ARM_ZeroCurve(ref.GetAsOfDate(),ref.GetCurrencyUnit()), itsCurve(ref.itsCurve) 
{}
//	-------------------------------------------------------------------------------------------------------------
ICM_ZeroCurveForward::~ICM_ZeroCurveForward()
{}
//	-------------------------------------------------------------------------------------------------------------
ICM_ZeroCurveForward& 
ICM_ZeroCurveForward::operator=(const ICM_ZeroCurveForward&ref)
{
	if (this!=&ref) 
	{
		this->~ICM_ZeroCurveForward(); 
		new(this)ICM_ZeroCurveForward(ref); 
	}
	return *this; 
}
//	-------------------------------------------------------------------------------------------------------------
//	virtual 
ARM_CLASS_NAME 
ICM_ZeroCurveForward::GetRootName(void)
{
	return unconst(itsCurve).GetRootName(); 
}
//	-------------------------------------------------------------------------------------------------------------
//	virtual 
void 
ICM_ZeroCurveForward::Copy(const ARM_Object* src)
{
	if (!src) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Cant copy from 0");
	*this=dynamic_cast<const ICM_ZeroCurveForward&>(*src); 
}
//	-------------------------------------------------------------------------------------------------------------
//	virtual 
ARM_Object* 
ICM_ZeroCurveForward::Clone(void)
{
	return new ICM_ZeroCurveForward(*this); 
}
//	-------------------------------------------------------------------------------------------------------------
//	virtual 
double ICM_ZeroCurveForward::DiscountFunction(double yearTerm)
{
	return unconst(itsCurve).DiscountPrice(yearTerm) / itsCorrectionFactor; 
}
 
