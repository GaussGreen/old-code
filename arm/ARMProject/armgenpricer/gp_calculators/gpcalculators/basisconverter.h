
#ifndef _INGPCALCULATORS_BASISCONVERTER_H
#define _INGPCALCULATORS_BASISCONVERTER_H

#include "gpbase/rootobject.h"
#include "gpbase/port.h"
#include "gpcalib/typedef.h"
#include "gpbase/assignop.h"
#include "gpbase/curve.h"


#include "gpbase/datestripcombiner.h"
#include <ccy/currency.h>
//#include <inst/forex.h>

class ARM_Currency;
class ARM_SwapLeg;
class ARM_Forex;

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
/// \class ARM_BasisConverter 
/// this class converts a  foriegn funding leg with margin  
///	from the first currency ( foreign) using the ZC Curve to evaluate forward rate and 
/// spreaded curve to  discount to the second curreny ( domestic). 
/// This conversion will be done flow by flow. 
/////////////////////////////////////////////////////////////////

class ARM_BasisConverter : public ARM_RootObject
{
public:
	/// constructor, copy constructor, assignment operator and destructor
	ARM_BasisConverter( const ARM_Currency& domesticCcy,
		const ARM_Currency& foreignCcy,
		const ARM_DateStripPtr&	domdatestrip,
		const ARM_DateStripPtr&	fordatestrip,
		const ARM_DateStripPtr&	funddatestrip,
		ARM_ZeroCurve* domesticZeroCurve,
		ARM_ZeroCurve* foreignZeroCurve,
		ARM_ZeroCurve* domesticDiscountZeroCurve,
		ARM_ZeroCurve* foreignDiscountZeroCurve,
		const ARM_Forex&        forex, 
		const std::vector<double>&	domesticNotional,
		const std::vector<double>&	foreignNotional,
		const std::vector<double>&	foreignSpread);

	ARM_BasisConverter( const ARM_BasisConverter& rhs );
	ASSIGN_OPERATOR(ARM_BasisConverter)
	virtual ~ARM_BasisConverter() {};

	virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_BasisConverter"; }
	//virtual ARM_CLASS_NAME GetRootName() { return ARM_BASISCONVETER; }
	virtual ARM_Object* Clone() const {return new ARM_BasisConverter(*this);};

	/// Function to calculate the domestic margin.
	std::vector<double>  ComputeDomMargin( ARM_SwapLeg* fundingLeg = NULL);

private:

	ARM_Date            itsAsOfDate; 
	ARM_DateStripPtr	itsDomesticDatestrip;
	ARM_DateStripPtr	itsForeignDatestrip;
	ARM_DateStripPtr	itsFundingDatestrip;
	std::vector<double>		itsDomesticNotional;
    std::vector<double>		itsForeignNotional;
	std::vector<double>		itsForeignSpread;

	/// the discount curve, Please don't delete the following attributs
	/// They are shared. I didn't use Ptr to avoid an unwanted clone.
	ARM_ZeroCurve* itsDomesticZeroCurve;
	ARM_ZeroCurve* itsForeignZeroCurve;
	ARM_ZeroCurve* itsDomesticDiscountZeroCurve;
	ARM_ZeroCurve* itsForeignDiscountZeroCurve;

	/// Forex
	ARM_Forex itsForex;

	ARM_Currency   itsDomesticCcy;
	ARM_Currency   itsForeignCcy;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

