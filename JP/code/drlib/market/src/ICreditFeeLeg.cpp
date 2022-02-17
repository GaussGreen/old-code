//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : ICreditFeeLeg.cpp
//
//   Description : Interface for credit fee leg
//
//   Author      : Antoine Gregoire
//
//   Date        : 3 Nov 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ICreditFeeLeg.hpp"

DRLIB_BEGIN_NAMESPACE

/** Loads ICreditFeeLeg interface */
static void loadICreditFeeLeg(CClassSP& clazz){
    REGISTER_INTERFACE(ICreditFeeLeg, clazz);
    EXTENDS(IGetMarket);
}

/** TYPE for ICreditFeeLeg */
CClassConstSP const ICreditFeeLeg::TYPE =
    CClass::registerInterfaceLoadMethod(
        "ICreditFeeLeg", typeid(ICreditFeeLeg), loadICreditFeeLeg);

bool ICreditFeeLegLoad() {
    return ICreditFeeLeg::TYPE != NULL;
}

/** Static utility method to be used for converting productTimeline to
	a map that facilitates quick conversion from any Date to a productTimeline index.
	if Date does not lie on the productTimeline, then it returns -1
*/
void
ICreditFeeLeg::buildDateToProductTimelineIndexMap(
	IntArray& dateToProductTimelineIndexMap,
	const DateTimeArray& productTimeline) //input
{
	dateToProductTimelineIndexMap.clear();
	int size = productTimeline.back().getDate() + 1;
	dateToProductTimelineIndexMap.resize(size, -1);

	for (int i=0; i<productTimeline.size(); ++i)
		dateToProductTimelineIndexMap[productTimeline[i].getDate()] = i;
}

/** Static utility method for creating scaling factors for the risky cashflows
	of the ICreditFeeLeg.
	scalingFactor = Coupon X YearFraction /lossConfigNotional * notional

*/
void
ICreditFeeLeg::buildScalingFactors(
	DoubleArray& scalingFactor, //output
	DoubleArray& historicalScalingFactor, //output
	const DateTime& valueDate, //input
	const AccrualPeriodArray& accrualPeriods, //input
	const DateTimeArray& payDates, //input
	const CouponNotionalTypesArray& couponNotionalTypes,	//input
	const DayCountConventionArray& dayCountConventions, //input
	const DoubleArray& coupons, //input
	const DoubleArray& notionals, //input
	double lossConfigNotional) //input
{
	int size = accrualPeriods.size();

	if (  ( size != dayCountConventions.size() ) ||
		( size != dayCountConventions.size() ) ||
		( size != coupons.size() ) ||
		( size != notionals.size() ) )
		throw ModelException("ICreditFeeLeg::buildScalingFactors",
			"Input arrays are off different sizes");

	scalingFactor.resize(size, 0);
	historicalScalingFactor.resize(size, 0);

	for (int i=0; i< size; ++i)
	{
		if (payDates[i] >= valueDate)
		{
			if ((couponNotionalTypes)[i] != AVERAGE)
			{
				double yearFraction =  dayCountConventions[i]->years(
											accrualPeriods[i]->startDate(),
											accrualPeriods[i]->endDate());
				
				scalingFactor[i] = (yearFraction *
									coupons[i] *
									notionals[i] /
									lossConfigNotional );
				
				continue; 
			}
			
			if (accrualPeriods[i]->startDate() < valueDate )
			{
				if (accrualPeriods[i]->endDate() <= valueDate )
				{
					double yearFraction =  dayCountConventions[i]->years(
												accrualPeriods[i]->startDate(),
												accrualPeriods[i]->endDate());

					historicalScalingFactor[i] = (yearFraction *
													coupons[i] *
													notionals[i] /
													lossConfigNotional );
				}
				else //we will have both scaling and historical scaling factors
				{
					double yearFraction = dayCountConventions[i]->years(
																		accrualPeriods[i]->startDate(),
																		valueDate);

					historicalScalingFactor[i] = (yearFraction *
													coupons[i] *
													notionals[i] /
													lossConfigNotional );

					yearFraction = dayCountConventions[i]->years(
																valueDate,
																accrualPeriods[i]->endDate());

					scalingFactor[i] = (yearFraction *
										coupons[i] *
										notionals[i] /
										lossConfigNotional );
				}
			}
			else
			{
				double yearFraction =  dayCountConventions[i]->years(
											accrualPeriods[i]->startDate(),
											accrualPeriods[i]->endDate());

				scalingFactor[i] = (yearFraction *
									coupons[i] *
									notionals[i] /
									lossConfigNotional );
			}
		}
	}
}


/*=============================================================================
 * IHasCreditFeeLeg
 *===========================================================================*/
IHasCreditFeeLeg::~IHasCreditFeeLeg(){}

void IHasCreditFeeLeg::load(CClassSP& clazz){
    REGISTER_INTERFACE(IHasCreditFeeLeg, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IHasCreditFeeLeg::TYPE =
CClass::registerInterfaceLoadMethod(
    "IHasCreditFeeLeg", typeid(IHasCreditFeeLeg), load);

/*=============================================================================
 * ICreditFeeLegGenerator
 *===========================================================================*/
ICreditFeeLegGenerator::~ICreditFeeLegGenerator(){}

void ICreditFeeLegGenerator::load(CClassSP& clazz){
    REGISTER_INTERFACE(ICreditFeeLegGenerator, clazz);
    EXTENDS(IObject);
}

CClassConstSP const ICreditFeeLegGenerator::TYPE =
CClass::registerInterfaceLoadMethod(
    "ICreditFeeLegGenerator", typeid(ICreditFeeLegGenerator), load);

DRLIB_END_NAMESPACE
