//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : KnockInFwd.cpp
//
//   Description : knock-in forward
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : October 28, 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KnockInFwd.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/CashSettleDate.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/InstrumentUtil.hpp"



DRLIB_BEGIN_NAMESPACE

const string KnockInFwd::SPREAD_NONE = "NONE";
const string KnockInFwd::SPREAD_INCREASING = "INCREASING";
const string KnockInFwd::SPREAD_DECREASING = "DECREASING";
const string KnockInFwd::SPREAD_SYMMETRIC = "SYMMETRIC";
const string KnockInFwd::SPREAD_SYMMETRICPCT = "SYMMETRICPCT";
const double KnockInFwd::MINIMUM_ALPHA = 0.001;


////////////////////
void KnockInFwd::Validate()
{
    static const string method = "KnockInFwd::Validate";
    int i;


    Generic1Factor::validate();
    
    // some schedules need to be the same length
    if (barrierLevels->size() != monitorDates->size()) {
        throw ModelException(method, "monitorDates and barrierLevels must be the same length");
    }
    if (numberOfShares->size() != monitorDates->size()) {
        throw ModelException(method, "monitorDates and numberOfShares must be the same length");
    }
    if (forwardStrikes->size() != monitorDates->size()) {
        throw ModelException(method, "monitorDates and forwardStrikes must be the same length");
    }
    if (usePutSchedule == true) {
        if (putStrikes->size() != monitorDates->size()) {
            throw ModelException(method, "monitorDates and putStrikes must be the same length");
        }
    }    

    if (monitorDates->size() < 1) {
        throw ModelException(method, "There must be at least one monitoring date");
    }

    if (forwardMatDates->size() < 1) {
        throw ModelException(method, "There must be at least one forward maturity date");
    }
    
    // historic monitor prices must be set.
    i = 0;
    while (i < monitorDates->size() && (*monitorDates)[i] < valueDate) {
        if (i >= histMonSamples->size() || (*histMonSamples)[i] <= 0.) {
            throw ModelException(method, "at least one historic monitoring sample has not been set");
        }
        i++;
    }
    
    // historic mat prices must be set. Only check for non-zero if in settlement period for them
    i = 0;
    while (i < fullForwardMatDates->size()) {
        if (valueDate >= (*fullForwardMatDates)[i] && valueDate < instSettle->settles((*fullForwardMatDates)[i], asset.get())) {
            if (i >= fullHistTradeSamples->size()) {
                throw ModelException(method, "at least one historic trade sample has not been set");
            }
            if ((*fullHistTradeSamples)[i] <= 0.) {
                throw ModelException(method, "historic trade samples must be set");
            }
        }
        i++;
    }
    
    
    // Monitor dates must be increasing
    DateTime::ensureIncreasing(*(monitorDates.get()),
                               "monitorDates",
                               true);

    // forwardMatDates must be increasing
    if (CString::equalsIgnoreCase(tradeDateStyle, "AllOnOrBefore")) {
        DateTime::ensureIncreasing(*(forwardMatDates.get()),
                                   "forwardMatDates",
                                   true);
    } 
    
    if (CString::equalsIgnoreCase(tradeDateStyle, "OneForOne")) {
        if (forwardMatDates->size() != monitorDates->size()) {
            throw ModelException(method, "monitorDates and forwardMatDates must be the same length");
        }
    }

    // forwardMatDates must be on or after corresponding monitor dates
    if (CString::equalsIgnoreCase(tradeDateStyle, "OneForOne")) {
        for (i=0; i<monitorDates->size(); i++) {
            if ((*monitorDates)[i] > (*forwardMatDates)[i]) {
                throw ModelException(method, "each monitorDate must be on or before respective forwardMatDate");
            }
        }
    } else {
        // just check the last date
        if ((*monitorDates)[monitorDates->size()-1] > (*forwardMatDates)[forwardMatDates->size()-1]) {
            throw ModelException(method, "last monitorDate must be on or before last forwardMatDate");
        }
    }
    
    // check spread types
    if (!(CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_NONE) ||
          CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_INCREASING) ||
          CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_DECREASING) ||
          CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_SYMMETRIC)  ||
          CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_SYMMETRICPCT))) {

        throw ModelException(method, "Unsupported spread type (" + spreadType +
             "). Valid spread types are NONE, INCREASING, DECREASING, SYMMETRIC, ant SYMMETRICPCT.");
    }

    if (!(CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_NONE)) && strikeSpread <= 0.) {
        throw ModelException(method, "strikeSpread must be > 0. for spread type " + spreadType + ".");
    }

    // check trade date style
    if (!(CString::equalsIgnoreCase(tradeDateStyle, "OneForOne") ||
          CString::equalsIgnoreCase(tradeDateStyle, "AllOnOrBefore"))) {

        throw ModelException(method,"Unsupported trade date style (" + tradeDateStyle +
                 "). Valid styles are OneForOne and AllOnOrBefore.");
    }

    // check upDown. This routine will fail if upDown is not valid  
    getUpDownBool(); 

    // check divTreatment
    if (!(CString::equalsIgnoreCase(divTreatment, "PassThroughPayDate") ||
          CString::equalsIgnoreCase(divTreatment, "PassThroughPayDateDollar") ||
          CString::equalsIgnoreCase(divTreatment, "None"))) {
        throw ModelException(method,"Unsupported divTreatment type (" + divTreatment + 
            "). Valid types are None, PassThroughPayDate, and PassThroughPayDateDollar.");
    } else if (CString::equalsIgnoreCase(divTreatment, "None") == false) {
        // divs must be dollar if passing them through    
        DividendListSP  divs = AssetUtil::getAllDivsBetweenDates(asset.get(), 
                                                                 (*monitorDates)[0], 
                                                                 (*forwardMatDates)[forwardMatDates->size()-1]);

        DateTime dummyDate;
        if (divs->hasYieldDividend() || divs->hasContinuousDividend(dummyDate)) {
            throw ModelException(method, "All divs must be dollar divs when they are passed through");
        }

        // no fwd starting for div pass through
//        if (startDate > valueDate) {
//            throw ModelException(method, "Forward starting not allowed if divs are passed through");
//        }

        // no div pass through for quantos or ccy struck
        if (!(ccyTreatment == CAsset::CCY_TREATMENT_NONE || ccyTreatment == CAsset::CCY_TREATMENT_VANILLA)) {
            throw ModelException(method, "No fancy ccy treatment if divs are passed through");
        }
    
    }

//    // physical or rolling cash settlement only
//    if (!(instSettle->isPhysical() || CashSettlePeriod::TYPE->isInstance(instSettle.get()))) {
//        throw ModelException(method, "Settlement must be physical or cash rolling");
//    }
    
    if (CashSettleDate::TYPE->isInstance(instSettle.get())) {
        DateTime cashSettleDate = instSettle->settles(valueDate, asset.get());
        if (cashSettleDate < (*forwardMatDates)[forwardMatDates->size()-1]) {
            throw ModelException(method, "Settlement cannot occur before last forward maturity date");
        }
    }

    // only rolling cash for quantos
    if (ccyTreatment == "P" && instSettle->isPhysical()) {
        throw ModelException(method, "Settlement can't be physical when ccy protected");
    }    
    
    // if fwd starting, can't be one contract
//    if (fwdStarting && oneContract) {
//        throw ModelException(method, "Can't be forward starting and one contract");
//    }

    if (fwdStarting) {
        if (startDate > (*monitorDates)[0]) {
            throw ModelException(method, "Start date can't be after 1st monitoring date");
        }

        if (valueDate >= startDate && initialSpot <= 0.) {
            throw ModelException(method, "Initial spot must be set if valuing on or after start date");
        }
    }
    
    return;
}

bool KnockInFwd::avoidVegaMatrix(const IModel* model) 
{


    // only do ccy struck if there's only one monitor date. 
    // Otherwise the strikes will vary by maturity even if the barrierLevels are flat.
    if (monitorDates->size() > 1 && ccyTreatment == CAsset::CCY_TREATMENT_STRUCK) {
        return true;
    }

    // return true if all of the options are historic. Note that there could still be 
    // vol sesitivity in the fwd for ccy prot. Ignore that for now.
    if ((*monitorDates)[monitorDates->size()-1] <= valueDate) {
        return true;
    }


    // ccy prot should be fine as the ccy prot strikes for the fwds will be automatically added
    
    // don't do for ccy protected either as the forwards also have a vol dependence
//    if (ccyTreatment == CAsset::CCY_TREATMENT_PROTECTED) {
//        return true;
//    }


    // the commented out code below checks that the barrier level and put schedules are flat.
    // reintroduce if the number of strikes for vega matrix becomes inconvenient.

    // do vega matrix if there's only one barrier level and at most one put strike
/*    int i;
    for (i=0; i<barrierLevels->size()-1; i++) {
        if (!Maths::equals((*barrierLevels)[i], (*barrierLevels)[i+1])) {
            return true;
        }
    }

    if (usePutSchedule == true) {
        for (i=0; i<putStrikes->size()-1; i++) {
            if (!Maths::equals((*putStrikes)[i], (*putStrikes)[i+1])) {
                return true;
            }
        }
    } */

    return false;
}


// Returns the strike for vega matrix. Throws if there's more than 1.
DoubleArraySP KnockInFwd::getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*       model)
{

    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));
    DoubleArraySP sensStrikesCall = DoubleArraySP(new DoubleArray(0));
    DoubleArraySP sensStrikesHigh = DoubleArraySP(new DoubleArray(0));
    DoubleArraySP sensStrikesLow = DoubleArraySP(new DoubleArray(0));
    DoubleArraySP sensStrikesPut = DoubleArraySP(new DoubleArray(0));

    double lowStrike;
    double highStrike;

    if (avoidVegaMatrix(model)) {
        throw ModelException("KnockInFwd::getSensitiveStrikes", 
                             "VEGA_MATRIX is not valid for this instrument");
    }

    // get start date for vol interpolation
    DateTime imntStartDate = fwdStarting ? startDate:valueDate;
    // get last exercise date in exercise schedule
    DateTime maturityDate = (*monitorDates)[monitorDates->size()-1];

    // create a vol request object to be passed on

    SensitiveStrikeDescriptor sensStrikeDesc;
    sensStrikeDesc.forwardOnly = false;

    int i;
    int j;
    for (i=0; i<monitorDates->size(); i++) {
        if ((*monitorDates)[i] > valueDate) {
            LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest((*barrierLevels)[i], 
                                                           imntStartDate, 
                                                           (*monitorDates)[i],
                                                           fwdStarting));
            asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                       sensStrikeDesc, sensStrikesCall);
            for (j=0; j<sensStrikesCall->size(); j++) {
                sensStrikes->push_back((*sensStrikesCall)[j]);
            }

            getDigStrikes(i, &lowStrike, &highStrike);
            volRequest->setStrike(highStrike);
            asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                       sensStrikeDesc, sensStrikesHigh);
            for (j=0; j<sensStrikesHigh->size(); j++) {
                sensStrikes->push_back((*sensStrikesHigh)[j]);
            }
            volRequest->setStrike(lowStrike);
            asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                       sensStrikeDesc, sensStrikesLow);
            for (j=0; j<sensStrikesLow->size(); j++) {
                sensStrikes->push_back((*sensStrikesLow)[j]);
            }

            if (usePutSchedule == true) {
                volRequest->setStrike((*putStrikes)[i]);
                asset->getSensitiveStrikes(volRequest.get(), outputName, 
                                           sensStrikeDesc, sensStrikesPut);
                for (j=0; j<sensStrikesPut->size(); j++) {
                    sensStrikes->push_back((*sensStrikesPut)[j]);
                }
            }
        }
    }

    return sensStrikes;
}

void KnockInFwd::validatePop2Object()
{
    static const string method("KnockInFwd::validatePop2Object");
    int i;
    int j;
    // one time initialization

    // I hate to throw errors in validatePopToObj but you can't initialize if these fail.
    // Could catch the error here and throw it in validate (like with the Bond object) but it doesn't
    // seem worth the effort for such a seldomly used product.

    // forwardMatDates must be increasing
    if (CString::equalsIgnoreCase(tradeDateStyle, "AllOnOrBefore")) {
        DateTime::ensureIncreasing(*(forwardMatDates.get()),
                                   "forwardMatDates",
                                   true);

    } 
    if (CString::equalsIgnoreCase(tradeDateStyle, "OneForOne")) {
        if (forwardMatDates->size() != monitorDates->size()) {
            throw ModelException(method, "monitorDates and forwardMatDates must be the same length");
        }
    }
    
    // resize the histMonSamples if they were passed in too small
    if (histMonSamples->size() != monitorDates->size()) {
        histMonSamples->resize(monitorDates->size());
    }

    fullForwardMatDates = DateTimeArraySP(new DateTimeArray());
    fullHistTradeSamples = DoubleArraySP(new DoubleArray());

    fullForwardMatDates->resize(monitorDates->size());
    fullHistTradeSamples->resize(monitorDates->size());

    if (CString::equalsIgnoreCase(tradeDateStyle, "OneForOne")) {
        for (i=0; i<monitorDates->size(); i++) {
            (*fullForwardMatDates)[i] = (*forwardMatDates)[i];
            if (i < histTradeSamples->size()) { // don't fill passed the end
                (*fullHistTradeSamples)[i] = (*histTradeSamples)[i];
            }
        }
    } else if (CString::equalsIgnoreCase(tradeDateStyle, "AllOnOrBefore")) {
        j = 0;
        for (i=0; i<monitorDates->size(); i++) {
            if ((*forwardMatDates)[j] < (*monitorDates)[i]) {
                j++;
                if (j >= forwardMatDates->size() || (*forwardMatDates)[j] < (*monitorDates)[i]) {
                    throw ModelException(method,
                             "Bad forwardMatDates. Can't skip a sample");
                }
            }

            (*fullForwardMatDates)[i] = (*forwardMatDates)[j];
            if (j < histTradeSamples->size()) { // don't fill passed the end
                (*fullHistTradeSamples)[i] = (*histTradeSamples)[j];
            }

        }
    } else {
        throw ModelException(method,
                 "Unsupported trade date style (" + tradeDateStyle +
                 "). Valid styles are OneForOne and AllOnOrBefore.");
    }
    
}

// up = true, down = false
bool KnockInFwd::getUpDownBool() const
{
    static const string method("KnockInFwd::getUpDownBool");
    bool upDownBool;

    if (CString::equalsIgnoreCase(upDown, "U") ||
        CString::equalsIgnoreCase(upDown, "up")) {
            upDownBool = true;
    } else if (CString::equalsIgnoreCase(upDown, "D") ||
        CString::equalsIgnoreCase(upDown, "Down")) {
            upDownBool = false;
    } else {
        throw ModelException(method, "Unsupported upDown type (" + upDown +
             "). Valid types are up and down.");
    }

    return upDownBool;
}

// 1 up, -1 down
double KnockInFwd::getUpDownDouble() const
{
    if (getUpDownBool() == true) {
        return 1.;
    } else {
        return -1.;
    }
}

// get the high and low strikes for the digital approximation. Also returns the strikeSpread for
// scaling the callSpread result. Not that for forward starting the strikes are percents but the
// strikeSpread is scaled by the spot at start so that it's always an absolute number.
void KnockInFwd::getDigStrikes(int index, double *lowStrike, double *highStrike) const
{
    static const string method("KnockInFwd::getDigStrikes");
    bool upDownBool = getUpDownBool();

    // just to be careful
    if (index < 0 || index >= monitorDates->size()) {
        throw ModelException(method, "index outside of range");
    }

    if (CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_NONE)) {
        *lowStrike = (*barrierLevels)[index]*(1. - KnockInFwd::MINIMUM_ALPHA);
        *highStrike = (*barrierLevels)[index]*(1. + KnockInFwd::MINIMUM_ALPHA);
    } else if (CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_INCREASING)) {
        if (upDownBool == true) {
            *lowStrike = (*barrierLevels)[index] - strikeSpread;
            *highStrike = (*barrierLevels)[index];
        } else {
            *lowStrike = (*barrierLevels)[index];
            *highStrike = (*barrierLevels)[index] + strikeSpread;
        }
    } else if (CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_DECREASING)) {
        if (upDownBool == true) {
            *lowStrike = (*barrierLevels)[index];
            *highStrike = (*barrierLevels)[index] + strikeSpread;
        } else {
            *lowStrike = (*barrierLevels)[index] - strikeSpread;
            *highStrike = (*barrierLevels)[index];
        }
    } else if (CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_SYMMETRIC)) {
        *lowStrike = (*barrierLevels)[index] - strikeSpread/2.;
        *highStrike = (*barrierLevels)[index] + strikeSpread/2.;
    } else if (CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_SYMMETRICPCT)) {
        *lowStrike = (*barrierLevels)[index]*(1. - strikeSpread/2.);
        *highStrike = (*barrierLevels)[index]*(1. + strikeSpread/2.);
    } else {
        throw ModelException(method, "Unsupported spread type (" + spreadType +
             "). Valid spread types are NONE, INCREASING, DECREASING, SYMMETRIC, and SYMMETRICPCT.");
    }

    return;
}

/** when to stop tweaking */
DateTime KnockInFwd::endDate(const Sensitivity* ) const {
    return (*forwardMatDates)[forwardMatDates->size()-1];
}

bool KnockInFwd::priceDeadInstrument(CControl* control, CResults * results) const
{   // to do
    return false;
}

/* true if we are at the end of a monitoring period */
bool KnockInFwd::isEndMonitoringPeriod(const DateTime& testDate) const
{
	int i = 0;
	bool isEndPeriod = false;	
	bool monitoringDateFound = false;
	bool nextMonitorDateFound = false;
	DateTime nextMonitoringDate;

	// checking whether testDate is a monitoring date
	while ( i<monitorDates->size() && !monitoringDateFound )
	{
		if( testDate.equals((*monitorDates)[i]) )
		{
			monitoringDateFound = true;
			if (i<monitorDates->size()-1)
			{ 
				nextMonitorDateFound = true;
				nextMonitoringDate = (*monitorDates)[i+1];
			}			
		}
		i++;
	}

	/* part commented out to return false if date is not a monitoring period
	if (!monitoringDateFound)
	{
		throw ModelException("KnockInFwd::isEndMonitoringPeriod", "date entered is not a monitoring date");
	}*/

	// now looking whether we are at a payment date
	if (monitoringDateFound)
	{
		if ( nextMonitorDateFound && ((*fullForwardMatDates)[i-1] < (*fullForwardMatDates)[i]) )
		{
			isEndPeriod = true;
		}
		if ( !nextMonitorDateFound)
		{
			isEndPeriod = true;
		}		
	}

	return isEndPeriod;
}

/** find pay date in PaymentDates */
const DateTime& KnockInFwd::getPayDate(const DateTime& today) const
{ 
	int i = 0;
	bool found = false;

	while (!found &&  i < monitorDates->size() )
	{
		if ( (*monitorDates)[i].equals(today) ) 
		{
			found = true;
		}

		i++;
	}

	if (found)
	{
		return (*fullForwardMatDates)[i - 1];
	}
	else
	{
		throw ModelException("KnockInFwd::getPayDate", "date supplied is not in monitoring dates!");
	}
}

// or the first after today            
DateTime KnockInFwd::getNextPaymentDate(const DateTime& today) const
{
	DateTime nextPayDate;

	bool found = false;
	int i = 0;
	while (!found && i < monitorDates->size())
	{
		if ( today <= (*monitorDates)[i] )
		{
			found = true;
			nextPayDate = (*fullForwardMatDates)[i];
		}
		i++;
	}
	if (!found) // in this case, we are after all monitoring dates so result here is last payment date
	{
		nextPayDate = (*forwardMatDates)[forwardMatDates->size()-1];
	}

	return nextPayDate;
}

// extract the first monitoring date corresponding to this settlement date
DateTime KnockInFwd::getFirstMonitoringDate(const DateTime& settleDate) const
{
	DateTime monitoringDate;

	bool foundMonitor = false;
	int i = 0;

	while (!foundMonitor && i < monitorDates->size())
	{
		if ( settleDate.equals((*fullForwardMatDates)[i]))
		{
			foundMonitor = true;
			monitoringDate = (*monitorDates)[i];
		}
		i++;
	}

	if (!foundMonitor)
	{
		throw ModelException("KnockInFwd::getFirstMonitoringDate related to ", " date " + settleDate.toString() + 
		" failed to find the first monitoring date");
	}

	return monitoringDate;
}


// gets relevant info for monitoring period containing evalDate
void KnockInFwd::getPeriodHistContractsAndStrikes(const DateTime& evalDate, 
                                      DoubleArray&                pdContracts, 
                                      DoubleArray&                pdStrikes, 
                                      double&                     putValue) const
{
    pdContracts.resize(0);
    pdStrikes.resize(0);
    putValue = 0.0;

	int i = 0;
	double pastSpotValue = 0.0;
	double upDownDouble = getUpDownDouble(); // up = 1, down = -1

	DateTime currentDate;
	DateTime settleDate;
	DateTime PayDate = getNextPaymentDate(evalDate);

	while ( i<monitorDates->size() && (*monitorDates)[i] <= evalDate) 
	{
		currentDate = (*monitorDates)[i];
		settleDate = getPayDate(currentDate);
		
		pastSpotValue = (*histMonSamples)[i];
		
		// in any case check this monitoring period's past values
		if (settleDate.equals(PayDate))
		{			
			double histMonSamp = (*histMonSamples)[i];
            //forwardStrike = (*instKIFwd->forwardStrikes)[i];
                
            // checking if barrier is crossed or not
			if (histMonSamp*upDownDouble > (*barrierLevels)[i]*upDownDouble)
			{
				pdContracts.push_back((*numberOfShares)[i]);
                pdStrikes.push_back((*forwardStrikes)[i]);
			}

            // always participate in the put
            if (usePutSchedule)
            {
                putValue += (*numberOfShares)[i] * Maths::max( ( (*putStrikes)[i] - histMonSamp)*upDownDouble, 0.);
            }

		}
		i++;	
	}	
}

// gets num contracts for period.  Required for computing num shares for divs.
double KnockInFwd::getPeriodHistContracts(const DateTime& evalDate) const
{
    DoubleArray pdContracts, pdStrikes;
    double putValue;
    
    getPeriodHistContractsAndStrikes(evalDate, pdContracts, pdStrikes, putValue);
    double nbContracts = 0.0;

    for (int i = 0; i < pdContracts.size(); i++)
    {
        nbContracts += pdContracts[i];
    }

    return nbContracts;
}

DateTimeArraySP KnockInFwd::getPaymentDates() const
{
	DateTimeArraySP payDates(new DateTimeArray(forwardMatDates->size()));
    
    for(int i=0; i<forwardMatDates->size(); i++)
	    (*payDates)[i] = instSettle->settles((*forwardMatDates)[i], asset.get());

	return payDates;
}

void KnockInFwd::getMATRIXinfo(CashFlowArray &cfl, PhysicalDeliveryArray &pda) const
{
	DateTime settleDate, pastSettleDate, monitorDate, divSettleDate, divDate;

	const DateTime endDate = (*forwardMatDates)[forwardMatDates->size() - 1];
	const DateTime matDate = (*monitorDates)[monitorDates->size() - 1];

    DateTimeArrayConstSP exDivDates;
    DoubleArrayConstSP exDivAmounts;

	cfl.resize(0);
	pda.resize(0);

    // return empty cfl if all monitor dates are in future
    if ( valueDate < monitorDates->front() ) return;

    // computing the scaling factor
    double fwdStrt = fwdStarting?asset->fwdValue(startDate):0.0;
    double scalingFactor = InstrumentUtil::scalePremium(oneContract,fwdStarting,notional,fwdStrt,initialSpot);

    // find the latest monitor date, can be pastDate
    int iEnd = monitorDates->size()-1, i, j; 
    while( iEnd >= 0 && (*monitorDates)[iEnd] > valueDate ) 
        iEnd--;
    if( iEnd < 0 )    // sanity check, should never happen
        throw ModelException("KnockInFwd::getMATRIXinfo", "Internal error.");
    
    // get all divs up front before we start looping through monitoring dates
    if( includeDivs()) 
    {
            DividendListSP exDivs = AssetUtil::getAllDivsBetweenDates(asset, startDate, valueDate);       
            exDivDates = exDivs->getExDivDates();
            exDivAmounts = exDivs->getDivAmounts();
    }

    j = 0;
    for(i=0; i<=iEnd; i++)
    {
        bool isEndPeriod = (i == (monitorDates->size()-1) || (*fullForwardMatDates)[i] != (*fullForwardMatDates)[i+1] );
        
        // for cash settle, add periodic payment
        // if end of monitor period, include current period accumulated contract.
        // otherwise, only include previous fully historical periods        
        if( isEndPeriod)
        {
            monitorDate = (*monitorDates)[i];
            settleDate = instSettle->settles((*fullForwardMatDates)[i], asset.get());

            // now retrieve any known num contracts, strikes, and put strikes for this period
            DoubleArray pdContracts(0), pdStrikes(0);
            double putValue = 0.0;
            getPeriodHistContractsAndStrikes(monitorDate, pdContracts, pdStrikes, putValue);

            // put value is always cash settled.
            double prodVal = putValue;
            
            // for cash settled instruments, add cash value of forward contract
            if (!instSettle->isPhysical())
            {
                prodVal += productValue((*fullHistTradeSamples)[i], pdContracts, pdStrikes);
            }
            
            if( !Maths::isZero(prodVal) )
                cfl.push_back(CashFlow(settleDate, scalingFactor*prodVal));

            if (instSettle->isPhysical())
            {
                double aggregateNumContracts = 0.0, aggregateStrike = 0.0;
                getDeliveryInfo(pdContracts, pdStrikes, aggregateNumContracts , aggregateStrike);
                PhysicalDelivery delivery(scalingFactor*aggregateNumContracts, 
                                          aggregateStrike, monitorDate, settleDate);        
                pda.push_back(delivery);
            }
        }
        
        // adding dividend payments as they become known
        if( includeDivs() && exDivDates->size() > 0 )
        {
            DateTime monitorStart = (*monitorDates)[i];
            DateTime monitorEnd = isEndPeriod?monitorStart.rollDate(1):(*monitorDates)[i+1];
            while( j < exDivDates->size() && 
                (*exDivDates)[j].getDate() >= monitorStart.getDate() &&
                (*exDivDates)[j].getDate() < monitorEnd.getDate() &&
                !Maths::isZero((*exDivAmounts)[j])) // if zero div do nothing
            {
                // contract does not include participation on ex div date
                DateTime divDate = DateTime( (*exDivDates)[j].getDate(), DateTime::START_OF_DAY_TIME);
                settleDate = instSettle->settles(divDate, asset.get());
                double divHistContracts = getPeriodHistContracts(divDate);
                double divAmount = divHistContracts*(*exDivAmounts)[j];
                if (!Maths::isZero(divAmount))
                    cfl.push_back(CashFlow(settleDate, scalingFactor*divAmount));
                j++;
            }
        }
    }

    return;
}

// sets aggregateNumContracts to the sum of pdContracts, aggregateStrike to the center of mass of strikes, weighted
// by pdContract.  
void KnockInFwd::getDeliveryInfo(const DoubleArray& pdContracts, 
                            const DoubleArray& pdStrikes, 
                            double &aggregateNumContracts, 
                            double &aggregateStrike) const
{
    if (pdContracts.size() != pdStrikes.size())
    {
        throw ModelException("KnockInFwd::getDeliveryInfo", "pdContracts and pdStrikes must have the same length!");
    }

    aggregateNumContracts = 0.0;
    aggregateStrike = 0.0;

    for (int i = 0; i < pdContracts.size(); i++)        
    {
        aggregateNumContracts += pdContracts[i];
        aggregateStrike += pdStrikes[i] * pdContracts[i];
    }
    
    if (!Maths::isZero(aggregateNumContracts))
    {
        aggregateStrike /= aggregateNumContracts;    
    }
    
}

// returns value of knocked in forward contracts for all known periods
double KnockInFwd::productValue(double spot,
                                const DoubleArray& pdContracts, 
                                const DoubleArray& pdStrikes) const
{
    double value = 0.0;
    for (int i = 0; i < pdContracts.size(); i++)
    {
        value += productValue(spot, pdContracts[i], pdStrikes[i]);
    }
    return value;
}

// returns value of knocked in forward contracts for one known period
double KnockInFwd::productValue(double spot,
                                double nbContracts,
                                double strike) const
{
    double forwardPrice = strike - spot;
    
    return nbContracts * forwardPrice;
}

// true if div treatment is "PassThroughOnPayDate" or "PassThroughOnPayDateDollar"
bool KnockInFwd::includeDivs() const
{
        return CString::equalsIgnoreCase(divTreatment, "PassThroughPayDate") ||
                    CString::equalsIgnoreCase(divTreatment, "PassThroughPayDateDollar");
}

//////////////////////////////////////////////////////////
/** product class for payoff */
//////////////////////////////////////////////////////////
class KnockInFwdClosedForm: public CClosedFormLN::IProduct{
private:
    const KnockInFwd*  instKIFwd; // a reference

public:
    KnockInFwdClosedForm(const KnockInFwd* instKIFwd): instKIFwd(instKIFwd){}

    void price(CClosedFormLN*   model,
               Control*         control, 
               CResults*        results) const;
protected:

};



// the main entry point for model
void KnockInFwdClosedForm::price(CClosedFormLN*   ,
                                 Control*        control, 
                                 CResults*       results) const
{
    static const string method = "KnockInFwdClosedForm::price";
    try {
        double futureSamplePrice = 0.;
        double digPrice;
        double callPrice;
        double putPrice;
        double discFactor;
        double fwdRatio;
        double barrier;
        double lowStrike;
        double highStrike;
        double numShares;
        double forwardStrike;
        double forwardPrice;
        double pastSamplePrice = 0.;
        double divPassThrough;
        double divCorrection;
        int i;
        int j;
        
        bool upDownBool = instKIFwd->getUpDownBool(); // up = true, down = false
        double upDownDouble = instKIFwd->getUpDownDouble(); // 1 = up, -1 = down
        
        HolidaySP hols(Holiday::noHolidays());
        InstrumentSettlementSP noSettle = InstrumentSettlementSP(new CashSettlePeriod(0, hols.get()));
        
        // Past knock ins
        for (i = 0; i<instKIFwd->monitorDates->size(); i++) {
            if ((*instKIFwd->monitorDates)[i] > instKIFwd->valueDate) {   
                 break;
            }
    
            double histMonSamp = (*instKIFwd->histMonSamples)[i];
            double histTradeSamp = (*instKIFwd->fullHistTradeSamples)[i];

            // we're before the forward maturity date
            if (instKIFwd->valueDate < (*instKIFwd->fullForwardMatDates)[i]) {
    
                numShares = (*instKIFwd->numberOfShares)[i];
                forwardStrike = (*instKIFwd->forwardStrikes)[i];
                
                if (histMonSamp*upDownDouble >= (*instKIFwd->barrierLevels)[i]*upDownDouble) {
                    forwardPrice = forwardStrike - instKIFwd->asset->fwdValue((*instKIFwd->fullForwardMatDates)[i]);
                } else {
                    forwardPrice = 0.;
                }

                if (instKIFwd->usePutSchedule == true) {
                    putPrice = Maths::max(((*instKIFwd->putStrikes)[i] - histMonSamp)*upDownDouble, 0.);
                } else {
                    putPrice = 0.;
                }


                discFactor = instKIFwd->instSettle->pv(
                    (*instKIFwd->fullForwardMatDates)[i],
                    instKIFwd->discount.get(), 
                    instKIFwd->asset.get());

                pastSamplePrice += discFactor*numShares*(forwardPrice + upDownDouble*putPrice);

            } else if (instKIFwd->valueDate < instKIFwd->instSettle->settles((*instKIFwd->fullForwardMatDates)[i], instKIFwd->asset.get())) { // <= ??
                // we're after the forward maturity date but we still need to include the value until the trade
                // settles.
                // But only if it cash settles. If it settles physically we need to drop it as a share trade is entered on
                // the trade date and not the settlement date as with cash.
                if (instKIFwd->instSettle->isPhysical() == false) 
                {
                    numShares = (*instKIFwd->numberOfShares)[i];
                    forwardStrike = (*instKIFwd->forwardStrikes)[i];
                
                    if (histMonSamp*upDownDouble >= (*instKIFwd->barrierLevels)[i]*upDownDouble) {
                        forwardPrice = forwardStrike - histTradeSamp;
                    } else {
                        forwardPrice = 0.;
                    }
                
                    if (instKIFwd->usePutSchedule == true) {
                        putPrice = Maths::max(((*instKIFwd->putStrikes)[i] - histMonSamp)*upDownDouble, 0.);
                    } else {
                        putPrice = 0.;
                    }
                
                    discFactor = instKIFwd->instSettle->pv(
                        (*instKIFwd->fullForwardMatDates)[i],
                        instKIFwd->discount.get(), 
                        instKIFwd->asset.get());
                    
                    pastSamplePrice += discFactor*numShares*(forwardPrice + upDownDouble*putPrice);
                }
            }   

            // dividend pass through. Do here as it can happen after settlement
            if (histMonSamp >= (*instKIFwd->barrierLevels)[i]) {
                if (CString::equalsIgnoreCase(instKIFwd->divTreatment, "PassThroughPayDate") ||
                    CString::equalsIgnoreCase(instKIFwd->divTreatment, "PassThroughPayDateDollar")) {
                    // any div that goes ex between monitor and forward mat date gets passed through on the
                    // dividends pay date (if the barrier was crossed, of course)
           
                    DividendListSP  divs = AssetUtil::getAllDivsBetweenDates(instKIFwd->asset.get(), 
                                                                             (*instKIFwd->monitorDates)[i], 
                                                                             (*instKIFwd->fullForwardMatDates)[i]);

                    DateTimeArrayConstSP payDates = divs->getPayDates();

                    DoubleArrayConstSP divAmounts = divs->getDivAmounts();

                    for (j=0; j<payDates->size(); j++) {
                        if (instKIFwd->valueDate < (*payDates)[j]) {
                            pastSamplePrice -= (*instKIFwd->numberOfShares)[i] * (*divAmounts)[j] * instKIFwd->discount->pv(instKIFwd->valueDate, (*payDates)[j]);
                        }
                    }
                }
            }
        }


        // Future knock ins
        for (i = instKIFwd->monitorDates->size()-1; i>=0; i--) {
            if ((*instKIFwd->monitorDates)[i] <= instKIFwd->valueDate) { 
                break;
            }
            
            barrier = (*instKIFwd->barrierLevels)[i];
            forwardStrike = (*instKIFwd->forwardStrikes)[i];
            numShares = (*instKIFwd->numberOfShares)[i];


            instKIFwd->getDigStrikes(i, &lowStrike, &highStrike);

            // the digital. Not a true dig as the strikes are percents for forward starting but that 
            // washes out as the barrier and fwd strike are in pcts below
            digPrice = CVanilla::priceSpread(
                instKIFwd->valueDate, // valueDate
                instKIFwd->startDate,  // startDate
                (*instKIFwd->monitorDates)[i], // matDate
                upDownBool, // isCall
                instKIFwd->fwdStarting, // fwdStarting,
                true, // oneContract,
                1, // notional,
                1, // initialSpot,
                lowStrike, // lowStrike,
                highStrike, // highStrike,
                noSettle.get(),
                instKIFwd->asset.get(),
                instKIFwd->discount.get()) / (highStrike - lowStrike);
            
            divPassThrough = 0;
            if (CString::equalsIgnoreCase(instKIFwd->divTreatment, "PassThroughPayDate") ||
                CString::equalsIgnoreCase(instKIFwd->divTreatment, "PassThroughPayDateDollar")) {
                // any div that goes ex between monitor and forward mat date gets passed through on the
                // dividends pay date (if the barrier was crossed, of course)
               
                DividendListSP  divs = AssetUtil::getAllDivsBetweenDates(instKIFwd->asset.get(), 
                                                                         (*instKIFwd->monitorDates)[i], 
                                                                         (*instKIFwd->fullForwardMatDates)[i]);

                DateTimeArrayConstSP payDates = divs->getPayDates();

                DoubleArrayConstSP divAmounts = divs->getDivAmounts();

                for (j=0; j<payDates->size(); j++) {
                    divPassThrough += (*divAmounts)[j] * instKIFwd->discount->pv((*instKIFwd->monitorDates)[i], (*payDates)[j]);
                }
            }

            callPrice = CVanilla::priceBS(
                instKIFwd->valueDate, // valueDate
                instKIFwd->startDate,  // startDate
                (*instKIFwd->monitorDates)[i], // matDate
                upDownBool,  // isCall
                instKIFwd->fwdStarting, // fwdStarting,
                true,  // oneContract,
                1,     // notional,
                1,     // initialSpot,
                barrier, // strike,
                noSettle.get(),
                instKIFwd->asset.get(),
                instKIFwd->discount.get());

            if (instKIFwd->usePutSchedule == true) {
                putPrice = CVanilla::priceBS(
                    instKIFwd->valueDate, // valueDate
                    instKIFwd->startDate,  // startDate
                    (*instKIFwd->monitorDates)[i], // matDate
                    !upDownBool,  // isCall
                    instKIFwd->fwdStarting, // fwdStarting,
                    true,  // oneContract,
                    1,     // notional,
                    1,     // initialSpot,
                    (*instKIFwd->putStrikes)[i], // strike,
                    noSettle.get(),
                    instKIFwd->asset.get(),
                    instKIFwd->discount.get());
            } else {
                putPrice = 0.;
            }

            // calculate the discount between the monitoring date and settlement of the futures
            discFactor = instKIFwd->instSettle->pv(
                (*instKIFwd->fullForwardMatDates)[i],
                instKIFwd->discount.get(), 
                instKIFwd->asset.get()) / instKIFwd->discount->pv((*instKIFwd->monitorDates)[i]);

            fwdRatio = instKIFwd->asset->fwdValue((*instKIFwd->fullForwardMatDates)[i]) / instKIFwd->asset->fwdValue((*instKIFwd->monitorDates)[i]);

            if (CString::equalsIgnoreCase(instKIFwd->divTreatment, "PassThroughPayDate")) {
                divCorrection = numShares*divPassThrough/instKIFwd->asset->fwdValue((*instKIFwd->monitorDates)[i])*(barrier*digPrice+callPrice);
            } else if (CString::equalsIgnoreCase(instKIFwd->divTreatment, "PassThroughPayDateDollar")){
                divCorrection = numShares*divPassThrough*(digPrice);
                if (instKIFwd->fwdStarting == true) { // digPrice is off by this factor
                    divCorrection /= instKIFwd->asset->fwdValue(instKIFwd->startDate);
                }
            } else {
                divCorrection = 0;
            }
            

            futureSamplePrice += discFactor*numShares*((forwardStrike - fwdRatio*barrier)*digPrice - upDownDouble*fwdRatio*callPrice + upDownDouble*putPrice) - divCorrection;
       
        }
        
        double returnPrice = futureSamplePrice+pastSamplePrice;

        if (instKIFwd->oneContract == false) {
            if (instKIFwd->fwdStarting == true) {
                returnPrice *= instKIFwd->notional / instKIFwd->asset->fwdValue(instKIFwd->startDate);
            } else {
                returnPrice *= instKIFwd->notional / instKIFwd->initialSpot;
            }
        }

        results->storePrice(returnPrice, instKIFwd->discount->getCcy());

        // add any requested outputs such as FWD_AT_MAT, DELAY_PRICE
        instKIFwd->addRequests(control,
                               results,
                               futureSamplePrice+pastSamplePrice,
                               instKIFwd->endDate(NULL));

		// record KNOWN_CASHFLOWS and PHYSICAL_DELIVERY
		if ( control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS) || instKIFwd->instSettle->isPhysical()) {
			CashFlowArray allCFs(0);
            PhysicalDeliveryArray pda(0);
            
            // now populate the MATRIX info
            instKIFwd->getMATRIXinfo(allCFs, pda);
            
            // KNOWN_CASHFLOWS 
            OutputRequestUtil::recordKnownCashflows(control,
                results,
                instKIFwd->discount->getCcy(),
                &allCFs); 
            
            // PHYSICAL_DELIVERY
            if (instKIFwd->instSettle->isPhysical()) 
            {
                if ( control->requestsOutput(OutputRequest::PHYSICAL_DELIVERY) ) 
                {
                    PhysicalDelivery::recordPhysicalDelivery(control,
                        results,
                        instKIFwd->asset->getTrueName(),
                        &pda);
                }                
            }

		}
		
		// PAYMENT_DATES
		if ( control->requestsOutput(OutputRequest::PAYMENT_DATES) ) {
			DateTimeArraySP dates = instKIFwd->getPaymentDates();
			
			OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
		}

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* KnockInFwd::createProduct(
    CClosedFormLN* ) const{

    return new KnockInFwdClosedForm(this);
}

/** Rolls the value date for theta */
bool KnockInFwd::sensShift(Theta* shift)
{    
    int i;
    DateTime newDate = shift->rollDate(valueDate);

    if (fwdStarting && newDate.isGreaterOrEqual(startDate) &&
        startDate.isGreaterOrEqual(valueDate))
    {
        fwdStarting = false;
        initialSpot = asset->getThetaSpotOnDate(shift, startDate);

        // scale the schedules
        // double numOptions = oneContract? 1. : notional/initialSpot;

        // don't scale if it's always pct
        if (CString::equalsIgnoreCase(spreadType, KnockInFwd::SPREAD_SYMMETRICPCT) == false) {
            strikeSpread *= initialSpot;
        }
        
        for (i=0; i<monitorDates->size(); i++) {
            (*barrierLevels)[i] = (*barrierLevels)[i]*initialSpot; 
  //          (*numberOfShares)[i] = (*numberOfShares)[i]*numOptions; 
            (*forwardStrikes)[i] = (*forwardStrikes)[i]*initialSpot; 
            if (usePutSchedule == true) {
                (*putStrikes)[i] = (*putStrikes)[i]*initialSpot; 
            }
        }
    }


    if (valueDate == newDate) {
        // only set stuff if it was not set before
        for (i=0; i<monitorDates->size(); i++) {
            if (newDate == (*monitorDates)[i] && Maths::isZero((*histMonSamples)[i])) {
                (*histMonSamples)[i] = asset->getThetaSpotOnDate(shift, (*monitorDates)[i]);
            }

            if (newDate == (*fullForwardMatDates)[i] && Maths::isZero((*fullHistTradeSamples)[i])) {
                (*fullHistTradeSamples)[i] = asset->getThetaSpotOnDate(shift, (*fullForwardMatDates)[i]);
            }
        }

    } else {

        for (i=0; i<monitorDates->size(); i++) {
            if (newDate.isGreaterOrEqual((*monitorDates)[i]) && (*monitorDates)[i].isGreaterOrEqual(valueDate)) {
                (*histMonSamples)[i] = asset->getThetaSpotOnDate(shift, (*monitorDates)[i]);
            }

            if (newDate.isGreaterOrEqual((*fullForwardMatDates)[i]) && (*fullForwardMatDates)[i].isGreaterOrEqual(valueDate)) {
                (*fullHistTradeSamples)[i] = asset->getThetaSpotOnDate(shift, (*fullForwardMatDates)[i]);
            }
        }
    }

    // roll today 
    valueDate = newDate;
    
    return true;
};

// for reflection
KnockInFwd::KnockInFwd(): Generic1Factor(TYPE){
}
    

class KnockInFwdHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(KnockInFwd, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        EMPTY_SHELL_METHOD(defaultKnockInFwd);
        FIELD(monitorDates, "one date for eack ki fwd");
        FIELD(barrierLevels, "knock-in level");
        FIELD(numberOfShares, "number of shares underlying ki");
        FIELD(forwardMatDates, "forward trade dates");
        FIELD(forwardStrikes, "strikes of the forwards");
        FIELD(histMonSamples, "level of u/l on hist sample dates");
        FIELD(histTradeSamples, "level of u/l on hist forward maturity dates");
        FIELD(spreadType, "increasing, decreasing, or none");
        FIELD(strikeSpread, "spread for barrier call spread approx");
        FIELD(tradeDateStyle, "OneForOne, AllOnOrBefore, AllOnOrBeforeNoSettlement");
        FIELD(upDown, "u/d up and in or down and out");
        FIELD(usePutSchedule, "include the puts");
        FIELD(putStrikes, "put strikes");
        FIELD_MAKE_OPTIONAL(putStrikes);
        FIELD(divTreatment, "none, passThroughEx, passThroughPay, passThroughAtSettle");
        FIELD(fullForwardMatDates, "private");
        FIELD_MAKE_TRANSIENT(fullForwardMatDates);
        FIELD(fullHistTradeSamples, "private");
        FIELD_MAKE_TRANSIENT(fullHistTradeSamples);    
    }

    static IObject* defaultKnockInFwd(){
        return new KnockInFwd();
    }
};

CClassConstSP const KnockInFwd::TYPE = CClass::registerClassLoadMethod(
    "KnockInFwd", typeid(KnockInFwd), KnockInFwdHelper::load);
   
bool KnockInFwdLoad()
{
    return (KnockInFwd::TYPE != 0);
}




DRLIB_END_NAMESPACE

