//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : BondParams.cpp
//
//   Description : Full Bond implementation from parameters
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : December 12, 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BondParams.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE
/** Full Bond implementation from parameters */

BondParams::~BondParams() {
    // empty
}

BondParams::BondParams(): Bond(TYPE){
    isValid = false;
    errorMessage = "bond not initialized";

    redemptionPct = 1.;
    dayCountConvString = "30/360";
    oddLastShort = false;
    endOfMonthAdj = true;
    eomIgnoreLeapYear = false;
    badDayConvString = "N";
    exDivDays = 0;
    exDivRule = "N";
    taxRate = 0.;
    
    isAnOID = false;
    yieldForOID = 0.;
 
    oddFirst = false;
    oddLast = false;
    oddLastLong = false;
    oddFirstLong = false;
}


/** Change the settlement of the bond if applicable */
void BondParams::setSettlement(SettlementSP newSettle) {
    settle = SettlementSP(newSettle.get());
}

CashFlowArraySP BondParams::getCashFlows() const
{
    static const string method = "BondParams::getCashFlows";
    try {  
        throwIfNotValid();

        return CashFlowArraySP(cashFlows.clone());
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

CashFlowArrayConstSP BondParams::getCashFlowRef() const
{
    static const string method = "BondParams::getCashFlowsRef";
    try {
        throwIfNotValid();

        return cashFlows;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Pull out the component assets & correlations from the market data */
void BondParams::getMarket(const IModel* model, const MarketData* market){
    try {
        hols.getData(model, market);
        settle->getMarket(model, market);

        initialize();
        isValid = true;
    } 
    catch (exception& e){
        throw ModelException(e, "BondParams::getMarket", "Failed for BondParams");
    }
}


CashFlowArraySP BondParams::getCashFlows(const DateTime &startDate) const
{
    static const string method = "BondParams::getCashFlows_DateTime";
    try {  
        throwIfNotValid();
        CashFlowArraySP myCashFlows(new CashFlowArray());
        int i;
        int numCFs=0;

        if (startDate <= getMaturityDate()) {

            if (startDate >= exCouponDates->back()) {
                // special case. return principal but no coupon
                myCashFlows->resize(1);
                (*myCashFlows)[0].date = getMaturityDate();
                (*myCashFlows)[0].amount = getRedemption();
            } else {
                // count the cashflows after the startDate
              
                while(numCFs < cashFlows->size() 
                      && (*exCouponDates)[cashFlows->size() - numCFs - 1] > startDate) {
                    numCFs += 1;
                }
                myCashFlows->resize(numCFs);
                for (i = 0; i<numCFs; i++) {
                    (*myCashFlows)[numCFs-i-1] = (*cashFlows)[cashFlows->size()-i-1];
                }
            }
        }

        return myCashFlows;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

}

double BondParams::getAccruedAtDate(const DateTime &aiDate) const 
{
    double accrued=0.0;
    int idxNextCF;
    int i;

    throwIfNotValid();

    if (aiDate <= datedDate || aiDate >= maturityDate) {
        accrued = 0.0;
    } else {
        // figure out where we fall in the schedule
        idxNextCF = 0;
        while((*schedCouponDates)[idxNextCF] <= aiDate) {
            idxNextCF += 1;
        }
        
        // Determine parameters needed by accruedFactor()
        // For now accruedFactor() is only used for REGULAR coupon periods because irregular periods need to agree with
        // BondParams::initialize() which has not yet been updated
        DateTime prevDate, nextDate;
        bool goneEx = false;
        bool regularCoupon = false;  
        if (aiDate >= (*exCouponDates)[idxNextCF]) {
            goneEx = true;
        }
        
        if (idxNextCF == 0 && cashFlows->size() != 1) {
             
            if (oddFirstLong == false) {
                prevDate = datedDate;
                nextDate = (*schedCouponDates)[idxNextCF];
                regularCoupon = true;
                
            } else {
                // Irregular period - figure out where we fall in the pseudo front schedule
                int idxNextPseudoCF;
                idxNextPseudoCF = 0;
                while((*pseudoFrontFlows)[idxNextPseudoCF].date <= aiDate) {
                    idxNextPseudoCF += 1;
                }
                
                if (idxNextPseudoCF > 0) {
                    accrued = 0;
                    for (i=0; i<idxNextPseudoCF; i++) {
                        accrued += (*pseudoFrontFlows)[i].amount;
                    }
                    accrued += (*pseudoFrontFlows)[idxNextPseudoCF].amount * 
                        (*accruedDCC).days((*pseudoFrontFlows)[idxNextPseudoCF-1].date, aiDate, endOfMonthAdj, eomIgnoreLeapYear)/
                        (*accruedDCC).days((*pseudoFrontFlows)[idxNextPseudoCF-1].date, (*pseudoFrontFlows)[idxNextPseudoCF].date, endOfMonthAdj, eomIgnoreLeapYear);
                } else {
                    accrued = (*pseudoFrontFlows)[0].amount * 
                        (*accruedDCC).days(datedDate, aiDate, endOfMonthAdj, eomIgnoreLeapYear)/
                        (*accruedDCC).days(datedDate, (*pseudoFrontFlows)[0].date, endOfMonthAdj, eomIgnoreLeapYear);
                }
                if (aiDate >= (*exCouponDates)[idxNextCF]) {
                   // subtract off the next coupon amount if it's already gone ex
                    accrued -= (*cashFlows)[idxNextCF].amount;
                }
            }
        
        } else if (idxNextCF == cashFlows->size()-1) {
            if (oddLastLong == false) {
                if (cashFlows->size() == 1) {
                    
                    if ((*accruedDCC).days(datedDate, (*schedCouponDates)[idxNextCF], endOfMonthAdj, eomIgnoreLeapYear) != 0 ) { 
                        prevDate = datedDate;
                        nextDate = (*schedCouponDates)[idxNextCF];
                        regularCoupon = true;
                    } else {
                        accrued = 0.0;
                    }
                } else {
                    if ((*accruedDCC).days((*schedCouponDates)[idxNextCF-1], (*schedCouponDates)[idxNextCF], endOfMonthAdj, eomIgnoreLeapYear) != 0 ) { 
                        prevDate = (*schedCouponDates)[idxNextCF-1];
                        nextDate = (*schedCouponDates)[idxNextCF];
                        regularCoupon = true;
                    } else {
                        accrued = 0.0;
                    }
                }
            } else {
                // Irregular period - figure out where we fall in the pseudo back schedule
                int idxNextPseudoCF;
                idxNextPseudoCF = 0;
                while((*pseudoBackFlows)[idxNextPseudoCF].date <= aiDate) {
                    idxNextPseudoCF += 1;
                }
                
                for (i=1; i<idxNextPseudoCF; i++) {
                    accrued += (*pseudoBackFlows)[i].amount;
                }
                accrued += (*pseudoBackFlows)[idxNextPseudoCF].amount * 
                    (*accruedDCC).days((*pseudoBackFlows)[idxNextPseudoCF-1].date, aiDate, endOfMonthAdj, eomIgnoreLeapYear)/
                    (*accruedDCC).days((*pseudoBackFlows)[idxNextPseudoCF-1].date, (*pseudoBackFlows)[idxNextPseudoCF].date, endOfMonthAdj, eomIgnoreLeapYear);
                
                if (aiDate >= (*exCouponDates)[idxNextCF]) {
                    // subtract off the next coupon amount if it's already gone ex
                    accrued -= ((*cashFlows)[idxNextCF].amount - getRedemption());
                }
            }

        } else {
            prevDate = (*schedCouponDates)[idxNextCF-1];
            nextDate = (*schedCouponDates)[idxNextCF];
            regularCoupon = true;
        }

        if (regularCoupon == true) {
            accrued = couponPct * faceValue * (1.-taxRate) * 
                            accruedDCC->accruedFactor(aiDate, 
                                                      prevDate, 
                                                      nextDate,
                                                      frequency,
                                                      "REGULAR",
                                                      endOfMonthAdj,
                                                      goneEx);
        }
    }        

    return accrued;
}

double BondParams::getFaceValue() const {
    return faceValue;
}

double BondParams::getNotional(const DateTime &ntlDate) const 
{
    double notional;

    throwIfNotValid();
 
    if (isAnOID == false){
        notional = getFaceValue();
    } else {
        if (ntlDate >= maturityDate) {
            notional = getFaceValue();
        } else if (ntlDate < datedDate) {
            notional = 0; // shouldn't matter what we set it to here
        } else {
            // notional = priceFromYield(yieldForOID, true, ntlDate);

            // calculate notional based on a zero coupon bond and semi-annual compounding
            double freqToUse = double(CompoundBasis::SEMI_ANNUAL);
            notional = faceValue;
            notional /= pow((1+yieldForOID/freqToUse), (*accruedDCC).years(ntlDate, maturityDate)*freqToUse);
        }
    }

    return notional;
}

double BondParams::getRedemption() const {
    return faceValue*redemptionPct;
}

double BondParams::getMaturityCoupon() const
{
    throwIfNotValid();

    return (cashFlows->back().amount - getRedemption());
}
  
DateTime BondParams::getMaturityDate() const
{
    throwIfNotValid();
    // return the last coupon date. Maturity can be adjusted for bad days.
    // Yes, people do schedule maturity on a bad day. 
    return cashFlows->back().date;
}

DateTime BondParams::getUnadjMaturityDate() const
{
    throwIfNotValid();
    return maturityDate;
}

DateTimeArraySP BondParams::getExCouponDates() const
{
    throwIfNotValid();
 
    return exCouponDates;
}


void BondParams::validatePop2Object()
{
    // convert the day count convention strings to objects - taken out of the try/catch block as this
    // was causing a crash downstream (in intialize())
    badDayConv = BadDayConventionSP(BadDayConventionFactory::make(badDayConvString));
    accruedDCC = DayCountConventionSP(DayCountConventionFactory::make(dayCountConvString));

    // if no settlement (as we only introduced it in May 2004), invent one
    if (!settle.get()) {
        settle = SettlementSP(new RollingSettlement);
    }

    try{
        validateInputs();

        // initialize();
        // isValid = true;
    }
    catch (exception& e) {
        char *tmpString = ModelException(e).stackTrace();
        errorMessage = tmpString;
        free(tmpString);
        errorMessage = errorMessage + "\n" + "BondParams::validatePop2Object: Error in object. Continuing but will fail if used";
        isValid = false;
    }
    
    return;
}
      

void BondParams::throwIfNotValid() const
{
    if (isValid == false) {
        throw ModelException("BondParams::throwIfNotValid", 
                             "attempting to use invalid bond (" + errorMessage + ")");
    }
}

void BondParams::validateInputs() const
{
    static const string method = "BondParams::validateInputs";
    try {
        if (!(frequency == 0 || frequency == 1 || frequency == 2 || frequency == 4 || frequency == 12)) {
            throw ModelException(method, "frequency (" + Format::toString(frequency) + ") must be 0, 1, 2, 4, or 12.");
        }
        if (faceValue < 0.) {
            throw ModelException(method, "face value (" + Format::toString(faceValue) + ") cannot be negative");
        } 
        /*
        if (couponPct < 0.) {
            throw ModelException(method, "coupon (" + Format::toString(couponPct) + ") cannot be negative");
        } 
        */
        if (!(firstCouponDate.empty())) {
            if (firstCouponDate <= datedDate) {
                throw ModelException(method, "firstCouponDate (" + firstCouponDate.toString() + 
                    ") must be greater than datedDate (" + datedDate.toString() + ")");
            }
            if (firstCouponDate > maturityDate) {
                throw ModelException(method, "firstCouponDate (" + firstCouponDate.toString() + 
                    ") must be less than or equal to maturityDate (" + maturityDate.toString() + ")");
            }
        }
        if (exDivDays < 0) {
            throw ModelException(method, "exDivDays (" + Format::toString(exDivDays) + ") must be >= 0");
        } 
        if (!(exDivRule == "N" || exDivRule == "B" || exDivRule == "C")) {
            throw ModelException(method, "exDivRule (" + exDivRule + ") must be N, B, or C");
        } 
        if (taxRate < 0.) {
            throw ModelException(method, "tax rate (" + Format::toString(taxRate) + ") cannot be negative");
        } else if (taxRate > 1.) {
            throw ModelException(method, "tax rate (" + Format::toString(taxRate) + ") cannot be > 1");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// This routine generates the complete cashFlow list from the bond parameters.
// The logic for building the unadjusted for bad day cashFlow dates is roughly as follows: 
// If the firstCouponDate is empty and oddLastShort is false:
//    buid from maturity back to datedDate
// Else if the firstCouponDate is empty and oddLastShort is true:
//    build from the datedDate up to maturity
// Else if the firstCouponDate is not empty build from the firstCouponDate to maturity
//    Also make the last coupon long or short depending on the oddLastShort flag
// FOR THE FUTURE it may be good to see if the bad day adjusted firstCouponDate is 
//    on cycle with the schedule made if it were not set. Could avoid some grief as
//    it is often set incorrectly to an adjusted date in IMS.
// If endOfMonthAdjust is true adjust dates to eom if any date falls at eom
// Set the coupon amounts
// Adjust the dates for bad days. 
// Set the principal date and amount
void BondParams::initialize()
{
    validateInputs();

    static const string method = "BondParams::initialize";

    try {  
        int i;
        int numCFs;
        int interval;

        schedCouponDates = DateTimeArraySP(new DateTimeArray());
        exCouponDates = DateTimeArraySP(new DateTimeArray());
        cashFlows = CashFlowArraySP(new CashFlowArray());
     
        // handle zeros separately
        if (frequency != 0 && !Maths::isZero(couponPct)) { // we are not a zero
            
            interval = 12/frequency ;

            if (firstCouponDate.empty() == true) {

                if (oddLastShort == false) {
                    // start at maturity and work back
                    numCFs = 1;
                    while( MaturityPeriod::toDate(-numCFs*interval, "M", maturityDate).getDate() > datedDate.getDate()) {
                        numCFs += 1;
                    }
                    schedCouponDates->resize(numCFs);
                    // set the dates. 
                    for (i = 0; i<numCFs; i++) {
                        (*schedCouponDates)[numCFs-1-i] = 
                            MaturityPeriod::toDate(-i*interval, "M", maturityDate);
                    }
                }
                else { // oddLastShort == true
                    // start at datedDate and work forward
                    numCFs = 1;
                    while(MaturityPeriod::toDate(numCFs*interval, "M", datedDate) < maturityDate) {
                        numCFs += 1;
                    }
                    
                    schedCouponDates->resize(numCFs);
                    // set the dates. 
                    for (i = 0; i<numCFs-1; i++) {
                        (*schedCouponDates)[i] = 
                            MaturityPeriod::toDate((i+1)*interval, "M", datedDate);
                    }
                    
                    (*schedCouponDates)[numCFs-1] = maturityDate;
                }
            } else { // firstCouponDate has been set
                numCFs = 1;
                while(MaturityPeriod::toDate(numCFs*interval, "M", firstCouponDate).getDate() < maturityDate.getDate()) {
                    numCFs += 1;
                }

                if (MaturityPeriod::toDate(numCFs*interval, "M", firstCouponDate).getDate() == maturityDate.getDate()) {
                    numCFs += 1;
                } else { // we've got an odd last coupon
                    if (oddLastShort == true) {
                        numCFs += 1;
                    }
                } 

                schedCouponDates->resize(numCFs);
                for (i = 0; i<numCFs-1; i++) {
                    (*schedCouponDates)[i] = 
                        MaturityPeriod::toDate((i)*interval, "M", firstCouponDate);
                }
                (*schedCouponDates)[numCFs-1] = maturityDate;
            } 
            
            // adjust for end of month
            bool adjEOM = false;
            if (endOfMonthAdj == true) {
                
                
                for (i = 0; i<numCFs-1; i++) { // never adjust maturity
                    if((*schedCouponDates)[i].isEndOfMonth()){
                        adjEOM = true;
                        break;
                    }
                }
                // handle case where only two coupons and last date falls at eom but first does not -- first should be moved to eom as long as last coupon isn't odd
                if (adjEOM == false && maturityDate.isEndOfMonth() && numCFs > 1) {// don't bother if it's already true
                    if (MaturityPeriod::toDate(interval, "M", (*schedCouponDates)[numCFs-2]).getDate() == maturityDate.getDate()) {
                        adjEOM = true;
                    }
                }
                        
                if (adjEOM == true){  // never adjust maturity
                    for (i = 0; i<numCFs-1; i++) {
                        (*schedCouponDates)[i] = (*schedCouponDates)[i].returnEndOfMonth(eomIgnoreLeapYear);
                    }
                }
                
                // we could have put an extra coupon at maturity if adjEOM is true as we compared
                // unadjusted dates above. Fix that here.
                if (numCFs > 1 && (*schedCouponDates)[numCFs-1].getDate() == (*schedCouponDates)[numCFs-2].getDate()) {
                    numCFs--;
                    schedCouponDates->pop_back(); // removes the redundant last date
                }
            }

            cashFlows->resize(numCFs);
            // set the coupons amounts
            if (numCFs >= 2)
            {
                // handle the first and last separately
                if ((*schedCouponDates)[0] == (adjEOM ? 
                    MaturityPeriod::toDate(interval, "M", datedDate).returnEndOfMonth(eomIgnoreLeapYear) :
                    MaturityPeriod::toDate(interval, "M", datedDate))) {
                    (*cashFlows)[0].amount = couponPct * faceValue * (1.-taxRate) / (double) frequency;
                    oddFirst = false;
                } else { // odd first coupon
                    oddFirst = true;            
                    int numPseudoCFs;
                    numPseudoCFs = 0;

                    while( (*schedCouponDates)[0] > MaturityPeriod::toDate(interval, "M", datedDate) && 
                        MaturityPeriod::toDate(-(numPseudoCFs+1)*interval, "M", (*schedCouponDates)[0]).isGreaterOrEqual(datedDate)) {
                        numPseudoCFs += 1;
                    }
                    DateTime pseudoFirstDate;


                    // the following line should get replaced by the commented out bit
                    pseudoStartDate = MaturityPeriod::toDate(-(numPseudoCFs+1)*interval, "M", (*schedCouponDates)[0]);

                    pseudoFirstDate = MaturityPeriod::toDate(-(numPseudoCFs)*interval, "M", (*schedCouponDates)[0]);
                    if (adjEOM == true) {
                        pseudoStartDate = pseudoStartDate.returnEndOfMonth(eomIgnoreLeapYear);
                        pseudoFirstDate = pseudoFirstDate.returnEndOfMonth(eomIgnoreLeapYear);
                    }
                    
                    (*cashFlows)[0].amount = couponPct * faceValue * (1.-taxRate) * 
                        ((double) numPseudoCFs / (double) frequency + 
                        (double)(*accruedDCC).days(datedDate, pseudoFirstDate, endOfMonthAdj, eomIgnoreLeapYear)/
                        (double)(*accruedDCC).days(pseudoStartDate, pseudoFirstDate, endOfMonthAdj, eomIgnoreLeapYear)/
                        (double) frequency);
                    
                    if (numPseudoCFs > 0) {
                        oddFirstLong = true;
                        pseudoFrontFlows = CashFlowArraySP(new CashFlowArray());
                        pseudoFrontFlows->resize(numPseudoCFs+1); // an extra one for the first real coupon
                        for (i = 0; i<numPseudoCFs+1; i++) {
                            (*pseudoFrontFlows)[numPseudoCFs-i].date = 
                             MaturityPeriod::toDate(-i*interval, "M", (*schedCouponDates)[0]);
                        }
                        if (adjEOM == true){  // never adjust first real coupon date
                            for (i = 0; i<numPseudoCFs; i++) {
                                (*pseudoFrontFlows)[i].date = (*pseudoFrontFlows)[i].date.returnEndOfMonth(eomIgnoreLeapYear);
                            }
                        }
                        for (i = numPseudoCFs; i>0; i--) {
                            (*pseudoFrontFlows)[i].amount = couponPct * faceValue * (1.-taxRate) / (double) frequency;
                        }

                        (*pseudoFrontFlows)[0].amount = couponPct * faceValue * (1.-taxRate) * 
                            (double)(*accruedDCC).days(datedDate, pseudoFirstDate, endOfMonthAdj, eomIgnoreLeapYear)/
                            (double)(*accruedDCC).days(pseudoStartDate, pseudoFirstDate, endOfMonthAdj, eomIgnoreLeapYear)/
                            (double) frequency;
                    }          
                }

                if ((adjEOM ? MaturityPeriod::toDate(interval, "M", (*schedCouponDates)[numCFs-2]).returnEndOfMonth(eomIgnoreLeapYear) :
                    MaturityPeriod::toDate(interval, "M", (*schedCouponDates)[numCFs-2])).getDate() == maturityDate.getDate()) {
                    (*cashFlows)[numCFs-1].amount = couponPct * faceValue * (1.-taxRate) / (double) frequency;
                    oddLast = false;
                } else { // odd last coupon
                    oddLast = true;
                    int numPseudoCFs;
                    numPseudoCFs = 0;

                    while(MaturityPeriod::toDate((numPseudoCFs+1)*interval, "M", (*schedCouponDates)[numCFs-2]).getDate() <= maturityDate.getDate()) {
                        numPseudoCFs += 1;
                    }

                    DateTime pseudoPenultDate;
                    pseudoMatDate = MaturityPeriod::toDate((numPseudoCFs+1)*interval, "M", (*schedCouponDates)[numCFs-2]);
                    pseudoPenultDate = MaturityPeriod::toDate(numPseudoCFs*interval, "M", (*schedCouponDates)[numCFs-2]);
                    if (adjEOM == true) {
                        pseudoMatDate = pseudoMatDate.returnEndOfMonth(eomIgnoreLeapYear);
                        pseudoPenultDate = pseudoPenultDate.returnEndOfMonth(eomIgnoreLeapYear);
                    }

                    (*cashFlows)[numCFs-1].amount = couponPct * faceValue * (1.-taxRate) * 
                        ((double) numPseudoCFs + 
                        (double) (*accruedDCC).days(pseudoPenultDate, maturityDate, endOfMonthAdj, eomIgnoreLeapYear)/
                        (double) (*accruedDCC).days(pseudoPenultDate, pseudoMatDate, endOfMonthAdj, eomIgnoreLeapYear))/
                        (double) frequency;
                    
                    if (numPseudoCFs > 0) {
                        oddLastLong = true;
                        pseudoBackFlows = CashFlowArraySP(new CashFlowArray());
                        pseudoBackFlows->resize(numPseudoCFs+2); // extras for true penultimate coupon and maturity coupon
                        for (i = 0; i<numPseudoCFs+1; i++) {
                            (*pseudoBackFlows)[i].date = 
                             MaturityPeriod::toDate(i*interval, "M", (*schedCouponDates)[numCFs-2]);
                        }
                        (*pseudoBackFlows)[numPseudoCFs+1].date = maturityDate;

                        if (adjEOM == true){  // don't adjust real penultimate or maturity coupon dates
                            for (i = 1; i<numPseudoCFs+1; i++) {
                                (*pseudoBackFlows)[i].date = (*pseudoBackFlows)[i].date.returnEndOfMonth(eomIgnoreLeapYear);
                            }
                        }
                        (*pseudoBackFlows)[0].amount = -999.; //should never be used
                        for (i = 1; i<numPseudoCFs+1; i++) {
                            (*pseudoBackFlows)[i].amount = couponPct * faceValue * (1.-taxRate) / (double) frequency;
                        }
                        (*pseudoBackFlows)[numPseudoCFs+1].amount = couponPct * faceValue * (1.-taxRate) * 
                            (double) (*accruedDCC).days(pseudoPenultDate, maturityDate, endOfMonthAdj, eomIgnoreLeapYear)/
                            (double) (*accruedDCC).days(pseudoPenultDate, pseudoMatDate, endOfMonthAdj, eomIgnoreLeapYear)/
                            (double) frequency;
                    }   
                }

                for (i = 1; i<numCFs-1; i++) {
                    (*cashFlows)[i].amount = couponPct * faceValue * (1.-taxRate) / (double) frequency;
                }
            } else if (numCFs == 1) { // only one coupon
                // could make pseudo CFs here -- but what's the use as I've never seen an odd long 1 coupon bond
                //     that follows corporate conventions. Error is minimal and alib would fail with bizarre 
                //     messages if you even tried it
                (*cashFlows)[numCFs-1].amount = couponPct * faceValue * (1.-taxRate) * 
                    (*accruedDCC).years(datedDate, maturityDate);
            }
            
            // set the cashFlow dates -- adjusted for bad days
            for (i = 0; i<numCFs; i++) {
                (*cashFlows)[i].date = badDayConv->adjust((*schedCouponDates)[i], hols.get());
            }
                    
        } else { // a zero
            numCFs = 1;
            schedCouponDates->resize(1);
            cashFlows->resize(numCFs);
            (*schedCouponDates)[0] = maturityDate;
            if (!Maths::isZero(couponPct)) {
                (*cashFlows)[0].amount = couponPct * faceValue * (1.-taxRate) * 
                    (*accruedDCC).years(datedDate, maturityDate);
            }
            (*cashFlows)[0].date = badDayConv->adjust((*schedCouponDates)[0], hols.get());
        }
 
        // get the exCouponDates
        exCouponDates->resize(numCFs);
        if (exDivRule == "N" || exDivDays == 0) {
             for (i = 0; i<numCFs; i++) {
                (*exCouponDates)[i] = (*schedCouponDates)[i];
            }
        } else if (exDivRule == "C") {
            for (i = 0; i<numCFs; i++) {
                (*exCouponDates)[i] = (*schedCouponDates)[i].rollDate(-1*exDivDays);
            }
        } else if (exDivRule == "B") {
            for (i = 0; i<numCFs; i++) {
                (*exCouponDates)[i] = hols->addBusinessDays((*schedCouponDates)[i], -1*exDivDays);
            }
        } else {
            throw ModelException(method, "exDivRule (" + exDivRule + ") must be N, B, or C");
        }
        // if the last ex-div date after the last cashFlow date (which can happen when the bad day 
        // convention is previous) move the last ex date to the effective maturity
        if ((*exCouponDates)[numCFs-1] > (*cashFlows)[numCFs-1].date) {
            (*exCouponDates)[numCFs-1] = (*cashFlows)[numCFs-1].date;
        }

        // set the principal
        (*cashFlows)[numCFs-1].amount += faceValue * redemptionPct;
        isValid = true;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

    return;
}

double BondParams::priceFromYield(double yield, bool clean, DateTime yDate) const
{
    throwIfNotValid();
    static const string method = "BondParams::priceFromYield";
    try {
        double myPrice;
        if (yDate > getMaturityDate()) {
            throw ModelException("startDate cannot be after maturityDate");
        }
        CashFlowArraySP cfs = getCashFlows(yDate);
        double freqToUse;
        int i;

        if (frequency == 0) {
            freqToUse = 2;
        } else {    
            freqToUse = (double) frequency;
        }

        if (frequency != 0 && !Maths::isZero(couponPct) && cashFlows->size() > 1) { // not a zero
            // handle the last one specially in case it's odd
            myPrice = (*cfs)[i=cfs->size()-1].amount;
            if (oddLast == true) {
                if (oddLastLong == true) {
                    if (yDate.getDate() <= (*pseudoBackFlows)[pseudoBackFlows->size()-2].date.getDate()) {
                        myPrice /= pow((1+yield/freqToUse),
                            (double) (*accruedDCC).days((*pseudoBackFlows)[pseudoBackFlows->size()-2].date, maturityDate, endOfMonthAdj, eomIgnoreLeapYear)/
                            (double) (*accruedDCC).days((*pseudoBackFlows)[pseudoBackFlows->size()-2].date, pseudoMatDate, endOfMonthAdj, eomIgnoreLeapYear));
                    } else {
                        myPrice /= pow((1+yield/freqToUse),
                            (double) (*accruedDCC).days(yDate, maturityDate, endOfMonthAdj, eomIgnoreLeapYear)/
                            (double) (*accruedDCC).days((*pseudoBackFlows)[pseudoBackFlows->size()-2].date, pseudoMatDate, endOfMonthAdj, eomIgnoreLeapYear));
                    } 
                        
                    for (i=pseudoBackFlows->size()-3; i>=0; i--) { // you're limited to one pseudo coupon at present, but we may want to add more later
                        if (yDate.getDate() <= (*pseudoBackFlows)[i].date.getDate()) {
                            myPrice /= (1+yield/freqToUse);
                        } else {
                            myPrice /= pow((1+yield/freqToUse),
                                (double) (*accruedDCC).days(yDate, (*pseudoBackFlows)[i+1].date, endOfMonthAdj, eomIgnoreLeapYear)/
                                (double) (*accruedDCC).days((*pseudoBackFlows)[i].date, (*pseudoBackFlows)[i+1].date, endOfMonthAdj, eomIgnoreLeapYear));
                        }
                    }
                } else {
                    if (cfs->size() > 1) {
                        myPrice /= pow((1+yield/freqToUse),
                            (double) (*accruedDCC).days((*cfs)[cfs->size()-2].date, maturityDate, endOfMonthAdj, eomIgnoreLeapYear)/
                            (double) (*accruedDCC).days((*cfs)[cfs->size()-2].date, pseudoMatDate, endOfMonthAdj, eomIgnoreLeapYear));
                    } else {
                        myPrice /= pow((1+yield/freqToUse),
                            (double) (*accruedDCC).days(yDate, maturityDate, endOfMonthAdj, eomIgnoreLeapYear)/
                            (double) (*accruedDCC).days((*cashFlows)[cashFlows->size()-2].date, pseudoMatDate, endOfMonthAdj, eomIgnoreLeapYear));
                    }
                }
            } else {
                // on schedule
                if (cfs->size() > 1) {
                    myPrice /= (1+yield/freqToUse);
                } else {
                    myPrice /= pow((1+yield/freqToUse),
                        (double) (*accruedDCC).days(yDate, maturityDate, endOfMonthAdj, eomIgnoreLeapYear)/
                        (double) (*accruedDCC).days((*cashFlows)[cashFlows->size()-2].date, maturityDate, endOfMonthAdj, eomIgnoreLeapYear));
                }
            }

            // the middle coupons
            for (i=cfs->size()-2; i>0; i--) {
                myPrice += (*cfs)[i].amount;
                myPrice /= (1+yield/freqToUse);
            }

            // the first coupon if it's not also the last
            if (cfs->size() > 1) {
                myPrice += (*cfs)[0].amount;
                if (oddFirst == true && yDate.getDate() < (*exCouponDates)[0].getDate()) {
                    // we're in the odd front portion
                    if (oddFirstLong == true) {
                        for (i=pseudoFrontFlows->size()-2; i>=0 && yDate.getDate() < (*pseudoFrontFlows)[i].date.getDate(); i--) {
                            myPrice /= (1+yield/freqToUse);
                        }
                        if (i >= 0) {
                            myPrice /= pow((1+yield/freqToUse),
                                (double) (*accruedDCC).days(yDate, (*pseudoFrontFlows)[i+1].date, endOfMonthAdj, eomIgnoreLeapYear)/
                                (double) (*accruedDCC).days((*pseudoFrontFlows)[i].date, (*pseudoFrontFlows)[i+1].date, endOfMonthAdj, eomIgnoreLeapYear));
                        } else {
                            // we're between the pseudoStart and the first coupon
                            myPrice /= pow((1+yield/freqToUse),
                                (double) (*accruedDCC).days(yDate, (*pseudoFrontFlows)[i+1].date, endOfMonthAdj, eomIgnoreLeapYear)/
                                (double) (*accruedDCC).days(pseudoStartDate, (*pseudoFrontFlows)[i+1].date, endOfMonthAdj, eomIgnoreLeapYear));
                        }
                    } else {
                        myPrice /= pow((1+yield/freqToUse),
                            (double) (*accruedDCC).days(yDate, (*cfs)[0].date, endOfMonthAdj, eomIgnoreLeapYear)/
                            (double) (*accruedDCC).days(pseudoStartDate, (*cfs)[0].date, endOfMonthAdj, eomIgnoreLeapYear));
                    }
                } else {
                    myPrice /= pow((1+yield/freqToUse),((*cfs)[0].amount-getAccruedAtDate(yDate))/(*cfs)[0].amount); // not right for long front
                    // ((*cfs)[0].amount-getAccruedAtDate(yDate))/(*cfs)[0].amount = year fract remaining in first period
                }
            }
       } else {
            // a zero
            myPrice = (*cfs)[0].amount;
            myPrice /= pow((1+yield/freqToUse), (*accruedDCC).years(yDate, maturityDate)*freqToUse);
       }
       if (clean == true) {
           myPrice -= getAccruedAtDate(yDate);
       }
       return myPrice;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

Bond* BondParams::getBondToPut(const DateTime &putDate, double putLevel) const
{
    throwIfNotValid();
    static const string method = "BondParams::getBondToPut";
    try {
        BondParamsSP putBond(dynamic_cast<BondParams*>(clone()));
    
        putBond->maturityDate = putDate;
        putBond->redemptionPct = putLevel/faceValue;
        
        if (schedCouponDates->size() > 1) {
            if (firstCouponDate.empty() == true) {
                putBond->firstCouponDate = (*schedCouponDates)[0];
            }

            if (oddLastLong == true) {
                if (putDate <= (*pseudoBackFlows)[1].date) {
                    putBond->oddLastShort = true;
                } 
            } else {
                putBond->oddLastShort = true;
            }
        }
        // reset the derived members of putBond that need to be
        putBond->isValid = false;
        putBond->oddFirst = false;
        putBond->oddLast = false;
        putBond->oddLastLong = false;
        putBond->oddFirstLong = false;
        putBond->exCouponDates = DateTimeArraySP(   );
        putBond->schedCouponDates = DateTimeArraySP(   );
        putBond->cashFlows = CashFlowArraySP(   );
        putBond->pseudoFrontFlows = CashFlowArraySP(   );
        putBond->pseudoBackFlows = CashFlowArraySP(   );

        putBond->initialize();
        putBond->isValid = true;

        // need to subtract of the last coupon if it's not a full one.
        // the put should include accrued if it needs to
        (putBond->cashFlows)->back().amount -= getAccruedAtDate(putDate);

        return putBond.release();
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

double BondParams::getBondCouponRate() const
{
    return couponPct;
}

DateTime BondParams::getAccrualStartDate() const
{
    return datedDate;
}

CashFlowArraySP BondParams::getRedemptionPayments() const
{
    CashFlowArraySP redemptionPayments(new CashFlowArray(1));
    (*redemptionPayments)[0] = CashFlow(cashFlows->back().date,
                                        getRedemption());
    return redemptionPayments;
}

// when does bond settle?
DateTime BondParams::settles(const DateTime& tradeDate) const {
    return settle->settles(tradeDate);
}

class BondParamsHelper{
public:
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BondParams, clazz);
        SUPERCLASS(Bond);
        EMPTY_SHELL_METHOD(defaultBondParams);
        FIELD(faceValue, "face value");
        FIELD(redemptionPct, "redemption as a percent of accreated");
        FIELD_MAKE_OPTIONAL(redemptionPct)
        FIELD(couponPct, "coupon as percent of face value");
        FIELD(frequency, "coupon frequency (0,1,2,4,12)");
        FIELD(maturityDate, "maturity date");
        FIELD(datedDate, "date interest begins to accrue");
        FIELD(dayCountConvString, "day count convention");
        FIELD_MAKE_OPTIONAL(dayCountConvString);
        FIELD(firstCouponDate, "date first coupon is paid");
        FIELD_MAKE_OPTIONAL(firstCouponDate);
        FIELD(oddLastShort, "is the last coupon short")
        FIELD_MAKE_OPTIONAL(oddLastShort);
        FIELD(endOfMonthAdj, "put coupons on month ends if one falls on month end")
        FIELD_MAKE_OPTIONAL(endOfMonthAdj);
        FIELD(eomIgnoreLeapYear, "adjust coupons to 28th of Feb on leap years if true")
        FIELD_MAKE_OPTIONAL(eomIgnoreLeapYear);
        FIELD(badDayConvString, "bad day addjustment N, P, F, or M")
        FIELD_MAKE_OPTIONAL(badDayConvString);
        FIELD(hols, "holiday handle");
        FIELD(exDivDays, "number of days bond trades without coupon before coupon is paid")
        FIELD_MAKE_OPTIONAL(exDivDays);
        FIELD(exDivRule, "N = None, B = business days, C = CalenderDays")
        FIELD_MAKE_OPTIONAL(exDivRule);
        FIELD(taxRate, "tax rate")
        FIELD_MAKE_OPTIONAL(taxRate);
        FIELD(isAnOID, "was the bond originally issued at a discount")
        FIELD_MAKE_OPTIONAL(isAnOID);
        FIELD(yieldForOID, "yield for original issue discount")
        FIELD_MAKE_OPTIONAL(yieldForOID);
        FIELD(settle, "how bond settles");
        FIELD_MAKE_OPTIONAL(settle);
        // private members that need to be copied when bond is copied
        FIELD(isValid, "private member");
        FIELD_MAKE_TRANSIENT(isValid);
        FIELD(errorMessage, "private member");
        FIELD_MAKE_TRANSIENT(errorMessage);
        FIELD(cashFlows, "private member");
        FIELD_MAKE_TRANSIENT(cashFlows);
        FIELD(exCouponDates, "private member");
        FIELD_MAKE_TRANSIENT(exCouponDates);
        FIELD(schedCouponDates, "private member");
        FIELD_MAKE_TRANSIENT(schedCouponDates);
        FIELD(badDayConv, "private member");
        FIELD_MAKE_TRANSIENT(badDayConv);
        FIELD(accruedDCC, "private member");
        FIELD_MAKE_TRANSIENT(accruedDCC);
        FIELD(oddFirstLong, "private member");
        FIELD_MAKE_TRANSIENT(oddFirstLong);
        FIELD(oddLastLong, "private member");
        FIELD_MAKE_TRANSIENT(oddLastLong);
        FIELD(pseudoFrontFlows, "private member");
        FIELD_MAKE_TRANSIENT(pseudoFrontFlows);
        FIELD(pseudoBackFlows, "private member");
        FIELD_MAKE_TRANSIENT(pseudoBackFlows);
        FIELD(oddFirst, "private member");
        FIELD_MAKE_TRANSIENT(oddFirst);
        FIELD(oddLast, "private member");
        FIELD_MAKE_TRANSIENT(oddLast);
        FIELD(pseudoStartDate, "private member");
        FIELD_MAKE_TRANSIENT(pseudoStartDate);
        FIELD(pseudoMatDate, "private member");
        FIELD_MAKE_TRANSIENT(pseudoMatDate);

        // addin constructor
        Addin::registerConstructor("BOND_PARAMETERS",
                                   Addin::MARKET,
                                   "Creates a Bond specified via parameters",
                                   BondParams::TYPE);


    }

    static IObject* defaultBondParams(){
        return new BondParams();
    }

};

CClassConstSP const BondParams::TYPE = CClass::registerClassLoadMethod(
    "BondParams", typeid(BondParams), BondParamsHelper::load);
bool  BondParamsLoad() {
    return (BondParams::TYPE != 0);
   }




DRLIB_END_NAMESPACE
