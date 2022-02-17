//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DividendList.cpp
//
//   Description : Dividend List representation
//
//   Author      : Stephen Hope
//
//   Date        : 25 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_DIVIDENDLIST_CPP
#include "edginc/DividendList.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE
DividendList::IDivAdjuster::~IDivAdjuster(){}

/** Constructor */
DividendList::DividendList(const DividendArray&  divArray):
    CObject(TYPE), divArray(new DividendArray(divArray))
{
    validatePop2Object(); 
}

/** for reflection */
DividendList::DividendList(): CObject(TYPE)
{
    // empty
}

// PUBLIC METHODS

/** Validation checks + Re-ordering
    Conditions:
    1.  Increasing dates (not strictly)
    2.  Time of day == START_OF_DAY_TIME
    Re-ordering:
        All divs falling on the same day/time are reordered in the following order:
            1. PERCENT
            2. AMOUNT
            3. CONTINUOUS */
void DividendList::validatePop2Object()
{
    static const string routine = "DividendList::validatePop2Object";
    
    DividendArray& divs = *divArray;
    int NbDivs = divs.size();
    if (!NbDivs) return;

    /* Validate */
    DateTime firstContDivDate(0, 0);
    DateTime zeroDate(0, 0);
    int iDiv = 0;
    for (;iDiv < NbDivs; ++iDiv)
    {
        if ( divs[iDiv].getExDate().getTime() != Dividend::DIVIDEND_EXDIV_TIME )
        {
            throw ModelException(routine,
                                 "Dividend at ex div date " + 
                                 divs[iDiv].getExDate().toString() +
                                 " is not going ex-div at ex-div time");
        }

        if (iDiv < divs.size() - 1){
            if (divs[iDiv].getExDate().isGreater(divs[iDiv + 1].getExDate()))
            {
                throw ModelException(routine,
                                     "exDivDate " + 
                                     divs[iDiv].getExDate().toString() +
                                     " is greater than exDivDate " + 
                                     divs[iDiv+1].getExDate().toString());
            }
        }

        /* check dividend type */
        switch (divs[iDiv].getDivType())
        {
        case Dividend::AMOUNT:
        case Dividend::PERCENT:
            if (firstContDivDate.isGreater(zeroDate) &&
                firstContDivDate.isLess(divs[iDiv].getExDate())){
                /* All cont divs must fall after all discrete divs */
                throw ModelException(routine,
                                     "Discrete div paid " +
                                     divs[iDiv].getExDate().toString() + 
                                     " after continuous div " +
                                     firstContDivDate.toString());
            }
            break;
            
        case Dividend::CONTINUOUS:
            /* Record first cts div for check above. */
            if ( firstContDivDate.equals(zeroDate)){
                firstContDivDate = divs[iDiv].getExDate();
            }
            
            if (!(divs[iDiv].getExDate().equals(divs[iDiv].getPayDate()))){    
                throw ModelException(routine,
                                     "exDivDate " +
                                     divs[iDiv].getExDate().toString() +
                                     " is not equal to payDate " +
                                     divs[iDiv].getPayDate().toString() +
                                     ". Continuous dividends require ex-div "
                                     "and payment dates to be equal");
            }            
            break;            
        }
    }

    /* Do the sorting */
    sort();
}

/** Sorting method. 
    Re-order order of dividends falling on a same date / time in the following order
        1. PERCENT
        2. AMOUNT
        3. CONTINUOUS
    NB This routine does not change the order in which continuous dividends are inputted,
    so if several continuous dividends are inputted on the same date / time they will
    all but the last be ignored. */
void DividendList::sort()
{
    static const string routine = "DividendList::sort";
    
    DividendArray& divs = *divArray;
    int NbDivs = divs.size();
    if (!NbDivs) return;

    CIntArray divTypes(NbDivs);
    int iDiv = 0;
    int iDivBegin = iDiv;
    CDateTime exDivDateBegin(divs[iDiv].getExDate());
    for (;iDiv < NbDivs; ++iDiv)
    {
        switch (divs[iDiv].getDivType())
        {
        case Dividend::PERCENT:
            divTypes[iDiv] = 1;
            break;            
        case Dividend::AMOUNT:
            divTypes[iDiv] = 2;
            break;
        case Dividend::CONTINUOUS:
            divTypes[iDiv] = 3;
            break;            
        }

        if (iDiv > 0)
        {
            if (divs[iDiv].getExDate().isLess(exDivDateBegin))
            {
                throw ModelException(routine,
                                     "exDivDateBegin " + 
                                     exDivDateBegin.toString() +
                                     " is greater than exDivDate " + 
                                     divs[iDiv+1].getExDate().toString());
            }

            if(divs[iDiv].getExDate().isGreater(exDivDateBegin))
            {
                // does the sorting, if needed
                if( iDiv > iDivBegin + 1){
                    piksr2(iDivBegin,
                           iDiv, 
                           divTypes, 
                           divs);
                }
                iDivBegin = iDiv;
                exDivDateBegin = divs[iDiv].getExDate();        
            }
        }
    }
    if( iDiv > iDivBegin + 1){
        piksr2(iDivBegin,
               iDiv, 
               divTypes, 
               divs);
    }
}

void DividendList::piksr2(int iBegin, 
                          int iEnd, 
                          CIntArray& divTypes, 
                          CDividendArray& divs)
{
    int i, j;
    for (j = iBegin + 1; j < iEnd; j++) { // Pick out each element in turn.
        int divType = divTypes[j];
        Dividend div = divs[j];     // Should I use a ref ?
        i = j - 1;
        while (i >= iBegin && divTypes[i] > divType) { // Look for the place to insert it.
            divTypes[i + 1] = divTypes[i];
            divs[i + 1] = divs[i];
            i--;
        }
        divTypes[i + 1] = divType; // Insert it.
        divs[i + 1] = div;
    }
}

/** Given a dividend list, returns a dividend list where all dollar / yield 
    dividends falling on a same date / time are aggregated. If more than 1
    continuous dividends fall on the same date / time, only the last dividend
    is kept. */
DividendList* DividendList::aggregate(const DividendList& input,
                                      const IYieldCurve*  yc)
{
    try{
//        const DividendArray& inputArray = input.getArray();
        DividendArray inputArray(input.getArray());
        if (inputArray.size() == 0){
            return copy(&input);
        }

        DividendArray outputArray(0);
        outputArray.reserve(inputArray.size());

        auto_ptr<YieldCurve::IKey> key(yc->logOfDiscFactorKey());

        int iDiv = 0;
        while (iDiv < inputArray.size())
        {
            const DateTime& exDivDate = inputArray[iDiv].getExDate();

            double yield = 1.0;
            double absolute = 0.0;
            int ctsIdx = 0;
            bool hasYield = false;
            bool hasAbsolute = false;
            bool hasCts = false;

            while (iDiv < inputArray.size() && exDivDate.equals(inputArray[iDiv].getExDate())){
                const Dividend& dividend = inputArray[iDiv];

                double amount = dividend.getDivAmount();

                if (dividend.getPayDate().isGreater(exDivDate)) {
                    amount *= exp(key->calc(exDivDate, dividend.getPayDate()));
                }

                switch (dividend.getDivType())
                {
                case Dividend::PERCENT:
                    if (!Maths::isZero(amount)){
                        hasYield = true;
                    }
                    /* I don't understand this. If I have several yields at the same date / time,
                    they should sum up and the fwd should drop by "1.0 - sum yield". The fwd method 
                    however is implemented in such way that the fwd drops by "prod (1.0 - yield)".
                    I will stick to that for consistency reasons. */
                    yield *= 1.0 - amount;
                    break;            
                case Dividend::AMOUNT:
                    if (!Maths::isZero(amount)){
                        hasAbsolute = true;
                    }
                    absolute += amount;
                    break;
                /* The input list is assumed to be sorted, so we only need keep track of the last
                   continous dividend. All other continous dividends that are possibly present will be ignored. */
                case Dividend::CONTINUOUS:
                    hasCts = true;
                    ctsIdx = iDiv;
                    break;            
                }
                ++iDiv;
            }

            /* Watch out the order we output the divs when falling on a same date / time. */
            if (hasYield){
                outputArray.push_back(Dividend(exDivDate,
                                               exDivDate,          // discounting already taken care of
                                               Dividend::PERCENT,
                                               1.0 - yield));
            }

            if (hasAbsolute){
                outputArray.push_back(Dividend(exDivDate,
                                               exDivDate,          // discounting already taken care of
                                               Dividend::AMOUNT,
                                               absolute));
            }

            if (hasCts){
                outputArray.push_back(Dividend(exDivDate,
                                               inputArray[ctsIdx].getPayDate(),
                                               Dividend::CONTINUOUS,
                                               inputArray[ctsIdx].getDivAmount()));
            }
        }

        return new DividendList(outputArray);
    }
    catch(exception& e){
        throw ModelException(&e, "DividendList::aggregate");
    }
}

/** Create a DividendList full of all divs due to be paid > startDate, 
    and <= endDate. Returns an empty dividend list if no dividends match 
    this condition. It is the callers responsibility to free the memory 
    allocated by this function */
DividendList* DividendList::getAllDivsBetweenDates(
    const DateTime&   startDate,         // (I) No divs before this date 
    const DateTime&   endDate)const      // (I) No divs after this date 
{
    DividendListSP  resultDivs(new DividendList());
    resultDivs->divArray = DividendArraySP(new DividendArray(0));
    int startIdx;
    int  endIdx;
    int  idx;
    const DividendArray& divs = *divArray;
    int  numDivs = divs.size();
    if (numDivs > 0 )
    {
        startIdx = findDivArrayIndexFromDate(startDate);
        if (startIdx > -1)
        {
            for (endIdx = startIdx; endIdx < numDivs; endIdx++)
            {
                if (divs[endIdx].getExDate().isGreater(endDate))
                {
                    break;
                }
            }

            // We now know the end index.
            // Note that we include div start and exclude div end since 
            // findDivArrayIndexFromDate() gets idx for div due to be 
            // paid strictly after the date given
            // now copy the divs to the result list
            for (idx = startIdx; idx < endIdx; idx++)
            {
                resultDivs->divArray->push_back(divs[idx]);
            }
        }
    }

    return resultDivs.release();
}

/** same as getAllDivsBetweenDates, but ignores continuous dividends */
DividendList* DividendList::getAllDiscreteDivsBetweenDates(
    const DateTime&   startDate,         // (I) No divs before this date 
    const DateTime&   endDate) const     // (I) No divs after this date */
{
    DividendListSP  resultDivs(new DividendList());
    resultDivs->divArray = DividendArraySP(new DividendArray(0));
    int startIdx;
    int  endIdx;
    int  idx;
    const DividendArray& divs = *divArray;
    int  numDivs = divs.size();
    if (numDivs > 0 ) {
        startIdx = findDivArrayIndexFromDate(startDate);
        if (startIdx > -1) {
            for (endIdx = startIdx; endIdx < numDivs; endIdx++) {
                if (divs[endIdx].getExDate().isGreater(endDate)) {
                    break;
                }
            }

            // We now know the end index.
            // Note that we include div start and exclude div end since 
            // findDivArrayIndexFromDate() gets idx for div due to be 
            // paid strictly after the date given
            // now copy the divs to the result list
            for (idx = startIdx; idx < endIdx; idx++) {
                if ( divs[idx].getDivType() != Dividend::CONTINUOUS ) {
                    resultDivs->divArray->push_back(divs[idx]);
                }
            }
        }
    }

    return resultDivs.release();
}

/** return a DateTimeArray made up of all the pay dates in the dividend list */
DateTimeArrayConstSP DividendList::getPayDates()const
{
    DateTimeArraySP dates(new DateTimeArray(divArray->size()));
    const DividendArray& divs = *divArray;

    // copy the dates
    for (int idx = 0; idx < divs.size(); idx++)
    {
        (*dates)[idx] = divs[idx].getPayDate();
    }

    return dates;
}

/** return a DateTimeArray made up of all the pay dates in the dividend list */
DateTimeArrayConstSP DividendList::getExDivDates()const
{
    DateTimeArraySP dates(new DateTimeArray(divArray->size()));
    const DividendArray& divs = *divArray;

    // copy the dates
    for (int idx = 0; idx < divs.size(); idx++)
    {
        (*dates)[idx] = divs[idx].getExDate();
    }

    return dates;
}

/** return a doubleList made up of all the amounts in a div list */
    DoubleArrayConstSP  DividendList::getDivAmounts()const
{
    DoubleArraySP amounts(new DoubleArray(divArray->size()));
    const DividendArray& divs = *divArray;
    
    // copy the dividend amounts
    for (int idx = 0; idx < divs.size(); idx++)
    {
        (*amounts)[idx] = divs[idx].getDivAmount();
    }

    return amounts;
}

/** return a dividend list which is the result of merging the passed
    dividend list with 'this' dividend list. It is the callers 
    responsibility to free the memory allocated by the function. */
DividendList* DividendList::merge(const DividendList& otherList)const
{
    DividendListSP mergedDivs;
    const DividendArray& thisDivArray = *divArray;
    const DividendArray& otherDivArray = *otherList.divArray;

    // See if we can save ourselves some work ?

    // if the div list has requested to be merged with itself
    // or the passed div list is of zero length just return 'this'
    if (this == &otherList ||
        otherDivArray.size() == 0)
    {
        mergedDivs = DividendListSP(new DividendList(thisDivArray));
    }

    else if (thisDivArray.size() == 0)
    {
        mergedDivs = DividendListSP(new DividendList(otherDivArray));
    }

    else
    {
        // create an empty list and populate as we go along
        mergedDivs = DividendListSP(new DividendList());
        mergedDivs->divArray = DividendArraySP(new DividendArray(0));
        int idxThis;
        int idxPass;
        int idx;
        int numMergedDivs = thisDivArray.size() + otherDivArray.size();

        for (idx = 0, idxThis = 0, idxPass = 0; idx < numMergedDivs; idx++)
        {
            if (idxThis == thisDivArray.size())
            {
                mergedDivs->divArray->push_back(otherDivArray[idxPass]);
                idxPass++;
            }
            else if (idxPass == otherDivArray.size())
            {
                mergedDivs->divArray->push_back(thisDivArray[idxThis]);
                idxThis++;
            }
            else if(otherDivArray[idxPass].getExDate().isGreater(
                thisDivArray[idxThis].getExDate()))
            {
                mergedDivs->divArray->push_back(thisDivArray[idxThis]);
                idxThis++;
            }
            else
            {
                mergedDivs->divArray->push_back(otherDivArray[idxPass]);
                idxPass++;
            }
        }
    }
    
    return mergedDivs.release();
}

/** scales all dividends by multFact ie amount -> amount*multFact */
void DividendList::scale(double multFact){
    DividendArray& divs = *divArray;
    for (int i = 0; i < divs.size(); i++){
        divs[i].scale(multFact); 
    }
}

/** Converts any yield dividends to absolute amounts using fwd price of
    asset. Validates against continuous dividends */
bool DividendList::convertYieldDivs(const Asset* asset) {
    bool hasYield = false;
    DividendArray& divs = *divArray;
    if (!divs.empty()){
        static const string method("DividendList::convertYieldDivs");
        try{
            // optimise by only doing fwd price for actual yield divs
            DateTimeArray dates;
            for (int i = 0; i < divs.size(); i++){
                Dividend::TDividendType divType = divs[i].getDivType();
                const DateTime& exDivDate = divs[i].getExDate();
#if 0
                // Simply let these through now
                if (divType == Dividend::CONTINUOUS){
                    throw ModelException(method,
                                         "Continuous dividends "+
                                         exDivDate.toString()+
                                         " not supported");
                }
#endif
                if (divType == Dividend::PERCENT){
                    // note must use before ex-div to get right price
                    dates.push_back(DateTime(exDivDate.getDate(), 
                                             DateTime::BEFORE_EX_DIV_TIME));
                    hasYield = true;
                }
            }
            if (!dates.empty()){
                DoubleArray fwds(dates.size());
                asset->fwdValue(dates, fwds);
                // then do scaling
                int fwdIdx = 0;
                for (int j = 0; j < divs.size(); j++){
                    if (divs[j].getDivType() == Dividend::PERCENT){
                        divs[j].convertToDollar(fwds[fwdIdx]);
                        fwdIdx++;
                    }
                }
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }
    return hasYield;
}


/** Converts any yield dividends to absolute amounts using fwd price of
    asset. Validates against continuous dividends */
bool DividendList::convertDollarDivs(const Asset* asset) {
    bool hasDollar = false;
    DividendArray& divs = *divArray;
    if (!divs.empty()){
        static const string method("DividendList::convertDollarDivs");
        try{
            // optimise by only doing fwd price for actual yield divs
            DateTimeArray dates;
            for (int i = 0; i < divs.size(); i++){
                Dividend::TDividendType divType = divs[i].getDivType();
                const DateTime& exDivDate = divs[i].getExDate();
                if (divType == Dividend::AMOUNT){
                    // note must use before ex-div to get right price
                    dates.push_back(DateTime(exDivDate.getDate(), 
                                             DateTime::BEFORE_EX_DIV_TIME));
                    hasDollar = true;
                }
            }
            if (!dates.empty()){
                DoubleArray fwds(dates.size());
                asset->fwdValue(dates, fwds);
                // then do scaling
                int fwdIdx = 0;
                for (int j = 0; j < divs.size(); j++){
                    if (divs[j].getDivType() == Dividend::AMOUNT) {
                        divs[j].convertToYield(1.0 / fwds[fwdIdx]);
                        fwdIdx++;
                    }
                }
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }
    return hasDollar;
}


/** Adjusts each of the dividends using the adjustDividend() method
    on adjuster. */
void DividendList::adjust(IDivAdjuster* adjuster){
    DividendArray& divs = *divArray;
    for (int i = 0; i < divs.size(); i++){
        Dividend& div = divs[i];
        adjuster->adjustDividend(div);
    }
}

/** Scales any dollar dividends in the future (based upon pay
    date) by fwd rate whilst those in the past are scaled by the
    spot. Yield divs are scaled by fwd rate(pay date)/fwd rate(ex div date).
    Continuous divs are left alone */
void DividendList::makeStruck(const CAsset*    fxAsset,
                              const DateTime&  valDate){
    DividendArray& divs = *divArray;
    if (!divs.empty()){
        double spot = fxAsset->getSpot();
        // scale historic dates by fx spot. Note: even if ex-div dates in order
        // doesn't mean pay dates are in order
        for (int i = 0; i < divs.size(); i++){
            Dividend::TDividendType divType = divs[i].getDivType();
            if (divType == Dividend::AMOUNT || divType == Dividend::PERCENT){
                const DateTime& payDate = divs[i].getPayDate();
                if (payDate.isGreater(valDate)){
                    divs[i].scale(payDate.isGreater(valDate)? 
                                  fxAsset->fwdValue(payDate): spot);
                    if (divType == Dividend::PERCENT){
                        // yield divs in the past will fail here
                        divs[i].scale(1.0/
                                      fxAsset->fwdValue(divs[i].getExDate()));
                    }
                }
            }
        }
    }
}

/** return the dividend array  */
const DividendArray& DividendList::getArray()const
{
    return *divArray;
} 

/** Update known cashflows object with details of dividends (dollar divs
    are known, whilst yields aren't */
void DividendList::updateKnownCashFlows(
    OutputRequestUtil::KnownCashFlows& knownCFs,
    const string&                      isoCode,
    const DateTime&                    valueDate){
    DividendArray& divs = *divArray;
    for (int i = 0; i < divs.size(); i++){
        const DateTime& payDate = divs[i].getPayDate();
        const Dividend::TDividendType& divType = divs[i].getDivType();
        // dollar divs only 'known' after ex-div date 
        bool known = divType == Dividend::AMOUNT && 
            valueDate.isGreaterOrEqual(divs[i].getExDate());
        if (known){
            knownCFs.addKnownCashFlow(isoCode, 
                                      payDate, divs[i].getDivAmount());
        } else if (divType == Dividend::PERCENT || 
                   divType == Dividend::AMOUNT){
            knownCFs.addUnknownCashFlowDate(isoCode, payDate);
        } else {
            throw ModelException("DividendList::updateKnownCashFlows",
                                 "Continuous Dividends not supported");
        }
    }
}

/** Given a date, find the corresponding dividend with exDate
    immediately after that date.
    Returns 0 for the following cases:
    DividendList is empty
    theDate >= last exDate */  
DividendConstSP DividendList::getNextDivFromDate(const DateTime& date) const{
    int idx = findDivArrayIndexFromDate(date);
    if (idx < 0){
        return DividendConstSP();
    }
    return DividendConstSP(copy(&(*divArray)[idx]));
}

void DividendList::setElement(int idx, const Dividend& newDividend)
{
    if (idx >= divArray->size())
    {
        throw ModelException("DividendList::setElement", "Trying to set element "
            + Format::toString(idx) + " of dividend list with only "
            + Format::toString(divArray->size())
            + " elements.\n");
    }

    (*divArray)[idx] = newDividend;
    
}


// PRIVATE METHODS

/** Given a date, find the corresponding dividend idx with exDate
    immediately after that date , exDivDate[idx] STRICT > date.
    Returns -1 for the following cases:
    DividendList is empty
    theDate >= last exDate */ 
int DividendList::findDivArrayIndexFromDate(
    const DateTime&   theDate)const       // (I) the date
{
    const DividendArray& divs = *divArray;
    int numDivs = divs.size();
    int idx = -1;

    // Check if the date is within the dividend list range
    if (numDivs != 0                                    && 
        theDate.isLess(divs[numDivs-1].getExDate()) )
    {
        // check whether the first dividend is after theDate
        if ( theDate.isGreaterOrEqual(divs[0].getExDate()))
        {
            idx = 1;
            // loop through the dividends to find the idx
            while( idx < numDivs )
            {
                if ( theDate.isGreaterOrEqual(divs[idx-1].getExDate()) &&
                     theDate.isLess(divs[idx].getExDate())         )
                {
                    break;
                }
                ++idx;
            }
        }
        else
        {
            idx = 0;
        }
    }
    
    return idx;
}

/** Same as findDivArrayIndexFromDate except returns an iterator. 
    Returns divArray->end() for the following cases:
    DividendList is empty
    theDate < 1st exDate 
    theDate > last exDate */
vector<Dividend>::iterator DividendList::findDivArrayIteratorFromDate(
    const DateTime&   theDate)       // (I) the date
{
    int numDivs = divArray->size();
    int idx;
    vector<Dividend>::iterator iter = divArray->begin();
    DividendArray& divs = *divArray;

    // Check if the date is within the dividend list range
    if (numDivs == 0                                    ||
        theDate.isGreaterOrEqual(divs[numDivs-1].getExDate()))
    {
        iter = divs.end();
    }
    else
    {
        // check whether the first dividend is after theDate
        if ( theDate.isGreaterOrEqual(divs[0].getExDate()))
        {
            idx = 1;
            // loop through the dividends to find the idx
            while( idx < numDivs )
            {
                // let iter to the latter of the two dividends
                iter++; 
                if ( theDate.isGreaterOrEqual(divs[idx-1].getExDate()) &&
                     theDate.isLess(divs[idx].getExDate())         )
                {
                    break;
                }
                idx++;
            }
        }
    }

    return iter;
}

/** finds any continuous dividend paying at theDate - if there is one
    returns true and the amount else returns false and 0.0 */
bool DividendList::ctsDivPayingAtDate(
     const DateTime& theDate,
     double&         ctsDiv) const  // (O) ctsDiv (cts basis) paying at theDate
{
    const DividendArray& divs = *divArray;
    if (divs.size() == 0) {
        ctsDiv = 0.0;
    } else {
        for(int cIdx = divs.size()-1; cIdx>=0; cIdx--) {
            const Dividend& dividend = divs[cIdx];
            int   date = theDate.getDate();
            if (dividend.getDivType() == dividend.CONTINUOUS &&
                dividend.getExDate().getDate() <= date) 
            {
                double amount = dividend.getDivAmount();
                if (!Maths::isZero(amount)){
                    ctsDiv = log(1.+ amount);
                    return true;
                } else {
                    ctsDiv = 0.0;
                    return false;
                }
            }
        }
    }
    return false;
}


/** MU_PARALLEL shift method. Note that the Equity class inherits
    from MuParallel and this method is wrapped in the equity
    shift method */
bool DividendList::sensShift(MuParallel* shift, const DateTime& valueDate)
{
    static const string method = "DividendList::sensShift";
    double shiftSize = shift->getShiftSize();
    DividendArray& divs = *divArray;
    if (!(divs.size() == 0))
    {
        bool ctsDivAtStart = false;
        bool ctsDivAtEnd = false;
        double ctsDivStart = 0;
        double ctsDivEnd = 0;
        vector<Dividend>::iterator startIdx;
        vector<Dividend>::iterator endIdx;

        // see if there is a cts div in effect at start date 
        findCtsDivsForBucket(valueDate,
                             valueDate,
                             ctsDivAtStart,
                             ctsDivStart,
                             startIdx,
                             ctsDivAtEnd,
                             ctsDivEnd,
                             endIdx);

       if ( ctsDivAtStart )
       {
            // insert cts div at value date - continuous div today affects
            // tonight's growth - hence insert today 
            startIdx = dividendInsert(startIdx,
                                      Dividend::CONTINUOUS,
                                      valueDate,
                                      ctsDivStart);
        }
        
        if (startIdx != divs.end())
        {
            dividendTweakBetweenIndexes(startIdx,
                                        divs.end(),
                                        shiftSize);
        }
    }
    return false;
}

/** MU_S shift method. Note that the Equity class inherits
    from MuSpecial and this method is wrapped in the equity
    shift method */
bool DividendList::sensShift(MuSpecial* shift, 
                             double spotPrice,
                             const DateTime & expiry)
{
    const static string method = "DividendList::sensShift";
    DateTime firstContDivDate;
    double shiftSize = shift->getShiftSize();

    bool contDivFound = hasContinuousDividend(firstContDivDate);

    /* only insert dividend if start date is before first continous 
       dividend date. Note that the Equity::sensShift method which
       wraps this method has already checked that the expiry date
       is on or after today */
    if (!contDivFound || expiry.isLess(firstContDivDate))
    {
        // calculate absolute dividend amount
        double dividendAmount = spotPrice * shiftSize;

        // find the index of the dividend which we want to switch 
         vector<Dividend>::iterator insertionIdx = 
             findDivArrayIteratorFromDate(expiry);

         // insert a new dividend 
         dividendInsert(insertionIdx,
                        Dividend::AMOUNT,
                        expiry,
                        dividendAmount);
    }
    return false;
}

/** MU_POINTWISE shift method. Note that the Equity class inherits
    from MuPointwise and this method is wrapped in the equity
    shift method */
bool DividendList::sensShift(MuPointwise*    shift, 
                             const DateTime& valueDate,
                             DateTime&       bucketStartDate,
                             DateTime&       bucketEndDate)
{
    static const string method = "DividendList::sensShift";
    double shiftSize = shift->getShiftSize();

    if (bucketStartDate.isGreater(bucketEndDate))
    {
        throw ModelException(method, "bucket start date (" +
                             bucketStartDate.toString() +
                             ") is after the bucket end date (" +
                             bucketEndDate.toString() + ")");
    }

    /* NOTE - The price of some instruments are sensitive to the value of
       dividends even if the ex div date is in the past - so only
       tweak divs after today - except for cts divs which affect tonight's
       growth  */
    
    bucketStartDate = bucketStartDate.isGreater(valueDate)? 
        bucketStartDate : valueDate;
    bucketEndDate = bucketStartDate.isGreater(bucketEndDate)?
        bucketStartDate : bucketEndDate;

    /* if start date is value date then we place continuous dividends
       today (rather than tomorrow) as they effect tonight's growth - this
       is bit odd but is a consequence of a bucket and which divs are in a
       bucket */

    // if no dividends then dont do anything
    if (divArray->size() > 0)
    {
        bool ctsDivAtStart = false;
        bool ctsDivAtEnd = false;
        double ctsDivStart = 0;
        double ctsDivEnd = 0;
        vector<Dividend>::iterator divStartIdx;
        vector<Dividend>::iterator divEndIdx;
        
        /*  we may have to insert up to 2 extra dividends. Therefore we must
            reserve extra space on the divArray otherwise iterators
            pointing at elements on this array before insertion would become
            invalid after an insertion. */
        divArray->reserve(divArray->size() + 2);  // Note: size has not changed
            
        // see if there is a cts div in effect at start date 
        findCtsDivsForBucket(bucketStartDate,
                             bucketEndDate,
                             ctsDivAtStart,
                             ctsDivStart,
                             divStartIdx,
                             ctsDivAtEnd,
                             ctsDivEnd,
                             divEndIdx);
            
        if (ctsDivAtStart)
        {
            if (!(bucketStartDate.equals(valueDate)))
            {
                bucketStartDate = bucketStartDate.rollDate(1);
            }
            // insert cts div at start date + 1 (start date if on value date)
            divStartIdx = dividendInsert(divStartIdx,
                                         Dividend::CONTINUOUS,
                                         bucketStartDate,
                                         ctsDivStart);

            divEndIdx++;  // we've inserted a div 
        }

        if (ctsDivAtEnd)
        {
            // insert cts div at end date + 1
            divEndIdx = dividendInsert(divEndIdx,
                                       Dividend::CONTINUOUS,
                                       bucketEndDate.rollDate(1),
                                       ctsDivEnd);
        }
        
        // tweak dividends between the lower and upper indexes 
        dividendTweakBetweenIndexes(divStartIdx,
                                    divEndIdx,
                                    shiftSize);
    }
    return false;
}

/** given start and end date, report cts divs in effect at those time 
    plus index of appropriate dividend which marks the first div after 
    start date and first div after end date. Note index is array index, 
    also indexes can be equal to divList->fNumItems (ie one past last div).
    Validation only for start date and end date*/
void DividendList::findCtsDivsForBucket(
    const DateTime&  startDate,     
    const DateTime&  endDate,       
    bool&            ctsDivAtStart,  
    double&          ctsDivStart,    
    vector<Dividend>::iterator& divStartIdx,    
    bool&            ctsDivAtEnd,    
    double&          ctsDivEnd,      
    vector<Dividend>::iterator& divEndIdx)
{
    static const string method = "DividendList::findCtsDivsForBucket";

    if (startDate.isGreater(endDate))
    {
        throw ModelException(method,
                             "start date (" +
                             startDate.toString() +
                             ") is after the end date (" +
                             endDate.toString() +
                             ")");
    }
        
    // default values 
    ctsDivAtStart = false;
    ctsDivAtEnd   = false;

    for (divStartIdx  = divArray->begin(); divStartIdx !=divArray->end() &&
             startDate.isGreaterOrEqual(divStartIdx->getExDate());
         ++divStartIdx)
    {
        if (divStartIdx->getDivType() == Dividend::CONTINUOUS)
        {
            ctsDivStart = divStartIdx->getDivAmount();
            ctsDivAtStart = true;
        }
    }
    
    // set up data for end date 
    if (ctsDivAtStart)
    {
        ctsDivAtEnd = true;
        ctsDivEnd = ctsDivStart;
    }

    // now for end date - start from where we were 
    for (divEndIdx = divStartIdx; divEndIdx != divArray->end() &&
             endDate.isGreaterOrEqual(divEndIdx->getExDate());
         ++divEndIdx)
    {
        if (divEndIdx->getDivType() == Dividend::CONTINUOUS)
        {
            ctsDivEnd = divEndIdx->getDivAmount();
            ctsDivAtEnd = true;
        }
    }
    
    return;
}


/** tweak all divs for idx1 <= div < idx2 - note array indices. 
    No validation */
void DividendList::dividendTweakBetweenIndexes(
    vector<Dividend>::iterator  idx1,         
    vector<Dividend>::iterator  idx2,         
    double  shiftSize)
{
    for (vector<Dividend>::iterator  idx = idx1; 
         idx != idx2 && idx != divArray->end(); ++idx)
    {
        idx->tweakDividend(shiftSize);
    }
} 

/** insert a dividend at a given index, ex-div date = pay date */
vector<Dividend>::iterator DividendList::dividendInsert(
    vector<Dividend>::iterator idx,         
    Dividend::TDividendType  divType,      
    const  DateTime&   exDivDate,    
    double divAmount)
{
    Dividend newDiv(exDivDate, exDivDate, divType, divAmount);
    return divArray->insert(idx, newDiv);
}

/** Returns true if this DividendList contains any CONTINUOUS dividends,
    along with the date of the first continuous dividend */
bool DividendList::hasContinuousDividend(
    DateTime&    firstContDivDate)const
{
    bool contDivFound = false;
    const DividendArray& divs = *divArray;

    // cycle through each dividend checking the type 
    for (int idx = 0; idx < divs.size(); idx++)
    {
        if (divs[idx].getDivType() == Dividend::CONTINUOUS)
        {
            contDivFound = true;
            firstContDivDate = divs[idx].getExDate();
            break;
        }
    }
    
    return contDivFound;
}

/** Returns true if this DividendList contains any AMOUNT dividends */
bool DividendList::hasDollarDividend() const
{
    bool foundDollarDiv = false;
    const DividendArray& divs = *divArray;

    // cycle through each dividend checking the type 
    for (int idx = 0; idx < divs.size(); idx++)
    {
        if (divs[idx].getDivType() == Dividend::AMOUNT)
        {
            foundDollarDiv = true;
            break;
        }
    }
    
    return foundDollarDiv;
}

/** Returns true if this DividendList contains any AMOUNT dividends */
bool DividendList::hasYieldDividend() const
{
    bool foundYieldDiv = false;
    const DividendArray& divs = *divArray;

    // cycle through each dividend checking the type 
    for (int idx = 0; idx < divs.size(); idx++)
    {
        if (divs[idx].getDivType() == Dividend::PERCENT)
        {
            foundYieldDiv = true;
            break;
        }
    }
    
    return foundYieldDiv;
}

class DividendListHelper{
public:
    /* Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(DividendList, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDividendList);
        FIELD(divArray, "array of dividends");
    }

    static IObject* defaultDividendList() {
        return new DividendList();
    }
 };


CClassConstSP const DividendList::TYPE = CClass::registerClassLoadMethod(
    "DividendList", typeid(DividendList), DividendListHelper::load);

class DivListCreateAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // addin parameters
    DateTimeArraySP exDates;        // ex-dividend dates
    DateTimeArraySP payDates;       // payment dates
    CIntArraySP     divTypes;       // dividend types
    CDoubleArraySP  divAmounts;     // dividend amounts

    // create a dividend list
    static IObjectSP create(DivListCreateAddin *params)
        {
            static const string method = "DivListCreateAddin::create";

            // do some validation on the input parameters
            if ((params->exDates->size() != params->payDates->size())  ||
                (params->exDates->size() > params->divTypes->size())   ||
                (params->exDates->size() > params->divAmounts->size()) )
            {
                throw ModelException(method,
                                     "Mismatch between number of ex-div "
                                     "dates, pay dates, div types and "
                                     "div amounts.");
            }
            DividendArray array(0);
            DividendListSP newDivList(new DividendList(array));

            // loop through the paramters creating new dividends
            for (int i = 0; i < params->exDates->size(); i++)
            {
                Dividend tempDiv((*params->exDates)[i],
                                 (*params->payDates)[i],
                                 static_cast<Dividend::TDividendType>(
                                     (*params->divTypes)[i]),
                                 (*params->divAmounts)[i]);
                
                newDivList->divArray->push_back(tempDiv);
            }

            // check that we have created a decent dividend list
            newDivList->validatePop2Object();
            return newDivList;
        }


    DivListCreateAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DivListCreateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDivListCreateAddin);
        // order of registration effects order of parameters in addin function
        FIELD(exDates, "ex-div dates");
        FIELD(divAmounts, "div amounts");
        FIELD(divTypes, "div types");
        FIELD(payDates, "pay dates");
        Addin::registerClassObjectMethod("DIVIDEND_LIST",
                                         Addin::MARKET,
                                         "Creates a handle to a dividend "
                                         "list",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)create);
    }
    
    static IObject* defaultDivListCreateAddin(){
        return new DivListCreateAddin();
    }
 
};

CClassConstSP const DivListCreateAddin::TYPE = CClass::registerClassLoadMethod(
    "DivListCreateAddin", typeid(DivListCreateAddin), load);

class GetDivListAddin: public CObject {
public:
    static CClassConstSP const TYPE;

    // addin parameters
    IObjectSP divList;        // ex-dividend dates
    // create a dividend list
    static IObjectSP getDivList(GetDivListAddin *params) {
        static const string method = "GetDivListAddin::create";
        const DividendArray* divArrayPtr;
        if (DividendList::TYPE->isInstance(params->divList)){
            DividendList& divList = 
                dynamic_cast<DividendList&>(*params->divList);
            divArrayPtr = &divList.getArray();
        } else if (DividendArray::TYPE->isInstance(params->divList)){
            divArrayPtr = 
                dynamic_cast<const DividendArray*>(params->divList.get());
        } else {
            throw ModelException("GetDivListAddin", "divList must be of "
                                 "type DividendArray or DividendList");
        }
        const DividendArray& divList =  *divArrayPtr;

        smartPtr<DivListCreateAddin> outputDivList(new DivListCreateAddin());

        outputDivList->exDates    = DateTimeArraySP(new DateTimeArray(0));
        outputDivList->payDates   = DateTimeArraySP(new DateTimeArray(0));
        outputDivList->divTypes   = CIntArraySP    (new CIntArray(0)    );
        outputDivList->divAmounts = CDoubleArraySP (new CDoubleArray(0) );

        // loop through the dividend array creating new dividend lists
        for (int i = 0; i < divList.size(); i++)
        {
            outputDivList->exDates->push_back(   divList[i].getExDate()   );
            outputDivList->payDates->push_back(  divList[i].getPayDate()  );
            outputDivList->divAmounts->push_back(divList[i].getDivAmount());
            outputDivList->divTypes->push_back(  divList[i].getDivType()  );
        }

        return outputDivList;
    }

    GetDivListAddin(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetDivListAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetDivListAddin);
        // order of registration effects order of parameters in addin function
        FIELD(divList, "Dividend array");
        Addin::registerClassObjectMethod("GET_DIVIDEND_LIST",
                                         Addin::MARKET,
                                         "converts dividend list or array to"
                                         " a set of separate list",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)getDivList);
    }
    
    static IObject* defaultGetDivListAddin(){
        return new GetDivListAddin();
    }
};

CClassConstSP const GetDivListAddin::TYPE = CClass::registerClassLoadMethod(
    "GetDivListAddin", typeid(GetDivListAddin), load);

DEFINE_TEMPLATE_TYPE(DividendListArray);


DRLIB_END_NAMESPACE




