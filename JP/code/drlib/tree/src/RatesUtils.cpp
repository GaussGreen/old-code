//----------------------------------------------------------------------------
//   Group       : QR Interest Rates
//
//   Filename    : RatesUtils.cpp
//   Description : rates utility functions
//
//   Author      : Steve Marks
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "esl_log.h"
#include "esl.h"

DRLIB_BEGIN_NAMESPACE


void RatesUtils::AddDateToList(LocalArray<CRIT_DATE, CRITDATE> &list, DateTime value) {
    if (Add_To_DateList(&list.nb, 
                        &list.array,
                        value.toIrDate(),
                        0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
    {
        throw ModelException(__FUNCTION__, IrConverter::slog.pop());
    }
}

void RatesUtils::cbkProcessDL(
    vector<IRDate> &ev,
    vector<IRDate> &er)
{
    LocalArray<long, LONG> evLoc(ev);
    LocalArray<long, LONG> erLoc(er);

    if (CbkProcessDL(&evLoc.nb,
                     &evLoc.array,
                     &erLoc.nb,
                     &erLoc.array) != SUCCESS)
    {
        throw ModelException(IrConverter::slog.pop());
    }

    evLoc.set(ev);
    erLoc.set(er);
}

void RatesUtils::zbkOptDates(
    vector<IRDate> const& matDL, /** (I) target mat date list        */
    vector<IRDate> const& useDL, /** (I) target use date list        */
    int nbCritDates,             /** (I) Nb of critical dates        */
    CRIT_DATE const *critDates,  /** (I) critical date structs       */
    DateTime valueDate,          /** (I) value date                  */
    vector<IRDate> &zbkMats,     /** (O) zerobank mats               */
    vector<IRDate> &zbkErs)      /** (O) zbank earliest use dates    */
{
   LocalArray<long, LONG> matDates;
   LocalArray<long, LONG> ersDates;

   if (ZbkOptDates(matDL.size(),
                    const_cast<IRDate*>(&matDL[0]),
                    useDL.size(),
                    const_cast<IRDate*>(&useDL[0]),
                    nbCritDates,
                    const_cast<CRIT_DATE*>(critDates),
                    valueDate.toIrDate(),
                    &(matDates.nb),
                    &(matDates.array),
                    &ersDates.nb,
                    &ersDates.array) != SUCCESS)
    {
        throw ModelException(IrConverter::slog.pop());
    }
    matDates.append(zbkMats);
    ersDates.append(zbkErs);
}


// computes the period between expiry and refDate. Useful if expiry is a BenchmarkDate (>refDate)
MaturityPeriodSP RatesUtils::getMaturityPeriod(ExpirySP expiry, DateTime refDate) {
    try {
        if (!expiry)
            return MaturityPeriodSP();

        if (dynamic_cast<MaturityPeriod*>(expiry.get()))
            return MaturityPeriodSP::dynamicCast(expiry);

        BenchmarkDate* bmd = dynamic_cast<BenchmarkDate*>(expiry.get());
        if (!bmd)
            throw ModelException("Unknown Expiry concrete type");

        return MaturityPeriodSP(MaturityPeriod::dateSubtract(bmd->toDate(), refDate));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void RatesUtils::zbkDLFromIdx(const DateTimeArray &reset, const DateTimeArray &swapSt, 
                             MaturityPeriodSP tenor, MaturityPeriodSP frequency, 
                             DateTimeArray &idxZMat, DateTimeArray &idxZUse) {

    try {
        vector<long> resetDLv, swapStDLv; 
        IrConverter::DateTimeArrayToIr(reset, resetDLv);
        IrConverter::DateTimeArrayToIr(swapSt, swapStDLv);
        
        long *idxZMatDL=NULL;
        long *idxZUseDL=NULL;
        int nbIdxZMat=0;
        int nbIdxZUse=0;
        IrConverter::AutoDeleteDRArray<long>
            idxZMatDLDelete(&idxZMatDL, LONG, &nbIdxZMat);
        IrConverter::AutoDeleteDRArray<long>
            idxZUseDLDelete(&idxZUseDL, LONG, &nbIdxZUse);


        if (ZbkDLFromIdx (
            reset.size(),
            &resetDLv[0],
            &swapStDLv[0],
            tenor->toMonths(),
            IrConverter::toEslFreq(frequency->annualFrequency()),
            &nbIdxZMat,
            &idxZMatDL,
            &nbIdxZUse,
            &idxZUseDL) != SUCCESS) 
        {
                throw ModelException(IrConverter::slog.pop());
        }
        IrConverter::IrToDateTimeArray(idxZMatDL, nbIdxZMat, idxZMat);
        IrConverter::IrToDateTimeArray(idxZUseDL, nbIdxZUse, idxZUse);
    } 
    catch(exception &e) {
        throw ModelException(e,__FUNCTION__);
    }
}

/* Jerry Cohen - It does not look like this is being used, and we are   *
**               phasing out the version of Par_Yield that takes        *
**               explicit dates and rates in favor of one that takes a  *
**               zero curve (T_CURVE*) directly.                        */
#if 0
void RatesUtils::parYield(
     double      &parYieldRate         /** (O) Par forward yield         */
    ,double      &annuity              /** (O) Annuity                   */
    ,YieldCurve  &yc                   /** (I) Yield Curve */
    ,DateTime    currentDate           /** (I) Current date              */
    ,DateTime    startDate             /** (I) Forward start date        */
    ,const MaturityPeriodSP &tenor     /** (I) Index maturity            */
    ,const DayCountConventionSP &dcc   /** (I) Index Day count convention*/
    ,const MaturityPeriodSP &frequency /** (I) Index payment frequency   */
    ) 
{
    long *dates=0;
    try {
        CashFlowArraySP ratesAndDates(yc.getRatesAndDates());

        DateTimeArray zeroDates(CashFlow::dates(*ratesAndDates));
        dates = new long[ratesAndDates->size()];
        for (int i=0; i < ratesAndDates->size(); ++i) {
            dates[i] = zeroDates[i].toIrDate();
        }

        if (Par_Yield( &parYieldRate,
                       &annuity,
                       ratesAndDates->size(),
                       &((*CashFlow::amounts(*ratesAndDates))[0]),
                       dates,
                       currentDate.toIrDate(),
                       startDate.toIrDate(),
                       IrConverter::convExpiryToMonths(tenor.get()),
                       IrConverter::dccTo035A(*dcc),
                       IrConverter::toEslFreq(frequency.get())
                       )!=SUCCESS)
        {
            throw ModelException(IrConverter::slog.pop());
        }
        if (dates) delete[] dates;
    }
    catch (exception& e) {
        if (dates) delete[] dates;
        throw ModelException(e, "RatesUtils::parYield");
    }
}
#endif

DateTime RatesUtils::nextMonth(DateTime date, int nbMonths, bool spillOver) {
    return DateTime::fromIrDate(
            Nxtmth( date.toIrDate(), nbMonths, spillOver)
    );
}

void RatesUtils::genDateList(DateTime startDate, DateTime matDate, int freq, 
    bool stubAtEnd, DateTimeArray& dateList) 
{
    long *dateListP =NULL;
    try {
        int numDates;

        /* Create a date list with the ZeroFreq */
        if (DateListFromFreq(startDate.toIrDate(),
                             matDate.toIrDate(),
                             IrConverter::toEslFreq(freq),
                             (stubAtEnd ? 'B' : 'F'),
                             &numDates,
                             &dateListP) != SUCCESS)
        {
            throw ModelException(IrConverter::slog.pop());
        }
        dateList.resize(numDates);
        for (int i=0; i<numDates; ++i) dateList[i]=DateTime::fromIrDate(dateListP[i]);
        free(dateListP);
    }
    catch (exception& e) {
        if (dateListP != NULL) free(dateListP);
        throw ModelException(e, "RatesUtils::genDateList2");
    }
}

void RatesUtils::genNewEvenDateList(DateTime startDate, DateTime matDate, int freq, bool stubAtEnd, DateTimeArray& dateList) {
    EVENT_LIST  *eventList = NULL; /* Event list for float pmts */
    try {
        long tempDates[2] = {startDate.toIrDate(), matDate.toIrDate()};
        eventList = DrNewEventListFromFreq(
            2,
           tempDates,
           IrConverter::toEslFreq(freq),
           (stubAtEnd ? 'B' : 'F'),
           'N',    /* 'Dates in' check not required */
           NULL, NULL, NULL, NULL, NULL);
        if (eventList == NULL) throw ModelException(IrConverter::slog.pop());

        dateList.resize(eventList->NbEntries);
        for (int i=0; i<eventList->NbEntries; ++i) {
            dateList[i]=DateTime::fromIrDate(eventList->Dates[i]);
        }
        DrFreeEventList(eventList);
    }
    catch (exception& e) {
        DrFreeEventList(eventList);
        throw ModelException(e, "RatesUtils::genDateList");
    }
}

void RatesUtils::genDateAndCurveList(const DateTimeArray& dateListIn, const DoubleArray& curve,
                                     DateTimeArray& dateListOut, DoubleArray& curveOut) {
    EVENT_LIST  *eventList = NULL;
    try {
        long dateL[MAXNBDATE];
        int i;

        for (i = 0; i < dateListIn.size(); i++)  dateL[i] = dateListIn[i].toIrDate();

        eventList = DrNewEventListFromFreq(
                             dateListIn.size(),
                             dateL,
                             'E', // review this opt->optionType[0],   // ??? make more robust
                             'N',            /* Stub not allowed            */
                             'Y',            /* Input dates must be in list */
                             const_cast<double*>(&curve[0]),
                             NULL, NULL, NULL, NULL);
        if (!eventList) throw ModelException(IrConverter::slog.pop());

        dateListOut.resize(eventList->NbEntries);
        curveOut.resize(eventList->NbEntries);

        for (i = 0; i < eventList->NbEntries; i++) {
            dateListOut[i] = DateTime::fromIrDate(eventList->Dates[i]);
            curveOut[i] = eventList->Curve[0][i];
        }

        DrFreeEventList(eventList);
    }
    catch (exception& e) {
        DrFreeEventList(eventList);
        throw ModelException(e, "RatesUtils::genDateList");
    }
}

DateTime RatesUtils::fwdByExpiry(DateTime dateIn, Expiry &expiry, bool forward) {
    try {
        MaturityPeriod *mp = dynamic_cast<MaturityPeriod *>(&expiry);
        if (!mp) throw ModelException("expiry must be of type MaturityPeriod");
        int count;
        string unit;
        IRDate dateOut;

        mp->decompose(count, unit);
        char charUnit = unit[0];
        if (charUnit=='Y') { 
            count *=12; 
            charUnit='M'; 
        }
        if (DrDateFwdAny (dateIn.toIrDate(),
                          count,
                          charUnit,
                          (forward?'F':'B'),
                          &dateOut) != SUCCESS)
        {
            throw ModelException(IrConverter::slog.pop());
        }
        return DateTime::fromIrDate(dateOut);
    }
    catch (exception& e) {
        throw ModelException(e, "RatesUtils::fwdByExpiry");
    }
}

double RatesUtils::calcDcf(DateTime startDate, DateTime endDate, const DayCountConvention & dcc) {
    double dcf;
    if (DrDayCountFraction (startDate.toIrDate(),
                            endDate.toIrDate(),
                            IrConverter::dccTo035A(dcc),
                            &dcf) != SUCCESS)
    {
        throw ModelException(IrConverter::slog.pop());
    }
    return dcf;
}

void RatesUtils::calcDcfFrac(DateTime startDate, DateTime endDate, const DayCountConvention & dcc, int &days /* OUT */, int &denom /* IN */) {
    days = dcc.days(startDate, endDate);
    double dcf = calcDcf(startDate, endDate, dcc);
    double denomDouble = days/dcf;
    denom = floor(denomDouble+.5);
    double diff = denomDouble-denom;
    if (!Maths::isZero(diff/(denomDouble*10 /* a bit of slack */))) {
        throw ModelException("Inconsistent dcf nominator and denominator");
    }
}


double RatesUtils::calcDcf(DateTime startDate, DateTime endDate, string dcc) {
    auto_ptr<DayCountConvention> dccObj(DayCountConventionFactory::make(dcc));
    return calcDcf(startDate, endDate, *dccObj);
}

// *** utils ////////////
bool RatesUtils::happensNow(const DateTimeArray& array, DateTime stepDate, int *pos) 
{
    for (int i=0; i < array.size(); ++i) {
        if (array[i] == stepDate) {
            if (pos) *pos = i;
            return true;
        }
    }
    return false;
}

bool RatesUtilsLoad(void) {
    return true;
}

DRLIB_END_NAMESPACE
