//----------------------------------------------------------------------------
//   Group       : QR&D Interest Rates
//
//   Filename    : RatesUtils.hpp
//   Description : rates utility functions
//
//   Author      : Steve Marks
//
//----------------------------------------------------------------------------

#ifndef HYMODELTURBO_HPP
#define HYMODELTURBO_HPP

#include "edginc/config.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/IrConverter.hpp"
#include "esl_types.h"
#include "esl_date.h"

DRLIB_BEGIN_NAMESPACE

class KBondSched;
class KFloatLegSched;
class KOptionSched;

class TREE_DLL RatesUtils {
    public:
    static MaturityPeriodSP getMaturityPeriod(ExpirySP expiry, DateTime refDate);
    static double calcDcf(DateTime start, DateTime end, string dcc);
    static double calcDcf(DateTime startDate, DateTime endDate, const DayCountConvention & dcc);
    static void   calcDcfFrac(DateTime startDate, DateTime endDate, const DayCountConvention & dcc, int &days /* OUT */, int &denom /* IN */);

    static void   genNewEvenDateList(DateTime startDate, DateTime matDate, int payFreq, bool stubAtEnd, DateTimeArray& list);
    static void   genDateAndCurveList(const DateTimeArray& dateListIn, const DoubleArray& curve,
                                     DateTimeArray& dateListOut, DoubleArray& curveOut);
    static void genDateList(DateTime startDate, DateTime matDate, int freq, bool stubAtEnd, DateTimeArray& dateList) ;
    static DateTime fwdByExpiry(DateTime dateIn, Expiry &expiry, bool forward);
    static DateTime nextMonth(DateTime date, int nbMonths, bool spillOver);

/* See comment in RatesUtils.cpp                                            */
#if 0
    static void parYield(
        double       &parYieldRate         /** (O) Par forward yield         */
        ,double      &annuity              /** (O) Annuity                   */
        ,YieldCurve  &yc                   /** (I) Yield Curve */
        ,DateTime    currentDate           /** (I) Current date              */
        ,DateTime    startDate             /** (I) Forward start date        */
        ,const MaturityPeriodSP &tenor     /** (I) Index maturity            */
        ,const DayCountConventionSP &dcc   /** (I) Index Day count convention*/
        ,const MaturityPeriodSP &frequency /** (I) Index payment frequency   */
        );
#endif
                                     
    static bool happensNow(const DateTimeArray& array, DateTime stepDate, int *pos=0);
    static void zbkDLFromIdx(const DateTimeArray &reset, const DateTimeArray &swapSt, 
        MaturityPeriodSP tenor, MaturityPeriodSP frequency, 
        DateTimeArray &idxZMat, DateTimeArray &idxZUse);

    static void zbkOptDates(
        vector<IRDate> const& matDL, /** (I) target mat date list        */
        vector<IRDate> const& useDL, /** (I) target use date list        */
        int nbCritDates,             /** (I) Nb of critical dates        */
        CRIT_DATE const *critDates,  /** (I) critical date structs       */
        DateTime valueDate,          /** (I) value date                  */
        vector<IRDate> &zbkMats,     /** (O) zerobank mats               */
        vector<IRDate> &zbkErs);     /** (O) zbank earliest use dates    */


    static void cbkProcessDL(
        vector<IRDate> &ev,
        vector<IRDate> &er);


    template<typename T, int U>
    struct LocalArray {
        int nb;
        T *array;
        IrConverter::AutoDeleteDRArray<T> arrayDel;

        LocalArray(vector<T> &init) : nb(init.size()), arrayDel(&array, U, &nb) {
            array = (T*) DR_Array (U, 0, init.size()-1);
            if (!array)
                throw ModelException(__FUNCTION__, IrConverter::slog.pop());
            for (size_t i=0; i<init.size(); ++i) {
                array[i] = init[i];
            }
        }
        LocalArray() : nb(0), array(0), arrayDel(&array, U, &nb) {}

        void set(vector<T> &result) {
            result.resize(nb);
            for (int i=0; i<nb; ++i) {
                result[i] = array[i];
            }
        }
        void append(vector<T> &result) {
            for (int i=0; i<nb; ++i) {
                result.push_back(array[i]);
            }
        }
    };

    static void AddDateToList(LocalArray<CRIT_DATE, CRITDATE> &list, DateTime value);

};

DRLIB_END_NAMESPACE
#endif
