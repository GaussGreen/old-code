#ifndef QLIB_IR_PRODUCTS_IRCONVERTER_HPP
#define QLIB_IR_PRODUCTS_IRCONVERTER_HPP

#include "esl_types.h"
#include "esl_market.h"
#include "irx/irxutils.h"
#include "fix123head.h"
#include <eslhead.h>
#include "cupslib.h"
#include "cupsmodl.h"
#include "edginc/RateTree.hpp"
#include "esl_log.h"

DRLIB_BEGIN_NAMESPACE

class MaturityPeriod;
class YieldCurve;
class Expiry;
class IRVolBase;
class IRVol;
class FXAsset;
class Equity;
class IRCalib;
class DayCountConvention;
class SimpleEquity;

class TREE_DLL IrConverter
{
    public:
    static StringLogger slog;
    static void checkSlog();

    static void DateTimeArrayToIr(const DateTimeArray &dates, vector<long> &irDates);
    static void IrToDateTimeArray(long *irDates, int nbIrDates, DateTimeArray &dates);
    static string StringFromIRDate(IRDate date);
    static const char* toIrIndexName(char freq);
    static void  to_T_CURVE(T_CURVE &tc, const YieldCurve *yc, bool growthCurve = false );
    static void  to_NEW_T_CURVE(T_CURVE &tc, const YieldCurve *yc, bool growthCurve = false );
    static double convExpiryToYears(const Expiry *expiry);
    static unsigned convExpiryToMonths(const Expiry *expiry);
    static void  to_SWAPVOL_DATA(SWAPVOL_DATA &sv, ExpiryArray &selectedTenors,
                                 ExpiryArray &selectedExpiries, const IRVolBase *vol);
    static void  to_BASEVOL_DATA(BASEVOL_DATA &bv, ExpiryArray &selectedTenors,
                                 ExpiryArray &selectedExpiries, DateTime& volBaseDate,
                                 const IRVolBase *vol, const string &calibIndex);
    static void  to_FXVOLATILITY_DATA(FXVOLATILITY_DATA &fv, const FXAsset *asset);
    static void  to_EQVOLATILITY_DATA(EQVOLATILITY_DATA &ev, const Equity *asset);
    static void  to_FXSMILE_DATA(FXSMILE_DATA &fs, const FXAsset *asset);
    static void  fill_FIX3_TREE_DATA(FIX3_TREE_DATA &tree_data);
    static void  to_MKTVOL_DATA(MKTVOL_DATA &mktvol_data, BASEVOL_EXPOSURE_DATA *selectedBV,
                                SWAPVOL_EXPOSURE_DATA *selectedSV, const string &calibIndex,
                                const T_CURVE &t_curve, const BASEVOL_DATA *bv,
                                const SWAPVOL_DATA *sv, bool smoothing);
    static void  to_MKTVOL_DATA(MKTVOL_DATA& mktvol_data, BASEVOL_EXPOSURE_DATA* selectedBV,
                                SWAPVOL_EXPOSURE_DATA* selectedSV,
                                const T_CURVE& t_curve, bool smoothing,
                                const IRVolBase* pVolBase,
                                IRVol::VolType calibVolType, IRVol::CalibType calibType,
                                ExpirySP calibExpiry, CDoubleSP calibVolOverride);
    static void to_MKTVOL_DATA_TMX_SMILE(MKTVOL_DATA& mktvol_data, RateTree::CcyIRParamsSP ccyIRParams);
    static void  to_MODELPARAMETERS_DATA(MODELPARAMETERS_DATA &mp, const string &smileStyle, const string &modelSet, IRCalib &irCalib);
    static void  to_MODELPARAMETERS_DATA(RateTree::IRModelParams &mp, const string &smileStyle, const string &modelSet, IRCalib &irCalib);
    static void  to_MODELPARAMETERS_DATA(RateTree::IRModelParams &mp, const string& smileKey, const string& modelKey, const string& engineKey, IRModelConfigTable& engineTable, IRExoticParamTable& smileTable, IRExoticParamTable& modelTable);
    static void  checkEnum(const string &oneCharEnum, const char *allowedValues, const char *varName);
    static const IRVolBase* getVol(const YieldCurve *yc, bool baseVol); // yc must be UntweakableYC
    static char  dccTo035A(const string &dcc);
    static char  dccTo035A(const DayCountConvention &dcc);
    static char  compoundBasisToCoS(int cb);
    static char  stubConvToBNS(const string &stub);
    static char  toEslFreq(int freq);
    static char  toEslFreq(Expiry *expiry);
    static int   eslToQlibFreq(char freq);
    static string dccFrom035A(char dcc);
    static void  to_EQ_DATA(EQ_DATA &ed, const SimpleEquity &se);
    static IrxTCalendar* makeCalendar(const Holiday *holidays) { return NULL; }
    static int findNamedElementIndex(const StringArray& elementArray, const string& elementName, const string &errMsg);
    static void q3GetZerosData(DateTimeArray &dates, DoubleArray &rates, const YieldCurve &yc);
    static void q3GetandExtendZerosData(DateTimeArray &dates, DoubleArray &rates, DateTime &spotDate, const YieldCurve &yc, const DateTime& today);



    /** This is the workhorse of the tree building method that is shared by
        all engines. See code for the description of the algorithm.
        @param vd           (I) value date
        @param cd           (I) vector of critical dates
        @param cm           (I) vector of critical zero maturity dates
        @param cu           (I) vector of critical zero use dates
        @param om           (I) vector of optional zero maturity dates
        @param ou           (I) vector of optional zero use dates
        @param zm           (O) combined critical/optional zero maturity dates
        @param zu           (O) combined critical/optional zero use dates
        @param optimize     (O) if false dates are simply combined

        @return             true on success, false on error
    */
    static bool  optimizeCriticalDates(
            IRDate vd,
            vector<IRDate> const& cd,
            vector<IRDate> const& cm, vector<IRDate> const& cu,
            vector<IRDate> const& om, vector<IRDate> const& ou,
            vector<IRDate>&       zm, vector<IRDate>&       zu,
            bool optimize = true);

    /******* AutoDeleteDRArray *******/

    template <typename T>
    class AutoDeleteDRArray {
        T **ptr;
        int type;
        int *nl;
        int nlOffset;
        int *nh;
        int nhOffset;
    public:
        AutoDeleteDRArray(T **ptr, int type, int *nbEl);
        AutoDeleteDRArray(T **ptr, int type, int *nl, int nlOffset, int *nh, int nhOffset);
        ~AutoDeleteDRArray();
    };
};

/******* template code *******/

template <typename T>
inline IrConverter::AutoDeleteDRArray<T>::AutoDeleteDRArray(
    T **ptr, int type, int *nbEl)
  : ptr(ptr), type(type), 
    nl(0), nlOffset(0), 
    nh(nbEl), nhOffset(-1) 
{
    ASSERT(ptr);    
}


template <typename T>
inline IrConverter::AutoDeleteDRArray<T>::AutoDeleteDRArray(
    T **ptr, int type, 
    int *nl, int nlOffset, 
    int *nh, int nhOffset)
  : ptr(ptr), type(type), 
    nl(nl), nlOffset(nlOffset), 
    nh(nh), nhOffset(nhOffset) 
{
    ASSERT(ptr);    
}


template <typename T>
IrConverter::AutoDeleteDRArray<T>::~AutoDeleteDRArray() {
    int nlLoc = (nl ? *nl : 0) + nlOffset;
    int nhLoc = (nh ? *nh : 0) + nhOffset;
    Free_DR_Array(*ptr, type, nlLoc, nhLoc);
    *ptr = NULL;
}

DRLIB_END_NAMESPACE

#endif
