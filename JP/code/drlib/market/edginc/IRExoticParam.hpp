//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IRExoticParam.hpp
//
//   Description : Abstract smile base for IR models.
//
//   Author      : Anwar E Sidat
//
//   Date        : 18-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IREXOTICPARAM_HPP
#define QLIB_IREXOTICPARAM_HPP

#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/MarketTable.hpp"


DRLIB_BEGIN_NAMESPACE

// To avoid redundant file includes.
FORWARD_DECLARE(CDoubleMatrix);
FORWARD_DECLARE(BadDayConvention);
FORWARD_DECLARE(Holiday);
FORWARD_DECLARE(Expiry);
FORWARD_DECLARE(TimeMetric);
FORWARD_DECLARE(IRVol);
FORWARD_DECLARE(Surface);
FORWARD_DECLARE_WRAPPER(Holiday);

/** a base class (name holder) for ir exotic parameters
 */
class MARKET_DLL IRExoticParam : public MarketObject
{
public:
    static CClassConstSP const TYPE;

    virtual ~IRExoticParam();

protected:
    IRExoticParam(const CClassConstSP& clazz);
    IRExoticParam(const IRExoticParam& irv);
    IRExoticParam& operator=(const IRExoticParam& irv);

    /** Checks matrix sizes (it can be just a scalar). */
    static void validateMatrixSize(
        const CDoubleMatrix& m,
        int                  numRows,
        int                  numCols,
        const string&        strObjName,
        const string&        strSmileName);

    /** Creates a surface internally (used for MQ parameters). */
    static SurfaceSP makeSurface(
        const string&           surfaceName,
        const DoubleArraySP&    expiryArray,
        const DoubleArraySP&    tenorArray,
        const CDoubleMatrixSP&  vol);

    /** Creates a vol surface internally. */
    static IRVolSP makeVolSurface(
        const MarketData*       pMarket,
        const IModel*           pModel,
        const string&           volName,
        const TimeMetric*       timeMetric,
        const CDoubleMatrix&    vol,
        const ExpiryArray&      expiryArray,
        const ExpiryArray&      tenorArray,
        const DateTime&         baseDate,
        const HolidayWrapper    holidayWrapper,
        const BadDayConvention& badDayConv);

    /** Interpolates the smile distributions. */
    static double getData(
        IRVolSP         ptrVol,
        const DateTime& baseDate,
        const DateTime& expiryDate,
        const Expiry*   pMaturityTenor);

    /** Swap Maturity Tenor to year fraction */
    static double maturityToYears(const DateTime& baseDate, const Expiry& maturityTenor);

    /** Option Expiry Tenor to year fraction */
    static double expiryToYears(const DateTime& baseDate, const Expiry& expiryTenor);

    /** Expiry Date to year fraction */
    static double expiryDateToYears(const DateTime& baseDate, const DateTime& expiryDate);

private:
    static void load(CClassSP& clazz);
};


// Support for smart pointer, wrapper and market table
typedef smartConstPtr<IRExoticParam> IRExoticParamConstSP;
typedef smartPtr<IRExoticParam> IRExoticParamSP;
typedef array<IRExoticParamSP, IRExoticParam> IRExoticParamArray;
typedef MarketWrapper<IRExoticParam> IRExoticParamWrapper;
typedef smartPtr<IRExoticParamWrapper> IRExoticParamWrapperSP;
typedef array<IRExoticParamWrapperSP,IRExoticParamWrapper> IRExoticParamWrapperArray;
typedef MarketTable<IRExoticParam> IRExoticParamTable;
typedef MarketWrapper<IRExoticParamTable> IRExoticParamTableWrapper;
#ifndef QLIB_IREXOTICPARAM_CPP
    EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IRExoticParam>);
    EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IRExoticParam>);
    EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IRExoticParam>);
    EXTERN_TEMPLATE(class MARKET_DLL_SP MarketTable<IRExoticParam>);
    EXTERN_TEMPLATE(class MARKET_DLL_SP MarketWrapper<IRExoticParamTable>);
#else
    INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IRExoticParam>);
    INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IRExoticParam>);
    INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRExoticParam>);
    INSTANTIATE_TEMPLATE(class MARKET_DLL MarketTable<IRExoticParam>);
    INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRExoticParamTable>);
#endif

DRLIB_END_NAMESPACE

#endif
