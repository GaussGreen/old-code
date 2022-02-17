//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IRModelVNFM.hpp
//
//   Description : VNFM - Vladimir's n-Factor Model
//
//   Author      : Anwar E Sidat
//
//   Date        : 15-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IRModelVNFM_HPP
#define QLIB_IRModelVNFM_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/IRExoticParam.hpp"

DRLIB_BEGIN_NAMESPACE

// To avoid redundant file includes.
FORWARD_DECLARE(Expiry);

/** VNFM Model Class.
 */
class MARKET_DLL IRModelVNFM : public IRExoticParam
{
public:
    static CClassConstSP const TYPE;

    IRModelVNFM();
    virtual ~IRModelVNFM();

    /** Returns name of model. */
    virtual string getName() const;

    /** Get methods for input parameters. */
    virtual int numFactors() const { return nbFactors; }
    virtual const ExpiryArray& dates() const { return Dates; }
    virtual const DoubleMatrix& meanReversion() const { return MeanReversion; }
    virtual const DoubleMatrix& weight() const { return Weight; }
    virtual const DoubleMatrix& correlation() const { return Correlation; }
    virtual double backbone() const { return Backbone; }
    
    /** overrides default */
    virtual void validatePop2Object();

protected:

    IRModelVNFM(const CClassConstSP& clazz);
    IRModelVNFM(const IRModelVNFM& irv);
    IRModelVNFM& operator=(const IRModelVNFM& irv);


    //Fields
    string        Key;           // handle name
    int           nbFactors;     // number of factors (number of cols on matrices below)
    ExpiryArray   Dates;         // array of dates
    DoubleMatrix  MeanReversion; // matrix of mean reversion (rows represent date, columns represent factor)
    DoubleMatrix  Weight;        // matrix of weights (rows represent date, columns represent factor)
    DoubleMatrix  Correlation;   // (optional) correlation matrix only required for 2 and 3 factors, (rows represent date, columns represent 1&2, 1&3 and 2&3 etc).
    double        Backbone;      // value for backbone

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new IRModelVNFM(); }
    
    int          i_numDates;    // number of dates (number of rows on matrices above)
};

typedef smartConstPtr<IRModelVNFM> IRModelVNFMConstSP;
typedef smartPtr<IRModelVNFM>      IRModelVNFMSP;
typedef MarketWrapper<IRModelVNFM> IRModelVNFMWrapper;

#ifndef QLIB_IRModelVNFM_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IRModelVNFM>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IRModelVNFM>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IRModelVNFM>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IRModelVNFM>);
#endif

// support for wrapper class
#ifndef QLIB_IRModelVNFM_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IRModelVNFM>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRModelVNFM>);
#endif

DRLIB_END_NAMESPACE

#endif
