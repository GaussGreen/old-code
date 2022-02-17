//----------------------------------------------------------------------------
//
//   Group       : Global QR
//
//   Filename    : IRVolPair.hpp
//
//   Description : holds swapVol and baseVol
//
//   Date        : 10 Jan 2006
//
//
//----------------------------------------------------------------------------

#ifndef _IRVOL_PAIR_HPP
#define _IRVOL_PAIR_HPP

#include "edginc/IRVol.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IRVolPair: public IRVolBase,
                 public virtual IVolatilityBS {
public:
	
    static CClassConstSP const TYPE;

    /** overrides default */
    virtual void validatePop2Object();

    /** Returns name of vol */
    virtual string getName() const;

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const MarketObject*  yc) const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    // set volToUse
    void setVolToUse(const string& volToUse);

protected:
    virtual ~IRVolPair() {}

    friend class IRVolPairHelper;
    IRVolPair();
    IRVolPair(const IRVolPair& rhs);
    IRVolPair& operator=(const IRVolPair& rhs);

    /************* exported fields *************/
    string          name;    // name of the vol
    MarketWrapper<IRVol>   baseVol;
    MarketWrapper<IRVol>   swapVol;
    /************* end of fields *************/

    // trnasient field
    string      volToUse;
private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
};

// ***** helper class just to hold data to be collected
// this is just for passing down to rates trees, to be reviewed
class MARKET_DLL IRVolRaw : public virtual IVolProcessed,
                 public CObject{

    IRVolRaw();
    static IObject* emptyOne();
    static void load(CClassSP& clazz);

public:
    static CClassConstSP const TYPE;
    IRVolRaw(IRVolConstSP baseV, IRVolConstSP swapV) : CObject(TYPE), baseVol(baseV), swapVol(swapV) {}
    virtual string getName() const {return  baseVol->getName();}
    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, const DateTime &date2) const {return 0.0;}
    /** retieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric()const {return TimeMetricConstSP(   );}

    IRVolConstSP getBaseVol()const {return baseVol;}
    IRVolConstSP getSwapVol() const {return swapVol;}

private:
    IRVolConstSP     baseVol; // $unregistered
    IRVolConstSP     swapVol; // $unregistered
};
typedef smartPtr<IRVolRaw> IRVolRawSP;
// ***** end helper class

DRLIB_END_NAMESPACE

#endif
