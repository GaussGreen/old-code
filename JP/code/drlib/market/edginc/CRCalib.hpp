//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CRCalib.hpp
//
//   Description : smile params for credit vol
//
//   Author      : Eva X Strasser 
//
//   Date        : 
//
//----------------------------------------------------------------------------

#ifndef _CRCALIB_HPP
#define _CRCALIB_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/ICDSSpotVol.hpp"
#include "edginc/VolRequest.hpp"

DRLIB_BEGIN_NAMESPACE

/** CR Vol Parameters consisting of smile info only 
    SmileInfo = qLeft, qRight, fwdShift
    ModelInfo = no alpha, only one beta (mean reversion) and no rho 
    => therefore only SmileRequest and no ModelRequest as in IRCalib */
class MARKET_DLL CRCalib: public MarketObject, 
               public virtual ICDSVol {
public:
    static CClassConstSP const TYPE;

    /* constructor from single set of params */
    CRCalib(const string& smileParamsStyle, double qLeft, double qRight, double fwdShift);

    /** What getProcessedVol() returns */
    class MARKET_DLL VolProcessed: public CObject,
                        public virtual IVolProcessed{
    public:
        static CClassConstSP const TYPE;
        virtual ~VolProcessed();

        /** identifies the market data name of the volatility */
        virtual string getName() const;
        /** calculates the trading time between two dates */
        virtual double calcTradingTime(const DateTime &date1, 
                                       const DateTime &date2) const;
        /** retrieve time measure for the vol */
        virtual TimeMetricConstSP GetTimeMetric()const;

        /** Returns the parameters as a double array. This is a bit poor */
        const DoubleArray& getParams() const;

        const double getQLeft() const;
        const double getQRight() const;
        const double getFwdShift() const;

        // Constructor
        VolProcessed(const string& name, const DoubleArray& params);
    private:
        static void load(CClassSP& clazz);
        // fields
        string      name; // $unregistered
        DoubleArray params; // $unregistered
    };

    /* How to ask for the smile parameters */
    class MARKET_DLL SmileRequest: public CVolRequest{
    public:
        static CClassConstSP const TYPE;
        virtual ~SmileRequest();
        SmileRequest(const string& style);
    private:
        static void load(CClassSP& clazz);
        friend class CRCalib;
        string style; // $unregistered
    };

    /** Returns name of vol */
    virtual string getName() const;

    /** Returns a CRCalib::VolProcessed containing smile parameters ICDSParSpreads* parameter is ignored */
    virtual CVolProcessed* getProcessedVol(const CVolRequest*       volRequest,
                                           const ICDSParSpreads*    cds) const;

    /** Fails. */
    /*virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;*/

    /** Abstract base class for old and new style way of specifying smile
        parameters */
    class MARKET_DLL Smile: public MarketObject {
        friend class CRCalib;
    public:
        virtual ~Smile();
        Smile(): MarketObject(TYPE){}
        
        /* constructor from single set of params */
        Smile(const string& smileParamsStyle, double qLeft, double qRight, double fwdShift);

        static CClassConstSP const TYPE;    

        string getName() const; // as per MarketObject (mandatory)
        void validatePop2Object(); // as per CObject (non-mandatory)
    
    protected:
        Smile(CClassConstSP clazz);
        Smile(const Smile& rhs);
        Smile& operator=(const Smile& rhs);
    
    private:
        static IObject* defaultConstructor();
        static void load(CClassSP& clazz);
        // fields
        string  name;
        StringArray style;
        DoubleArrayArray params;
    };
    
    typedef smartPtr<Smile> SmileSP;
    typedef smartConstPtr<Smile> SmileConstSP;
    // support for wrapper class
    typedef MarketWrapper<Smile> SmileWrapper;

      /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);
private:  
    friend class CRCalibHelper;
    CRCalib();
    CRCalib(const CRCalib& rhs);
    CRCalib& operator=(const CRCalib& rhs);
    // fields 
    string           name;    // name of this market data object
    SmileWrapper     smileCalib;
};

typedef smartPtr<CRCalib> CRCalibSP;
typedef smartConstPtr<CRCalib> CRCalibConstSP;

// support for wrapper class
typedef MarketWrapper<CRCalib> CRCalibWrapper;


DRLIB_END_NAMESPACE

#endif
