//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRCalib.hpp
//
//   Description : calibration params for interest rate vol
//
//   Author      : Andrew J Swain
//
//   Date        : 7 November 2001
//
//
//----------------------------------------------------------------------------

#ifndef _IRCALIB_HPP
#define _IRCALIB_HPP

#include "edginc/IRVolBase.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

/** IR Vol Parameters consisting of smile info together with the parameters for
    a 1, 2, or 3 factor IR model */
class MARKET_DLL IRCalib : public IRVolBase {
public:
    static CClassConstSP const TYPE;
    static CClassConstSP Smile2Q_TYPE();

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
        /** Return labels associated with the params array */
        const StringArray& getParamLabel() const;

        //// Constructor
        VolProcessed(const string& name, const DoubleArray& params);
        VolProcessed(const string& name, const DoubleArray& params, const StringArray& paramLabel);
    private:
        static void load(CClassSP& clazz);
        /// fields ///
        string      name; // $unregistered
        DoubleArray params; // $unregistered
        StringArray paramLabel; // $unregistered
    };

    /* How to ask for the smile parameters */
    class MARKET_DLL SmileRequest: public CVolRequest{
    public:
        static CClassConstSP const TYPE;
        virtual ~SmileRequest();
        SmileRequest(const string& style);
    private:
        static void load(CClassSP& clazz);
        friend class IRCalib;
        string style; // $unregistered
    };

    /* How to ask for the model parameters */
    class MARKET_DLL ModelRequest: public CVolRequest{
    public:
        static CClassConstSP const TYPE;
        virtual ~ModelRequest();
        ModelRequest(const string& style);
    private:
        static void load(CClassSP& clazz);
        friend class IRCalib;
        string style; // $unregistered
    };

    //// these are just the usual TYPE's for the nested classes. Doing this is
    //// easier than putting each one in the header file. 
    //// ie SMILE2Q_TYPE = Smile2Q::TYPE
    static CClassConstSP const SMILE2Q_TYPE;

    /** Returns the type of Model for the supplied number of factors */
    static CClassConstSP getModelType(int numFactors);

    /** Returns name of vol */
    virtual string getName() const;

    /** Returns a IRCalib::VolProcessed containing either smile or model
        parameters. */
    virtual IVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const MarketObject*  yc) const;

    /** Abstract base class for old and new style way of specifying smile
        parameters */
    class MARKET_DLL SmileBase: public MarketObject {
    public:
        static CClassConstSP const TYPE;
        virtual ~SmileBase();
    protected:
        SmileBase(CClassConstSP clazz);
    private:
        SmileBase(const SmileBase& rhs);
        SmileBase& operator=(const SmileBase& rhs);
        friend class IRCalib;
        static void load(CClassSP& clazz);
    };
    typedef smartPtr<SmileBase> SmileBaseSP;
    typedef smartConstPtr<SmileBase> SmileBaseConstSP;
    // support for wrapper class
    typedef MarketWrapper<SmileBase> SmileBaseWrapper;

    /** Abstract base class for old and new style way of specifying smile
        parameters */
    class MARKET_DLL Model: public MarketObject {
    public:
        static CClassConstSP const TYPE;
        virtual ~Model();
    protected:
        Model(CClassConstSP clazz);
    private:
        Model(const Model& rhs);
        Model& operator=(const Model& rhs);
        friend class IRCalib;
        static void load(CClassSP& clazz);
    };
    // support for wrapper class
    typedef MarketWrapper<Model> ModelWrapper;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);
private:
    class Vol; // to be removed - for backwards compatibility
    class Model1FL;
    class Model2FL;
    class Model3FL;
    class Smile; // to be removed - for backwards compatibility
    class Smile2Q;
    class ModelBase;
    friend class IRCalibHelper;
    IRCalib();
    IRCalib(const IRCalib& rhs);
    IRCalib& operator=(const IRCalib& rhs);
    /////// fields /////
    string           name;    // name of this market data object
    SmileBaseWrapper smileCalib;
    ModelWrapper     modelCalib;
};

typedef smartPtr<IRCalib> IRCalibSP;
typedef smartConstPtr<IRCalib> IRCalibConstSP;
#ifndef QLIB_IRCALIB_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IRCalib>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IRCalib>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IRCalib>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IRCalib>);
#endif

// support for wrapper class
typedef MarketWrapper<IRCalib> IRCalibWrapper;
#ifndef QLIB_IRCALIB_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IRCalib>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRCalib>);
#endif


DRLIB_END_NAMESPACE

#endif
