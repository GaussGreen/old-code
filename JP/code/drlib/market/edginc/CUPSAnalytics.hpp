//----------------------------------------------------------------------------
//
//   Group       : DR
//
//   Filename    : CUPSAnalytics.hpp
//
//   Description : Code from CMLib library initially used for quanto adjustment
//                 to CDS. This file needs revisiting since we need to review
//                 how to bootstrap the curve and we should split the files into
//                 its component classes and put them in files with the right
//                 name.
//
//   Author      : CMLib (From file with samename)
//
//   Date        : 10 Dec 2004
//
//

#ifndef EDR_CUPSANALYTICS_HPP
#define EDR_CUPSANALYTICS_HPP

DRLIB_BEGIN_NAMESPACE

// The state of par CDS pricing/bootstrapping functions
// I1 * Spar = (1-R) * I2
// DF is risky spot discount factor
// "index" is the index of current benchmark maturity on the timeline
struct MARKET_DLL CDSFunctionState {
    double I1;
    double I2;
    double DF; // risky discount factor on interval [t0, t]
    int    index;

    CDSFunctionState();
};

// interface to par CDS pricing/bootstrapping functions
class MARKET_DLL BaseCDSFunction {
public:
    virtual ~BaseCDSFunction();
    // computes CDS value from forward spread
    // this function is used in bootstrapping 
    // it is called by root solver
    // lambda is the instantaneous hazard rate, 
    // which is assumed to be constant between 2 CDS benchmark dates 
    // (flat forward assumption)
    double operator() ( double lambda ) {
        m_state = m_initialState;
        updateState( lambda );
        return m_parSpread * m_state.I1 - (1-m_recoveryRate) * m_state.I2;
    }

    /** The zbrent solver seems stuck in C... */ 
    static double solverFunction(double lambda, void* data){
        BaseCDSFunction* cdsFunction = (BaseCDSFunction*)data;
        return (*cdsFunction)(lambda);
    }

    double parSpread( double lambda ) {
        m_state = m_initialState;
        updateState( lambda );
        return (1-m_recoveryRate)*m_state.I2/m_state.I1;
    }

    // specifies maturity and optionally par spread for next interval
    // should be called before using root solver or calling parSpread
    void nextBenchmark( int maturityIndex, double parSpread = 0 ) {
        m_initialState = m_state;
        m_benchmarkMaturity = maturityIndex + 1;
        m_parSpread = parSpread;
    }
protected:
    BaseCDSFunction( double recRate = 0 ) { 
        m_recoveryRate = recRate;
    }
    
    virtual void updateState( double lambda ) = 0;
    CDSFunctionState m_state, m_initialState;
    int m_benchmarkMaturity;
    double m_recoveryRate;
    double m_parSpread;
};

class MARKET_DLL ICovarianceFunctions {
public:
    virtual ~ICovarianceFunctions();

    // returns covariance C^r_(t,t)
    // where t is specified by maturity
    // (see CDS convexity adjustment and calibration paper)
    // http://stowe.ny.jpmorgan.com:8000/docdb/tmp/cdscalibadjust.pdf
    virtual double SpreadIrCov( int maturity ) = 0;

    // returns covariance C^X_(t,t)
    // where t is specified by maturity
    // (see CDS convexity adjustment and calibration paper)
    // http://stowe.ny.jpmorgan.com:8000/docdb/tmp/cdscalibadjust.pdf
    virtual double SpreadFxCov(const DoubleArray&    fxVol,
                               double                spreadFXCorr,
                               int                   maturity ) = 0;
};

// This is the interface to CUPS adjustment analytics
// It computes spread/IR covariance, spread/FX covariance and
// prices CDS with CUPS adjustment
class MARKET_DLL CUPSAnalyticsBase : public BaseCDSFunction,
                          public virtual ICovarianceFunctions { 
protected:
    DoubleArray    m_tOffset;
    double         m_betaR; // interest rate mean reversion
    double         m_betaL; // spread mean reversion
    DoubleArray    m_rBpVol; // interest rate basis point vol
    DoubleArray    m_lBpVol; // spread basis point vol
    double         m_spreadIRCorr;
    DoubleArray    m_fwdRlDF;// forward riskless discount factor
    DoubleArray    m_rRate;  // forward interest rate

    virtual void finalize() = 0;
public:
    
    //// This function selects CUPS adjustment analytics.
    //// There are different versions of analytics for cases of 
    //// zero mean reversion
    static CUPSAnalyticsBase* create(double betaR, double betaL);

    void initialize(
        const DoubleArray&    tOffset,
        const DoubleArray&    fwdRlDF, // forward riskless discount factor
        const DoubleArray&    rRate,   // forward interest rate
        const DoubleArray&    rBpVol,  // interest rate basis point vol
        const DoubleArray&    lBpVol,  // spread basis point vol
        double                spreadIRCorr,
        double                recoveryRate = 0);
    
};

// Subclass of BaseCDSFunction that implements standard CDS pricing
class MARKET_DLL UnadjustedCDSFunction : public BaseCDSFunction {
public:
    UnadjustedCDSFunction(
        const DoubleArray&    tOffset, // time offset
        const DoubleArray&    fwdRlDF, // forward riskless discount factor
        const DoubleArray&    rRate,    // forward interest rate
        double recoveryRate = 0 );
private:
    void updateState( double lambda );

    DoubleArray    m_tOffset;
    DoubleArray    m_fwdRlDF;
    DoubleArray    m_rRate;
};

DRLIB_END_NAMESPACE
#endif
