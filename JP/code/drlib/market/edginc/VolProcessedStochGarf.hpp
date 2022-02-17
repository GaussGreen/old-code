//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedStochGarf.hpp
//
//   Description : interface for processed local vol object
//
//   Date        : 28 June 2005
//
//
//----------------------------------------------------------------------------

#ifndef VOL_PROCESSED_STOCH_GARF_HPP
#define VOL_PROCESSED_STOCH_GARF_HPP

#include "edginc/VolProcessed.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/VolStochGarf.hpp"

DRLIB_BEGIN_NAMESPACE

/** Type of CVolProcessedDVF. Returned following a LocalVolRequest */
class MARKET_DLL VolProcessedStochGarf: public CVolProcessedDVF{
public:
    static CClassConstSP const TYPE;
	static void load(CClassSP& clazz);
    static IObject* defaultVolProcessedStochGarf(){
        return new VolProcessedStochGarf();
    }

    /** create a volatility calculator with a caching of forward values
        not implemented here */
    virtual CVolProcessedDVF::IVolCalculator* CreateVolCalculator(
        const DateTimeArray&   maturity, // the time axis of the lattice
        bool                   isIntraDayInterp = true) const;
    
    /** constructor used but getProcessedVol() call in VolStochGarf */
    VolProcessedStochGarf(const VolStochGarfSP &);

    virtual ~VolProcessedStochGarf();

    // set the level of volatilty
    void setScalar(double level);
    
    // get garf function with different ATM vol level
//	IVolProcessed*	getVolProcessed(double scale);
    
    /** Calculates volatility */
    virtual void CalcLocVol(CLatticeDouble*     strikes,
                            DateTimeArray&      maturities, // the time axis of the lattice above
                            CLatticeDouble*	locVol,
                            bool		isIntraDayInterp=true) const;
    
    // not possible to return correct one.  Need it for FD set up9
    virtual double computeImpVol(const DateTime& maturity,
                                 double          strike) const;
    
    virtual double CalcLocVol(const DateTime& maturity,
                              double          strike,
                              bool	      isIntraDayInterp=true) const;
    
 
    /** Calculates variances */
    virtual void CalcLocVar(CLatticeDouble* 	strikes,
                            DateTimeArray& 	    maturities, // the time axis of the lattice above
                            CLatticeDouble*	    locVar,
                            bool			  isIntraDayInterp=true) const;
    
    void CalcLocVar(const DateTimeArray&  maturity,
                    const double*         strike,
                    const int             NStr,
                    vector<double>&	     locVar,
                    bool			  isIntraDayInterp=true) const{ }
    
    virtual string getName() const;
    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const;
    /** retieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric()const;

    /** calc expiry dates from bench mark labels */
    void createExpiryDates();

    // calculate local variance for one slice, diffusion scaled var supplied, tailored for Tree1fVolStochGarf
    void CalcLocVar(const double* spots,
                   const vector<double>& diff_drift,
                   double step_cev_power,
                   double step_jump_rate,
                   double step_dt,
                   double step_diff_var,
                   int bot,
                   int top,
                    vector<double>& locVar) const; // output
    
protected:
    VolProcessedStochGarf(const CClassConstSP& clazz);
    // all fields are transient
    VolStochGarfSP myVol;

private:
    friend class VolStochGarf;
    friend class CTree1fVolStochGarfCalib;
    VolProcessedStochGarf();
    VolProcessedStochGarf(const VolProcessedStochGarf &rhs);
    VolProcessedStochGarf& operator=(const VolProcessedStochGarf& rhs);

    // bench mark expiry dates
    DateTimeArray BMDates;
    double scalar;  // indicate volatility level; $unregistered

};

typedef smartConstPtr<VolProcessedStochGarf> VolProcessedStochGarfConstSP;
typedef smartPtr<VolProcessedStochGarf> VolProcessedStochGarfSP;

DRLIB_END_NAMESPACE
#endif
