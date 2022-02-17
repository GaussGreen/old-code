//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CEVProcessed.hpp
//
//   Description : interface for processed local vol object
//
//   Date        : 30 Oct 2001
//
//
//----------------------------------------------------------------------------

#ifndef CEV_PROCESSED_HPP
#define CEV_PROCESSED_HPP

#include "edginc/VolProcessed.hpp"
#include "edginc/VolProcessedDVF.hpp"

DRLIB_BEGIN_NAMESPACE

class CEVJ;

/** Type of CVolProcessedDVF. Returned following a LocalVolRequest */
class MARKET_DLL CEVJProcessed: public CVolProcessedDVF{
public:
    static CClassConstSP const TYPE;
	static void load(CClassSP& clazz);
    static IObject* defaultCEVJProcessed(){
        return new CEVJProcessed();
    }

    /** create a volatility calculator with a caching of forward values
        not implemented here */
    virtual CVolProcessedDVF::IVolCalculator* CreateVolCalculator(
        const DateTimeArray&   maturity, // the time axis of the lattice
        bool                   isIntraDayInterp = true) const;
    
    /** constructor used but getProcessedVol() call in CEVJ */
    CEVJProcessed(CEVJ*);

    virtual ~CEVJProcessed();

    /** decide if vol include jump or not */
    void SetVolCalc(bool diffVolOnly);
    
    /** Calculates volatility */
    virtual void CalcLocVol(CLatticeDouble*     strikes,
                            DateTimeArray&      maturities, // the time axis of the lattice above
                            CLatticeDouble*	locVol,
                            bool		isIntraDayInterp=true) const;
    
    virtual double computeImpVol(const DateTime& maturity,
                                 double          strike) const { return 0.0;}
    
    virtual double CalcLocVol(const DateTime& maturity,
                              double          strike,
                              bool	      isIntraDayInterp=true) const { return 0.0;}
    
 

    /** decide if var include jump or not */
    void SetVarCalc(bool diffVarOnly);
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

    // calculate local variance for one slice, diffusion scaled var supplied, tailored for Tree1fCEVJ
    void CalcLocVar(const double* spots,
                   const vector<double>& diff_drift,
                   double step_cev_power,
                   double step_jump_rate,
                   double step_dt,
                   double step_diff_var,
                   int bot,
                   int top,
                    vector<double>& locVar) const; // output


    void CalcJumpDrift(const double* spots,
                                      const DateTime& thisStepDate,
                                      vector<double>& jump_drift,
                                      vector<double>& d_node,
                                      vector<double>& u_node,
                                      int bot, int top) const;
    
    double CalcStepJumpRate(const DateTime& thisStepDate, const DateTime& nextStepDate) const;
    
    double InterpLinear(double x) const;

    double InterpLinear(const DateTime date, const DoubleArray valueArray) const;
    
    DoubleArray GetParamArr(const string& param_name) const;

    double GetSpotRef();

    DateTimeArray GetBMDates() const;

protected:
    CEVJProcessed(const CClassConstSP& clazz);
    // all fields are transient
    CEVJ* VolCEVJ;

private:
    friend class CEVJ;
    friend class CTree1fCEVJCalib;
    CEVJProcessed();
    CEVJProcessed(const CEVJProcessed &rhs);
    CEVJProcessed& operator=(const CEVJProcessed& rhs);

    // bench mark expiry dates
    DateTimeArray BMDates;

    // flags to indicate if jump is included in vol/var calculation
    bool DiffVolOnly;
    bool DiffVarOnly;

    vector<double>      X; // measured from RefDate if used for term structure $unregistered
    vector<double>      Y; // $unregistered
};

typedef smartConstPtr<CEVJProcessed> CEVJProcessedConstSP;
typedef smartPtr<CEVJProcessed> CEVJProcessedSP;

DRLIB_END_NAMESPACE
#endif
