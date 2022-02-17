//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolProcessedStochGarf.cpp
//
//   Description : implementation for processed VolStochGarf vol object
//
//   Date        : 28 June 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolProcessedStochGarf.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const VolProcessedStochGarf::TYPE = CClass::registerClassLoadMethod(
    "VolProcessedStochGarf", typeid(VolProcessedStochGarf), load);

VolProcessedStochGarf::VolProcessedStochGarf(): CVolProcessedDVF(TYPE)
{
}

VolProcessedStochGarf::VolProcessedStochGarf(const VolStochGarfSP & vol): CVolProcessedDVF(TYPE) , myVol(vol)
{
}

// set the level of volatilty
void VolProcessedStochGarf::setScalar(double level)
{
    scalar = exp(level * myVol->volVol);
}


void VolProcessedStochGarf::load(CClassSP& clazz){
        REGISTER(VolProcessedStochGarf, clazz);
        SUPERCLASS(CVolProcessedDVF);
        EMPTY_SHELL_METHOD(defaultVolProcessedStochGarf);
        FIELD(myVol, "ptr to VolStochGarf input");
        FIELD_MAKE_TRANSIENT(myVol);
        FIELD(BMDates, "vol bench mark dates");
        FIELD_MAKE_TRANSIENT(BMDates);
}

VolProcessedStochGarf::~VolProcessedStochGarf()
{
}

/** calc expiry dates from bench mark labels */
void VolProcessedStochGarf::createExpiryDates()
{
    BMDates.resize(myVol->ATMVolBM->size());
    for (int i=0; i<BMDates.size(); i++)
    {// is this baseDate checked vs expiry baseDate ???
        BMDates[i] = (*(myVol->ATMVolBM))[i]->toDate(myVol->baseDate); 
    }
}

/** volatility calculator object using a cache of forward values */
class IVolCalculatorVolStochGarf : public CVolProcessedDVF::IVolCalculator {
public:
    /** initialise the caching of forward value storing the maturities, 
        the forwardmid and growthrate */
    IVolCalculatorVolStochGarf(const VolProcessedStochGarf* vol, 
                                const DateTimeArray&   Maturities // the time axis of the lattice above
                                ): vol(vol),
        maturity(Maturities.size())
        {
            static const string method("CVolProcessedStochGarf::VolCalculatorDVF");
            try {
                int iDate;                
                /* create a copy of the maturity dates */
                for(iDate=0; iDate<maturity.size(); iDate++) {
                    maturity[iDate] = Maturities[iDate];
                }                                                    
            }
            
            catch (exception& e) {
                throw ModelException(e, method);
            }
        }

    /** calculate local volatility */
    virtual void CalcLocVol(
            const CSliceDouble& strikes,        // Slice of strikes
            int                 iDate,          // index in the array of maturities
            CSliceDouble&       locVol)         // output
    {
        const static string method ="IVolCalculatorVolStochGarf::CalcLocVol";
        try {
            if (maturity.size()==0) {
                throw ModelException(method, "the cache for maturity has not been initialized");
            }
            else {
        
                const double* strike = &strikes[0];
                int sizeStrikeArray = strikes.size();
        
                vector<SImpV>  impV(sizeStrikeArray);
        
        
                locVol[0] = 0.0;
                for (int iStrike = 0; iStrike < sizeStrikeArray; ++iStrike) {
                    locVol[iStrike] = vol->CalcLocVol(maturity[iDate], strike[iStrike], true);
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    
private:
    const VolProcessedStochGarf* vol;
    DateTimeArray                maturity;
};

/** create a volatility calculator with a cache of forward values, storing the array of maturity dates */
CVolProcessedDVF::IVolCalculator* VolProcessedStochGarf::CreateVolCalculator(
    const DateTimeArray&   maturity,              // the time axis of the lattice
    bool                   isIntraDayInterp) const
{
    return new IVolCalculatorVolStochGarf(this, maturity);
}


// Main body of calculate local volatility
double VolProcessedStochGarf::CalcLocVol(const DateTime&  maturity,
                                         double    strike,
										 bool		isIntraDayInterp) const
{
    return myVol->computeLocVol(strike, maturity)*scalar;
}

// not possible to return correct one.  Need it for FD set up
double VolProcessedStochGarf::computeImpVol(const DateTime& maturity,
                                 double          strike) const { 
    return myVol->getATMVol(0.0); // just return ATMVol....
}

// calculate local volatility, note that jump contribution depends on DiffVolOnly flag
void VolProcessedStochGarf::CalcLocVol(CLatticeDouble*     strikes,
                               DateTimeArray&      maturities,
                               CLatticeDouble*	   locVol,
                               bool                isIntraDayInterp) const

{// to do:  This is not used at all .  Return has no meaning.
    const static string method ="VolProcessedStochGarf::CalcLocVol";
    throw ModelException(method, "CalcLocVol(CLatticeDouble....) is under construction");
}

// calculate local variance , note that jump contribution depends on DiffVarOnly flag
void VolProcessedStochGarf::CalcLocVar(CLatticeDouble*     strikes,
                               DateTimeArray&      maturities,
                               CLatticeDouble*	   locVar,
                               bool                isIntraDayInterp) const
{/*to do:  Nothing retuirn yet.
*/
}

// calculate local variance for one slice, diffusion scaled var supplied, tailored for Tree1fVolStochGarf
void VolProcessedStochGarf::CalcLocVar(const double* spots,
                               const vector<double>& diff_drift,
                               double step_cev_power,
                               double step_jump_rate,
                               double step_dt,
                               double step_diff_var,
                               int bot,
                               int top,
				               vector<double>& locVar) const // output
{
    const static string method ="VolProcessedStochGarf::CalcLocVar";

    // nothing to do
}


/** identifies the market data name of the volatility */
string VolProcessedStochGarf::getName() const 
{
    return myVol->getName();
}

/** calculates the trading time between two dates */
double VolProcessedStochGarf::calcTradingTime(const DateTime &date1, 
                               const DateTime &date2) const 
{
    return myVol->timeMetric->yearFrac(date1, date2);
}

/** retieve time measure for the vol */
TimeMetricConstSP VolProcessedStochGarf::GetTimeMetric() const
{
    return myVol->timeMetric;
}


DRLIB_END_NAMESPACE

