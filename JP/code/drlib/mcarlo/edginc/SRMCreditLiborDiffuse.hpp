//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditLiborDiffuse.hpp
//
//   Description : CreditLibor path generation
//
//
//----------------------------------------------------------------------------

#ifndef SRMCRLIBORDIFFUSE_HPP
#define SRMCRLIBORDIFFUSE_HPP

#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/SRMCreditDiffuse.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/DECLARE.hpp"
#include <set>

#include "edginc/SRMCreditLiborUtil.hpp"
#include "edginc/SRMCreditUtil.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"

DRLIB_BEGIN_NAMESPACE

class SRMCreditLiborDiffuse: public SRMCreditDiffuse 
{
public:

	// constructor
    SRMCreditLiborDiffuse(IQMCDiffusibleInterestRateSP srmRatesDiffuse):SRMCreditDiffuse(srmRatesDiffuse),srmCreditLiborUtil(NULL){}

	// populates class data
    void setSRMCreditLiborDiffuse(
        int                    randomIndex,
        const DateTime&        today,
        SRMCreditLiborUtilSP   srmCreditLiborUtil,
        double                 crFxCorr,
        const vector<double>&  prob); // historic dates

	// sets CR-IR correlation
	//virtual void setCrIrCorr(const vector<double>& corrCRIR);

	// finalazes the timeline, allocates memory
	virtual void finalizePathGenerator(DateTimeArrayConstSP allDates);

	// generate path across all dates
    virtual void generatePathAndSurvivalRates(      
        IQMCRNGManagerSP rngMgr);

	// contains the state variables of the LMM model
	struct DiffusedState
	{
		int             frozenIntensityIdx;			// index of the most recent frozen intensity
		int             firstIntAliveIdx;	    	// fist discrete intensity still alive at event date i
		double          frozenIntensity;			// most recent frozen intensity
		vector<double>  simulIntensity;	            // discrete intensities recorded at event dates
		vector<double>  survProb;	            // forward survival probabilities Q(t,T_k) recorded at event dates t for every k
	};

    typedef vector<DiffusedState>::const_iterator DiffusedStateIter;


	 /** accessing the expected value ExpSDF(md, fd) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
    virtual double getExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j);

    /** accessing the natural log of the expected value ExpSDF(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j);

protected:

	// populates esdfForwardYearFracs and expHesdfForwardYearFracs
    void storeT(const DateTime& today, const DateTimeArray& dates);


private:
	// model specific members
	int mNumInt;								// number of discrete intensities :=H_k, i.e. state variables

	double mYearFracToFirstReset;               // time in year frac from today to the first reset date

	vector<double> mInitialIntensities,			// H_k[0], for all k;
                   mResetDates,                 // := T_1, T-2,...T_n; the last intensity to evolve refers to the period [T_n, T_{n+1}]
				   mAccruedFrac,    		    	// := (0,T_1), (T_1, T_2),..., (T_n, T_{n+1)}
				   mSpotSurvProb,               // Q(0,T_j) for j=1,...n+1 ?? DO WE NEED IT????
				   mIntensityVols;				// volatilies to the state variables H_k

	// member common to other models
    //double rho03;								// correlation between CR and IR


	vector<double>  esdfForwardYearFracs,       // year fracs from today to measurement dates
                 	sqrtYearFrac;               // year fracs between simulation dates

	vector<DiffusedState>   expProb;            // for computing expected survival probabilit stores surv prob
												// length is equal to required nb where expProb has to be stored

   SRMCreditLiborUtilSP srmCreditLiborUtil;	    // util object
};

//declare smartPtr versions
DECLARE(SRMCreditLiborDiffuse);

DRLIB_END_NAMESPACE

#endif // SRMCRLIBORDIFFUSE_HPP
