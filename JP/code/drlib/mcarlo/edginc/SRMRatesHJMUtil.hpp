//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMRatesHJMUtil.hpp (from SRMRatesUtil.hpp)
//
//   Description : base class for SRM IR util classes (calibration etc)
//
//   Date        : 14 June 2004
//
//
//----------------------------------------------------------------------------

#ifndef SRMRATESHJMUTIL_HPP
#define SRMRATESHJMUTIL_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"


#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMRatesFactor.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE
class FXAsset;
DECLARE(FXAsset);

/** Helper class for models to get data out of market cache */        
class MCARLO_DLL SRMRatesHJMUtil : public SRMRatesUtil
{

public:
    
    ///// constructs and populates SRMRatesUtil object
    SRMRatesHJMUtil(
		    const DateTime&      baseDate,
            int                  numFactors,
            const string&        modelParams,
            const string&        smileParams,
            CVolProcessedSP      processedVol,
            IYieldCurveConstSP   discYC,
            IYieldCurveConstSP   diffYC,
            bool                 skipFlag,
            double               flatVolIr,
            const string&        cutoffChoice,
            double               constCutoffValue,
            const string&        corrSwapStart, // eg 1Y  (offset to today)
            const string&        corrSwapMat,   // eg 10Y (offset to start)
            const string&        corrSwapDCC,   // eg Act/365F
            const string&        corrSwapFreq); // eg 6M
           
    /** Constructs and populates SRMRatesUtil object suitable using the 'VF' style
        calibration. qLeft and qRight are hard coded to 1, as is alpha. Note
        that only single factor IR is supported */
    SRMRatesHJMUtil(
        const DateTime&      baseDate,
        FXAssetConstSP       fx, /* optional - if specified will add FX smile
                                    dates to timeline */
        const string&        modelParamsKey,
        VolProcessedBSIRSP   processedVol,
        IYieldCurveConstSP   discYC,
        IYieldCurveConstSP   diffYC,
        bool                 skipFlag,
        const string&        corrSwapStart, // eg 1Y  (offset to today)
        const string&        corrSwapMat,   // eg 10Y (offset to start)
        const string&        corrSwapDCC,   // eg Act/365F
        const string&        corrSwapFreq); // eg 6M

    virtual ~SRMRatesHJMUtil(){}

    virtual int numFactors() const;

    // called when we know simDates
    virtual void setTimeLine(DateTimeArrayConstSP simDates);
    /** this function returns product of spot vols and forward rates
        "interest rate basis point vol". The array returned is of the same
        length as the supplied array */
    DoubleArray basisPointVol(
        const DateTimeArray& dates) const; /* excludes today */

	/** Returns rho[rhoIndex] where rho are model parameters */
    virtual double getRho(int rhoIndex) const;  
    /** Returns the model rho parameters - length corresponds to off diagonal
        elements of a symmetric matrix (size of which is number of factors) */
    virtual const vector<double>& getRho() const; 
    // must keep in general interface, for now
    /** vol param */
    virtual double getAlpha(int factorIdx) const; 
    virtual void getAlpha(vector<double>& alphaByFactor) const;  // TODO: Move to SRMRatesHJMUtil?
    virtual void getBeta(vector<double>& betaByFactor) const;
    /** Returns beta for specified factor */
    virtual double getBeta(int factorIdx) const;
    /** Calls bFactor for the swap used for correlation purposes */
    virtual vector<double> bFactor() const;
    /** Calls bFactor for the swap given by the specified index */
    virtual vector<double> bFactor(int swaptionIndex) const;
    /*  Determine the value of the B coefficient (see Vladimir's memo)
        this function is Flat Forward ready */
    virtual vector<double> bFactor(const DateTime&      swapStart,
                           const DateTime&      swapMat,
                           DayCountConventionSP swapDayCC,
                           MaturityPeriodSP     swapFreq) const;
    // Calculates 'A factor' at next date using 'A factor' from previous date
    // and supplied parameters 
    // Note that fwdRateDt = fwdRate * deltaTime
    virtual void aFactor(double          fwdRateDt,
                         double          dt,
                         vector<double>& a) const; // (M)

	/** vol param */
    double getQRight() const{ return qRight; }
    /** vol param */
    double getQLeft() const{ return qLeft; }
    /** vol param */
    double getFwdShift() const{ return fwdShift; }

    /* Exponential decay function for specified factor */
    virtual double expDecay(int factor, double t) const;

    /* Produces the "usual" Lower triangular matrix for orthogonalising
       corrolated factors.  
       NOTE: If Nbfac < 3 then unused matrix elements are set to zero. */
    virtual DoubleMatrix triangulation() const;

    /*****  From util_s::Get3A **********************************************/
    /*
     *       Calculate the decomposition of the Correlation in a 2,3-factor case
     *       given the swap start date, a tenor and a frequency
     */
    virtual DoubleMatrix get3A() const;

	// returns the instantaneous vol by factor in the form (factor, time point)
	virtual void instFactorVol(
		vector< vector<double> >& vol,				  
		const vector<double>& DeltaTime,  // (I)  passed for convenience 
		const DateTimeArray& TPDate,      // (I)
        const vector<double>& IrFwdRate,  // (I)  pre-computed for performance
		int expiryIndex,                  // (I)  index in TPDate
		int fwdMatIndex) const;           // (I)  in case, the underlying has mat > expiry

    // Calculates a 'factor variance' as used by equity/FX vol calibration 
	virtual void irAssetVariance(
		vector<double>& irVariance,        // (O)  variance over [T(i-1),  T(i)] for all i <= N
		vector<double>& irAssetCovar,      // (O)  covariance in [T(i-1),  T(i)] for all i <= N
		const vector<double>& rhoEqIr,     // (I)
		const vector<double>& DeltaTime,   // (I)  passed for convenience 
		const DateTimeArray& TPDate,       // (I)  needed for simpleFwdCurve
        const vector<double>& IrFwdRate,   // (I)  pre-computed for performance
		const vector<int>& periodIndex,    // (I)  indices in TPDate
		int fwdMatIndex) const;            // (I)  indices in TPDate

    /** Returns 'variance' as determined by model parameters */
    virtual double modelVariance(double delt_1, double delt_2) const;

    /** Returns 'mean reversion integral' as determined by model parameters */
    virtual void modelMeanReversionIntegral(
        double          delt_1, 
        double          delt_2,
        vector<double>& meanReversionIntegrals) const; //O

    /** Utility routine to calculate the cutoff rate and populate
        the cutoff rate arrays */
    virtual void calcEffRateLimit(
        double          NbSigmasMax, // (I) Number or sigmas to cut at
        double          NbSigmasMin, // (I) Number or sigmas to cut at
        vector<double>& MaxRate,     /* (O) Cutoff forward rates  */
        vector<double>& MinRate) const;  /* (O) Cutoff forward rates      */

    /*  Calibration routine : populates ir->SpotVol and ir->tau array */
    virtual void spotVolBM(bool skipFlag);

    /** Calculates a 'factor variance' as used by FX vol calibration */
    virtual vector<double>::iterator factorVariance(
        vector<double>::iterator integral,  // (M)
        const vector<double>&    aFactor,
        double                   spotVol,
        double                   deltaTime) const;
	
    /** Calculates a 'factor covariance' as used by FX vol calibration */
    virtual vector<double>::iterator factorCovariance(
        vector<double>::iterator integral,  // (M)
        const vector<double>&    aFactor,
        double                   spotVol,
        double                   deltaTime) const;

    /** Calculates a 'factor covariance' as used by FX vol calibration */
    virtual vector<double>::iterator factorFXCovariance(
        bool                     isDomestic,
        vector<double>::iterator integral,  // (M)
        const vector<double>&    aFactor,
        const vector<double>&    rhoFxIRFac, /* (I) correl IR/FX */
        double                   spotVol,
        double                   deltaTime) const;
       
    void computeKFactor(vector<double>& kFactor, int factorIdx) const;
    void computeGFactor(vector<double>& gFactor, int factorIdx) const;
    
	double rBar(const DateTime& date) const;
    double rBar(int futureDatesIndex) const;

    /** computes KFactor between the 2 dates supplied. From irdiffuse::Kfactor */
    void kFactor(const DateTime& dateFrom,
                 const DateTime& dateTo,
                 vector<double>& k) const;
    /** computes KFactor between the 2 dates supplied  for specified model
        factor index */
    double kFactor(const DateTime& dateFrom, 
                   const DateTime& dateTo,
                   int   factorIdx) const;
    /** computes GFactor between the 2 dates supplied. From irdiffuse::Gfactor */
    void gFactor(const DateTime& dateFrom,
                 const DateTime& dateTo,
                 vector<double>& g) const;
    /** computes GFactor between the 2 dates supplied for specified model
        factor index */
    double gFactor(const DateTime& dateFrom,
                   const DateTime& dateTo,
                   int   factorIdx) const;

    /*  Calculate the 'stepped' tau given a date. */
    const DateTime& sTau(const DateTime& date) const;

	/*
	*  Compute the xi factor as explained in the doc (section 3.1.2)
	*  Xi is the deterministic part of the measure change from the 
	*  tau-forward measure
	*  to the annuity measure. 
	*  From swapvol::XiFactor
	*/
	void xiFactor(
		double&               xi,    /* (O) XiFactor */
		double&               der,   /* (O) derivative of xi with respect to tau */
		double&               maxXi, /* (O) maxXi is the vol of f(t,tau) */
		const vector<double>& weights,
		const DateTimeArray&  CpnPayDates,
		int                   NumCpns,
		const DateTime&       Tau,   /*(I) Tau value in the Newton-Raphson loop */
		const DateTime&       baseDate) const;       /* (I) ZeroBaseDate */
 
	// this function is called on the edfForwardDates dates
    void populatePartialIntIR(const DateTime& today, 
                              const DateTimeArray& myDatesIn, 
                              vector<double> partialIntegralOut[]);

    void populatePartialZeta(const DateTime& today, 
		                     const DateTimeArray& myDates, 
							 vector<double> zeta[]);

private :

    void calcFwdRateAndRBar();
    
	void populateGfactorArrayIR();

    void calcEffRateLimitIR_BM(double          NbSigmasMax, 
                               double          NbSigmasMin,
                               vector<double>& MaxEffRateIR,
                               vector<double>& MinEffRateIR) const;
    
	void calcEffRateLimitIR(
        double    NbSigmasMax,   // (I) Number or sigmas to cut at
        double    NbSigmasMin,   // (I) Number or sigmas to cut at
        vector<double>& MaxRate, /* (O) Cutoff forward rates  */
        vector<double>& MinRate) const;  /* (O) Cutoff forward rates      */
    
    vector<double>  RBar;    // on extendedTimeLine

    DateTimeArray   tau;

    SRMRatesFactorModelSP model;
    double          qLeft;
    double          qRight;
    double          fwdShift;
};

DECLARE(SRMRatesHJMUtil);

DRLIB_END_NAMESPACE
#endif
