//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditCIRDiffuse.hpp
//
//   Description : CIR credit path generation
//
//
//----------------------------------------------------------------------------

#ifndef SRMCREDITCIRDIFFUSE_HPP
#define SRMCREDITCIRDIFFUSE_HPP

#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/DECLARE.hpp"
#include <set>

#include "edginc/SRMCreditCIRUtil.hpp"
#include "edginc/SRMCreditUtil.hpp"
#include "edginc/SRMCreditDiffuse.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"

DRLIB_BEGIN_NAMESPACE

/** class for CIR credit path generation */
class SRMCreditCIRDiffuse: public SRMCreditDiffuse, public IQMCCreditSpreadAccessingFractional 
 {
public:

    /** constructor */
    SRMCreditCIRDiffuse(IQMCDiffusibleInterestRateSP srmRatesDiffuse, bool isFullMC);

    /** destructor */
    virtual ~SRMCreditCIRDiffuse();

    /** initialization */
    void setSRMCreditCIRDiffuse(
        int                    _randomIndex,
        const DateTime&        _today,
        SRMCreditCIRUtilSP     _srmCreditCIRUtil,
        double                 _crFxCorr,
        const vector<double>&  _prob);

    /** finalizes the timeline, allocates memory */
    virtual void finalizePathGenerator(DateTimeArrayConstSP allDates);

    /** generates path across all dates */
    virtual void generatePathAndSurvivalRates(
        IQMCRNGManagerSP rngMgr);

    /** for capturing state of diffusion */
    struct DiffusedState {
        double intensity;
    };

    /** for capturing diffused state on each date on which calculations are performed */
    typedef vector<DiffusedState>::const_iterator DiffusedStateIter;

    /** accesses the expected value ExpSDF(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double getExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j);

    /** accesses the natural log of the expected value ExpSDF(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j);

    /** Accessing the expected value E[ SDF(md, fd)^fraction] where md is a
        simulated measurement date and fd is some future date after the
        measurement is made, and a fraction is a power. */
    virtual double  getExpectedFractionalSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx,
                                        double fraction);

    /** Accessing the natural log of the expected value E[ SDF(md, fd)^fraction]
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedFractionalSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx,
                                        double fraction);


protected:

	/** fix all the model parameters except theta relatd (includes initial values */
	void setModelParametersExceptTheta(double _c, double _d, double _beta);

	/** initialize what is related to theta dates */
	void setThetaDates(const DateTimeArray& dates);
	/** from a series of dates and surv proba, find theta and initial value */
	void calibrateTheta(const DateTimeArray& dates, const vector<double> &logSurvProb);
	/** force theta and initial value to be a certain value */
	void setTheta(const DateTimeArray& dates, double _initialIntensity, const vector<double> &_theta);




    /** populates simulationTheta by looking up in piecewise constant theta */
    void storeSimulationTheta();

    /** populates esdfForwardYearFracs and expHesdfForwardYearFracs */
    void storeMeasurementTAndExpHT(const DateTimeArray& dates);

    /** returns B(x,y) - see model documentation */
    double b2(double x, double y, double expHx, double expHy);

    /** core function used in a2 - see model documentation */
    double g2(double x, double y, double expX, double expY);

    /** core function used in a2 - see model documentation */
    double f2(double x, double y, double expX, double expY);

    /** returns A(x,y) - see model documentation */
    double a2(double x, double y, double expX, double expY);

	/** returns B(x,y) for the process frac*CIR */
	double b2frac(double x, double y, double frac, double hfrac);

	/** core function used in a2 for the process frac*CIR */
	double g2frac(double x, double y, double frac, double hfrac);

	/** core function used in a2 for the process frac*CIR */
	double f2frac(double x, double y, double frac, double hfrac);

	/** returns A(x,y) for the process frac*CIR */
	double a2frac(double x, double y, double frac, double hfrac);

    /** a function to compute expected survival probabilities as of t=0 */
    void computeLogDetermSurvivalProb(const DateTimeArray& dates, vector<double> &logSurvProbs);


private:

    vector<DiffusedState>   expProb;            // stores state variables

    double          beta;                       // mean reversion speed - see model documentation
    double          c;                          // vol parameter - see model documentation
    double          d;                          // vol parameter - see model documentation
    double          h;                          // sqrt(beta^2 + 2 c)
    double          initialIntensity;           // initial value for diffusion
//    double          normalIntensity;            // used in vol parameters - see model documentation
    bool			cIsZero;                    // TRUE if c = 0 (with tolerance)
    bool			hIsZero;                    // TRUE if h = 0 (with tolerance)
    double          minIntensity;               // smallest allowable intensity (-d/c)

    bool            fullMC;                     // false means "no diffusion, take expected values instead"

 //   double qLeft;
 //   double qRight;
 //   double vol;     // sigma

    vector<double>  theta;                      // mean reversion levels - see model documentation
    int             numThetas;                  // size of theta
    vector<double>  yearFracsForTheta;          // theta is piecewise constant b/w these times
    vector<double>  expHYearFracsForTheta;      // exp(h*yearFracsForTheta) - stored for speed

    vector<double>  simulationTheta /* FIXME: calc on fly*/;            // values of theta stored at simulation dates
    vector<double>  sqrtYearFrac;  /*FIXME: move to a common class*/             // year fracs between simulation dates

    vector<double>  esdfForwardYearFracs;       // year fracs from today to measurement dates
    vector<double>  expHesdfForwardYearFracs;   // exp(h*esdfForwardYearFracs) - stored for speed

    static const int matYearsForNormalIntensity = 5;  // hard-coded maturity for normalIntensity

    SRMCreditCIRUtilSP  srmCreditCIRUtil;       // util object
    void            trimToDiffusion();          // aux function to save memory
};

/** declare smartPtr versions */
DECLARE(SRMCreditCIRDiffuse);

DRLIB_END_NAMESPACE

#endif // SRMCREDITCIRDIFFUSE_HPP
