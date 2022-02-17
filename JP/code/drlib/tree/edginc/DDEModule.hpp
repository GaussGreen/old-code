//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : DDEModule.hpp
//
//   Description : utility class for DDE calculations: calibration etc
//
//   Author      : Qing Hou
//
//----------------------------------------------------------------------------

#ifndef EDG_DDE_MODULE_H
#define EDG_DDE_MODULE_H

#include "edginc/Class.hpp"
#include "edginc/CleanSpreadVolCurve.hpp"
#include "edginc/EquitySpreadCorrCurve.hpp"
#include "edginc/DDEParams.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Theta.hpp"
#include "edginc/HaveNonEqSpread.hpp"
#include "edginc/MaturityPeriod.hpp"

DRLIB_BEGIN_NAMESPACE

class CAsset;
class AssetDDE;
class IVolProcessed;
class CVolBase;
class CVolRequest;
class TimeMetric;
class IDDEInitiator;
class SpreadEquityFunc;
class Control;
class Results;
class DDEOutSens;
class DDECalib;
class IDDEInitiator;

/** collection of utility function used to calibrate generic (2 factor) DDE equity 
	bare vol and related calculations */
class TREE_DLL DDEModule: public CObject, virtual public Theta::IShift, virtual public IHaveNonEqSpread
{
public:
// not done yet
    static CClassConstSP const TYPE;
    friend class DDEModuleHelper;
	friend class DDEModuleFD1FPricer;
	friend class DDEModuleFD1FPricerBase;
	friend class FD1FDDECalib;

	static const double DEF_IMPVAR_EPSILON;

	//typedef enum _DDEType { RG = 0, LN1F, LN2F } DDEType;

	static smartPtr<DDEModule> createModule(const string &ddeType, const DDECalib *calibParams, const AssetDDE *asset, const IDDEInitiator *initiator);

	static CleanSpreadVolCurveSP createCleanSpreadVolCurve(
		const CleanSpreadCurve *cleanSprds, const DDEParams *ddeParams, const TimeMetric *timeMetric);

	static EquitySpreadCorrCurveSP createEquitySpreadCorrCurve(
		const CleanSpreadCurve *cleanSprds, const DDEParams *ddeParams);

	static CleanSpreadCurveSP calibLNCleanSpreadCurve(
		const CleanSpreadCurve *cleanSprds, const CleanSpreadVolCurve *vols, const DDECalib *calib); 

	// return the clean spread curve for equity drift, due to participation
	// may be different from the name's original clean spread curve
	// the non-participation portion is returned by getNonEqSpreads
	const CleanSpreadCurve *getCleanSpreadCurve() const;

	const CleanSpreadCurve *getNonEqSpreads() const;

	IVolProcessed *getProcessedVol(const CVolRequest* volRequest) const;
	
	void recordOutputRequests(Control* control,
							Results* results) const;

	bool sensShift(Theta *shift);

	void update(const AssetDDE *asset, const IDDEInitiator *initiator);

	// get a spread as function of equity
	const SpreadEquityFunc *getSpreadFunc() const;
	
	void getSpreadCurves(const CleanSpreadCurve *&bbSprds,
						const CleanSpreadVolCurve *&vols,
						const EquitySpreadCorrCurve *&corrs) const;

    // routines for visualization of DDE model
    smartPtr<DDEOutSens> calcDDEVisual(const OutputRequest *request) const;

    smartPtr<CVolBase> calcDDEImpVol(const OutputRequest *request) const;

protected:
    DDEModule();
	DDEModule(CClassConstSP clazz, bool isOneFactor);
	
	void createCleanSpreadCurve();

	void createNonEqCleanSpreadCurve();

	// build template to calibrate to vanilla market and do integration
	void buildTemplate();

    // interp the template to get the strike
    double getCalibStrike(const DateTime &date) const;

	void smartCalib();

	// compute riskyness adj spread, as well as sprd vol and eq-sprd corr
	void preprocessCredit();

	virtual void calib();

	// calibrate risky zero and vanilla prices. return the risky equity var. 
	// create spread function if is 1 factor, or spread backbone if 2 factor
	virtual void calib1F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps, const DoubleArray &sprdVars,
			const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars);
	
	virtual void calib2F(const DateTimeArray &dates, const DoubleArray &t, const DoubleArray &pvs, const DoubleArray &ndps,
			const DoubleArray &fwds, const DoubleArray &strikes, const DoubleArray &bsVars, const DoubleArray &calls, DoubleArray &riskyVars);

	void store();

	void restore();

	void restoreCreditOnly();

	bool riskyVolOW() const; // if has risky vol overwrite

protected:
	bool		isOneFactor;
    DateTime	valueDate;
	DateTime	lastDate;
	const AssetDDE *asset; // $unregistered

	// internal transient members
	ScheduleSP					calibSchd;			// contain time/strkTemplate

	CleanSpreadCurveSP			cleanSprds;			// eq participation portion
	CleanSpreadCurveSP			nonEqCleanSprds;	// eq non participation portion
	CleanSpreadCurveSP			bbCleanSprds;
	CleanSpreadVolCurveSP		cleanSprdVols;
	EquitySpreadCorrCurveSP		corrs;
	smartPtr<CVolBase>			riskyEqVols;
	smartPtr<SpreadEquityFunc>	spreadFunc;

	DateTime					lastGoodCalibDate; // $unregistered

	// temp mem for sensitivity
	CleanSpreadCurveSP			cleanSprdsBak; // $unregistered
	CleanSpreadCurveSP			nonEqCleanSprdsBak; // $unregistered
	CleanSpreadCurveSP			bbCleanSprdsBak; // $unregistered
	CleanSpreadVolCurveSP		cleanSprdVolsBak; // $unregistered
	EquitySpreadCorrCurveSP		corrsBak; // $unregistered
	smartPtr<CVolBase>			riskyEqVolsBak; // $unregistered
	smartPtr<SpreadEquityFunc>	spreadFuncBak; // $unregistered

	// calib date/strike required by product
    const IDDEInitiator			*initiator; // $unregistered
    bool                        calibStrkIsPct;

    smartConstPtr<DDECalib>      calibParams; // calib parameters normally inside a model
};

typedef smartConstPtr<DDEModule> DDEModuleConstSP;
typedef smartPtr<DDEModule>      DDEModuleSP;

//////////////////////////////////////////////
//											//
//			IDDEInitiator					//
//											//
//////////////////////////////////////////////

/** classes that hold information to customize the dde module */
class TREE_DLL IDDEInitiator {
public:
	virtual DateTime maxMaturity() const=0;

	virtual void sensitiveDates(DateTimeArray &dates) const;

	// given dates, get strikes. dates may include model required dates plus the sensitiveDates.
	virtual void sensitiveStrikes(	const DateTimeArray		dates,
									DoubleArray				&strikes,	// same dimension as dates
									bool					&strikeIsPct) const;

	static void load(CClassSP& clazz);

	virtual ~IDDEInitiator();
};


//////////////////////////////////////////////
//                                          //
//          DDECalib                        //
//                                          //
//////////////////////////////////////////////

/** class to hold calibration info/method */

/** base class */
class TREE_DLL DDECalib: public CObject
{
public:
    static CClassConstSP const TYPE;
    friend class DDECalibHelper;

    static const string DEF_TEMPLATE_FREQ;
    static const int DEF_NB_STEPPERYEAR;
    static const int DEF_NB_STEPPERBENCHMK;
    static const double DEF_MAX_STEP_SIZE;
    static const double DEF_TOLERANCE_SPRD;
    static const double DEF_TOLERANCE_VOL;
    static const string DEF_MAX_HORIZON;
    static const string DEF_MIN_HORIZON;

    static void load(CClassSP& clazz);
    
    DDECalib();

    virtual void validatePop2Object();

    // create date template
    DateTimeArraySP getTemplateDates(DateTime valueDate,
                                     DateTime lastDate) const;

    static ExpiryArraySP createDefExpiries();

    DDECalib(const CClassConstSP &clazz);

public:
	// parameters to generate calibration benchmark
    bool                useProdDateOnly;// 1 if don't insert dates using freq/expiries
    bool                useEqualStep;   // 1 if use template freq, 0 if use expiries
    MaturityPeriodSP    templateFreq;   // benchmark for calib of vol (eq/sprd) etc
    ExpiryArraySP       expiries;       // benchmark for calib of vol (eq/sprd) etc
    
    // parameters to add finer steps between benchmarks for pricing
    int                 nStepPerBenchMk; // finer step between benchmarks for calibration
    double              maxStepSize;    // max size of finer step, in yrs
    int                 nStepPerYear;   // finer step between benchmarks for calibration
    
	double              tolSprd;        // spread tolerance
    double              tolEqVol;       // equity vol tolerance
    bool                useAtm;         // whether always calibrate to atm vol
    bool                atmIsFwd;       // 1 if ATM means atm fwd, 0 if atm spot

    // when calib strk too low & off surface, whether use the lowest strk on surface
    bool                strkOffSurfaceAdj; 
	
	bool				useRiskyVolOW;	// if use risky vol overwrite when provided

	// throw error if fail to calibrate, esp useful if calibration is successful for pricing, 
	// but fail to calibrate to maturity for tweaks, which cause sensitivity number to be wrong
	bool				throwIfFail;	
    // max horizon for calibration, in case deal has very long maturity
    // similarly have min horizon
    ExpirySP            maxHorizon;
    ExpirySP            minHorizon;
};

typedef smartConstPtr<DDECalib> DDECalibConstSP;
typedef smartPtr<DDECalib>      DDECalibSP;

class TREE_DLL DDECalibMC: public DDECalib // Monte Carlo calibration
{
public:
// not done yet
    static CClassConstSP const TYPE;
    static const int DEF_NB_PATH;
    friend class DDECalibMCHelper;

    static void load(CClassSP& clazz);

    DDECalibMC();

public:
    int             nbPath;
    IRandomSP       rand;               // random number generator for calibration
};

typedef smartConstPtr<DDECalibMC> DDECalibMCConstSP;
typedef smartPtr<DDECalibMC>      DDECalibMCSP;


class TREE_DLL DDECalibFD: public DDECalib // FiniteDiff calibration
{
public:
// not done yet
    static CClassConstSP const TYPE;
    static const int DEF_NB_STOCKSTEP;
    static const int DEF_NB_MAX_ITER;
    static const int DEF_NB_JACOBI;
    static const double DEF_TRUNCATION;
    static const double DEF_CREDIT_STRK;
    static const int DEF_NB_SEGMENT_ENDS;
    friend class DDECalibFDHelper;

    virtual void validatePop2Object();

    static void load(CClassSP& clazz);

	static MaturityPeriodArraySP createDefSegExpiries();
    static IntArraySP createDefSegStepPerYears();

    DDECalibFD(bool isForward=true);

public:
    int         nStockSteps;
    int         gridType;
    int         maxIter;
    double      TruncationStd;
    int         maxNbJacobi;    // max nb of jacobi calculation in root finding
    MaturityPeriodArraySP segExpiries;  // expiries for FD segments, limit to MaturityPeriod since we don't want to worry about historical dates
    IntArraySP  segStepPerYears;    // segment specific steps per year    

    // parameters used for forward calibration
    bool        isForward;
    bool        firstStepSmooth;// true if do smoothing for 1st time step using Black price
    bool        useDensity;     // true if forward calc on function C - K * dC/dK
    double      creditStrike;   // strike relative to ATM. used to calib spread from 
                                // vanilla call's sensitivity to strike
};

typedef smartConstPtr<DDECalibFD> DDECalibFDConstSP;
typedef smartPtr<DDECalibFD>      DDECalibFDSP;


DRLIB_END_NAMESPACE
#endif
