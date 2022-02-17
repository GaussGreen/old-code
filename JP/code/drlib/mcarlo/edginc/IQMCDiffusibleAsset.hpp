//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IQMCDiffusibleAsset.hpp
//
//   Description : An interface of hooks for asset diffusion
//
//
//----------------------------------------------------------------------------

#ifndef EDR_IQMCDIFFUSIBLEASSET_HPP
#define EDR_IQMCDIFFUSIBLEASSET_HPP

#include "edginc/IQMCDiffusibleAssetBase.hpp"
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

/************************************************************************
 * These classes form the interface for all the diffusible assets       *
 * they are convenient for a few generic methods acting on all the      *
 * different assets uniformly and for model-independent access to the   *
 * data the models produce                                              *
 ************************************************************************/

/***************************************************************************
 *  The model-independent interfaces to specific diffusible asset types    *
 ***************************************************************************/

/** The model-independent interface to a diffusible Interest Rate asset :
    This class should be a base class for any Rates diffusion model, e.g. SRMRatesHJM or SRMRatesLibor, etc */
class SVQmcDiscFactor;
class SVExpectedDiscFactor;
class IQMCDiffusibleFX;

class MCARLO_DLL IQMCDiffusibleInterestRate : public IQMCDiffusibleAsset
{
public:

    /**   Declaration of pure virtual methods to be overridden by actual models */

    /** (1) methods for engine control */

        /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) = 0;
        /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates) = 0;
    /** getting the simulation start date */
    virtual DateTime getBaseDate() = 0;


    /** (2) Getting the simulated DiscountingFactor between [t_0, measDate] */

    /** Accessing the diffused DF, i.e. an inverse of a value of the
        accumulator account. */
    virtual double  getDiscFactor(SpotIdx measurementDateIdx)           =0;


    /** (3) Getting the simulated ExpectedDF between [measDate, futureDate] */

    /** Accessing the expected value ExpDF(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedDiscFactor(size_t ycIdx,
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx)           =0;

    /** Accessing the natural log of the expected value ExpDF(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedDiscFactor(size_t ycIdx,
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx)           =0;

    /** getting the original price of zero bond maturing at a one of the future dates */
    virtual double getOriginalLnDiscFactor(SpotIdx futureDateIdx)   =0;

    virtual double getOriginalLnExpectedDiscFactor(size_t ycIdx,
                                                    FwdIdx measurementDateIdx,
                                                    FwdIdx futureDateIdx) = 0;


	/** Getting indicator on whether the yc is a diffusion curve(true)
	    or discounting curve (false) */
	virtual size_t registerYCFlavor(IYieldCurveConstSP ) = 0;//{ return size_t(-1); }
    virtual size_t getDiscYCIdx(void)= 0;// { return size_t(-1);} // returns Idx of Discount YC
    virtual size_t getDiffYCIdx(void)= 0;// { return size_t(-1);} // returns Idx of Diffusion YC
    typedef pair<IYieldCurveConstSP, vector<double> > YCForwards; // we pair YC and precalculated Forwards
    typedef vector<YCForwards> YCForwardsDB; // type to keep all seen YCs and storage for the relevant Forwards





//------------------------------------------------------------------------------
    /** (4) creating state variables */

    /** creating a state variable to access the DiscFactor information */
    virtual SVQmcDiscFactor* createDFSV(const DateTimeArray& measurementDates, bool mm);

    /** creating a state variable to access the ExpDiscFactor information */
    virtual SVExpectedDiscFactor* createExpDiscFactorSV(
                                const DateTime&         measurementDate,
                                const DateTimeArray&    futureDates,
                                bool                    doLog,
                                YieldCurveConstSP       yc,
                                bool mm)    ;


//------------------------------------------------------------------------------
    /** this is for efficient access to DiscFactor values, but optional */
    virtual double* getInternalDFArray() { return NULL; }

    /** Allow access to domLnMONEY variable */
    virtual vector<double>& getDomLnMONEY() = 0;

    virtual const vector<double>* getSigmaFX() = 0;
    virtual const vector<double>* getSpotFXFullPath() = 0;

    // TODO: should make pure abstract
    virtual IQMCDiffusibleFX* getFXAsset() {return 0;}

};

DECLARE(IQMCDiffusibleInterestRate);




/** The model-independent interfaces to a diffusible Credit Spread asset*/
class SVSurvivalDiscFactor;
class SVExpSurvDiscFactor;
class SVAggregatedSurvDiscFactor;
class SVDateOfDefault;

class MCARLO_DLL IQMCDiffusibleCreditSpreadBase : public IQMCDiffusibleAsset
{
public:
    /**   Declaration of pure virtual methods to be overridden by actual models */

    /** (1) methods for engine control */

        /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) = 0;
        /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates) = 0;
    /** getting the simulation start date */
    virtual DateTime getBaseDate() = 0;


    /** (2) Getting the simulated DiscountingFactor between [t_0, measDate] */

    /** Accessing the diffused SDF, i.e. a non-default probability by
        a given date */
    virtual double  getSurvivalDiscFactor(SpotIdx measurementDateIdx)   =0;
    /** getting the simulated recovery rates (same idx as for SDF) */
    virtual double getRecoveryRate(SpotIdx idx) = 0;

    virtual double getOriginalLnSurvivalDiscFactor(SpotIdx measurementDateIdx) = 0;
    virtual double getOriginalLnExpectedSurvivalDiscFactor(FwdIdx measurementDateIdx, FwdIdx futureDateIdx) = 0;


    /** (3) Getting the simulated ExpectedSDF between [measDate, futureDate] */

    /** Accessing the expected value ExpSDF(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx)           =0;

    /** Accessing the natural log of the expected value ExpSDF(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx)           =0;

    /**  Getting expected recovery rate for date futureDateIdx conditional 
         on the measurementDateIdx filtration */
    virtual double getExpRecoveryRate(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx)           =0;
//------------------------------------------------------------------------------
    /** (4) creating state variables */

    /** creating a state variable to access the SurvivalDiscFactor information */
    virtual SVSurvivalDiscFactor* createSDFSV(
                                const DateTimeArray& measurementDates, bool mm);

    /** creating a state variable to access the ExpDiscFactor information */
    virtual SVExpSurvDiscFactor* createExpSurvDiscFactorSV(
        const DateTime&         measurementDate,
        const DateTimeArray&    futureDates,
        bool                    doLog, bool mm);

    /** creating a state variable to gain compact access to a set of SDFs and ESDFs */
    SVAggregatedSurvDiscFactor*   createAggregatedSurvDiscFactorSV(
            DateTimeArraySP   sdfDates,            // set of dates for discount factors
            SpotIdxArraySP    sdfIdxSP,
            DateTimeArraySP   esdfRequestedDates,  //  union of {t_i} for all expected factors
            FwdIdxArraySP     esdfReqIdxSP,
            DateTimeArraySP   esdfForwardDates,    // union of all {T_j} for all expected factors
            FwdIdxArraySP     esdfForIdxSP,
            const DateTime&   maxDiffDate,
            const DateTime&   maxCurveMat,
            bool              doLog,
            bool              mm);

//------------------------------------------------------------------------------
    /** this is for efficient access to SurvivalDiscFactor values, but optional */
    virtual double* getInternalSDFArray() { return NULL; }

    // these two are necessary for handling date of default for simple and composite spreads. 
    // Default implementation - to do noting.
    virtual void setWholeTimelineSurvProbRequest() { }
    virtual double getWholeTimelineLogSurvProb(size_t idx) { return 0.0; } 

    // has default implementation of returning nothing
    virtual IQMCDiffusibleInterestRateSP getUnderlyingIRAsset() {return IQMCDiffusibleInterestRateSP();}
};

DECLARE(IQMCDiffusibleCreditSpreadBase);

/** additional optional pure interface that some credit models might choose to implement -- 
    it is used only under special composition rules -- when composition weight is fractional.
    On the other hand, if the model did not implement this interface - it cannot be used with 
    a fractional coefficient in the CreditComposite. */

class MCARLO_DLL IQMCCreditSpreadAccessingFractional : public virtual IQMCDiffusibleCreditSpreadBase
{
public:
    /** Accessing the expected value E[ SDF(md, fd)^fraction] where md is a
        simulated measurement date and fd is some future date after the
        measurement is made, and a fraction is a power. */
    virtual double  getExpectedFractionalSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx,
                                        double fraction)           =0;

    /** Accessing the natural log of the expected value E[ SDF(md, fd)^fraction]
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedFractionalSurvivalDiscFactor(
                                        FwdIdx measurementDateIdx,
                                        FwdIdx futureDateIdx,
                                        double fraction)           =0;
};

DECLARE(IQMCCreditSpreadAccessingFractional);

/** this class is here to enforce the virtual inheritance from IQMCDiffusibleDefaultableCreditSpreadBase
    which is crucial for correct implementation of IQMCCreditSpreadAccessingFractional, i.e.
    it completes the "diamond" of multiple inherited class */

class MCARLO_DLL IQMCDiffusibleCreditSpread : public virtual IQMCDiffusibleCreditSpreadBase {};
DECLARE(IQMCDiffusibleCreditSpread);




class MCARLO_DLL IQMCDiffusibleDefaultableCreditSpread : public IQMCDiffusibleCreditSpread
{
public:


    /** this is an extension of the Credit Spread Interface to include those
        objects that support DateOfDefault and RecoveryRate enquiries */

    /** Informs asset that the date of default will be asked for a path */
    /** using this method shall make it illegal to call SDF for this asset */
    virtual void setDateOfDefaultEnquiry() = 0;

    /** Retrieves the simulated date of default */
    virtual DateTime getDateOfDefault() = 0;

    /** Retrieves the recovery rate as of date of default */
    virtual double getRecoveryRateAtDefault() = 0;

    /** creating a state variable to access the DateOfDefault information */
    virtual SVDateOfDefault* createSVDateOfDefault();
};

DECLARE(IQMCDiffusibleDefaultableCreditSpread);








/** The model-independent interfaces to a diffusible FX Rate asset*/
class SVSpotFX;
class SVExpectedFX;

class MCARLO_DLL IQMCDiffusibleFX : public IQMCDiffusibleAsset
{
public:
    /**   Declaration of pure virtual methods to be overridden by actual models */

    /** (1) methods for engine control */

        /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) = 0;
        /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates) = 0;
    /** getting the simulation start date */
    virtual DateTime getBaseDate() = 0;

    /** (2) Getting the simulated DiscountingFactor between [t_0, measDate] */

    /** Accessing the diffused SpotFX, i.e. a non-default probability by
        a given date */
    virtual double  getSpotFX(SpotIdx measurementDateIdx)               =0;


    /** (3) Getting the simulated ExpectedFX between [measDate, futureDate] */

    /** Accessing the ExpectedFX(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedFX(
                                FwdIdx measurementDateIdx,
                                FwdIdx futureDateIdx)                   =0;

    /** Accessing the natural log of the ExpectedFX(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedFX(
                                FwdIdx measurementDateIdx,
                                FwdIdx futureDateIdx)                   =0;

//------------------------------------------------------------------------------
    /** (4) creating state variables */

    /** creating a state variable to access the SurvivalDiscFactor information */
    virtual SVSpotFX* createSpotFXSV(const DateTimeArray& measurementDates, bool mm);

    /** creating a state variable to access the ExpDiscFactor information */
    virtual SVExpectedFX* createExpectedFXSV(
        const DateTime&         measurementDate,
        const DateTimeArray&    futureDates,
        bool                    doLog, 
        bool                    mm);


//------------------------------------------------------------------------------
    /** this is for efficient access to SpotFX values, but optional */
    virtual double* getInternalSpotFXArray() { return NULL; }

    /** Few auxiliary methods that might be implemented in the model
        specific way if they are necessary for efficient use, otherwise
        a default implementation will be used instead */

    virtual IQMCDiffusibleInterestRateSP getDomesticIRAsset() = 0 ;
    virtual IQMCDiffusibleInterestRateSP getForeignIRAsset() = 0 ;

    // for moment matching
    virtual double getOriginalLnFwdSpot(SpotIdx dateIdx) = 0; // where date is from getSpotDateIdx list
    virtual double getOriginalLnExpectedFwdSpot(FwdIdx dateIdx) = 0; // where date is from getForwardForwardDateIdx list
 
};

DECLARE(IQMCDiffusibleFX);

/** The model-independent interfaces to a diffusible FX Rate asset*/
class SVSpotEQ;
class SVExpectedEQ;

class MCARLO_DLL IQMCDiffusibleEQ : public IQMCDiffusibleAsset
{
public:
    /**   Declaration of pure virtual methods to be overridden by actual models */

    /** (1) methods for engine control */

        /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) = 0;
        /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates) = 0;
    /** getting the simulation start date */
    virtual DateTime getBaseDate() = 0;
    
    /** (2) Getting the simulated DiscountingFactor between [t_0, measDate] */

    /** Accessing the diffused SpotFX, i.e. a non-default probability by
        a given date */
    virtual double  getSpotPrice(SpotIdx measurementDateIdx)            =0;

    /** (3) Getting the simulated ExpectedFX between [measDate, futureDate] */

    /** Accessing the ExpectedFX(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double  getExpectedPrice(
                                FwdIdx measurementDateIdx,
                                FwdIdx futureDateIdx)                   =0;

    /** Accessing the natural log of the ExpectedFX(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedPrice(
                                FwdIdx measurementDateIdx,
                                FwdIdx futureDateIdx)                   =0;

//------------------------------------------------------------------------------
    /** (4) creating state variables */

    /** creating a state variable to access the SurvivalDiscFactor information */
    virtual SVSpotEQ* createSpotEQSV(const DateTimeArray& measurementDates, bool mm);

    /** creating a state variable to access the ExpDiscFactor information */
    virtual SVExpectedEQ* createExpectedEQSV(
        const DateTime&         measurementDate,
        const DateTimeArray&    futureDates,
        bool                    doLog, 
        bool                    mm);


//------------------------------------------------------------------------------
    /** this is for efficient access to SpotFX values, but optional */
    virtual double* getInternalSpotPriceArray() { return NULL; }

    virtual IQMCDiffusibleInterestRateSP getDomesticIRAsset() = 0 ;
   

    // for moment matching
    virtual double getOriginalLnFwdSpot(SpotIdx dateIdx) = 0;  // where date is from getSpotDateIdx list
    virtual double getOriginalLnExpectedFwdSpot(FwdIdx dateIdx) = 0; // where date is from getForwardForwardDateIdx list
};

DECLARE(IQMCDiffusibleEQ);


///// For energy assets ////////////////////////////////////////////////
class SVExpEnergyFuture;

class MCARLO_DLL IQMCDiffusibleEnergy : public IQMCDiffusibleAsset
{
public:

    /**   Declaration of pure virtual methods to be overridden by actual models */

    /** (1) methods for engine control */

        /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) = 0;
        /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates) = 0;
	/** getting the simulation start date */
	virtual DateTime getBaseDate() = 0;

    /** (2) Getting the simulated DiscountingFactor between [t_0, measDate] */
    // Accessing the expected value ExpFP(md, fp) where md is a simulated measurement date and fd is
	// some future date after the measurement is made.
    virtual double  getExpectedPrice(	FwdIdx measurementDateIdx,
										FwdIdx futureDateIdx)           =0;

	/** Accessing the natural log of the ExpectedFP(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
    virtual double	getLnExpectedPrice( FwdIdx measurementDateIdx,
										 FwdIdx futureDateIdx)                   =0;

//------------------------------------------------------------------------------
    /** (4) creating state variables */

    // creating a state variable to access the expected energy future information
    virtual SVExpEnergyFuture* createExpEnergyFutureSV(
        const DateTime&         measurementDate,
        const DateTimeArray&    futureDates,
        bool                    doLog,
        bool                    mm);

//------------------------------------------------------------------------------
    /** this is for efficient access to DiscFactor values, but optional */
    virtual double* getInternalFuturePriceArray() { return NULL; }


};

DECLARE(IQMCDiffusibleEnergy);

//////////////////////////////////////////////////////////////////////////

/** The model-independent interface to a diffusible basis coupon asset :
    This class should be a base class for any basis index diffusion model,
    e.g. SRMSP or BGMSP, etc */
class SVExpectedBasisFwdSpread;

class MCARLO_DLL IQMCDiffusibleBasisIndex : public IQMCDiffusibleAsset
{
public:

    /**   Declaration of pure virtual methods to be overridden by actual models */

    /** (1) methods for engine control */

    /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/) = 0;

    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates) = 0;

    /** getting the simulation start date */
    virtual DateTime getBaseDate() = 0;

    /** (2) Getting the simulated expected Fws spread between
            [measDate, futureDate] */

    /** Accessing the expected basis coupon between md and fd where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    // TO DO:  FwdIdx should be interface specific so IR ones can't be passed
    //         to CR or Basis ones.
    virtual double getExpectedBasisFwdSpread(
        FwdIdx measurementDateIdx,  // Idx of measDate
        FwdIdx resetDatesIdx ) = 0;  // Idx of forward requested

//------------------------------------------------------------------------------
    /** (4) creating state variables */

    /** creating a state variable to access the ExpectedBasisFwdSpread information */
    // TO DO:  Why not wrap in a smart pointer?
    // RESOLUTION:  Drop circular dependence of headers and then change to SP
    virtual SVExpectedBasisFwdSpread* createExpBasisFwdSpreadSV(
        const DateTime& measurementDate,
        const DateTimeArray& resetDates );


    virtual IQMCDiffusibleInterestRateSP getUnderlyingIRAsset() = 0 ;

};

DECLARE(IQMCDiffusibleBasisIndex);

class MCARLO_DLL DiffusionDriver
{
public:
    DiffusionDriver(long val=0) : driver(val) {}
    ~DiffusionDriver() {}
    void setDriver(long val) { driver = val;  }
    long setDriver() const   { return driver; }

protected:
    long    driver;
};


DRLIB_END_NAMESPACE
#endif //EDR_IQMCDIFFUSIBLEASSET_HPP

