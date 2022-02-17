//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ImpliedSample.hpp
//
//   Description : Class that can provide a sample from the implied distribution
//
//   Date        : September 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/Class.hpp"
#include "edginc/Object.hpp"
#include "edginc/Format.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/LinearInterpolator.hpp"

#include "edginc/Asset.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/VolRequestLNStrike.hpp"
#include "edginc/PDFRequestLNStrike.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/Lattice.hpp"
#include "edginc/Function.hpp"


#ifndef EDR_IMPLIEDSAMPLE_HPP
#define EDR_IMPLIEDSAMPLE_HPP

DRLIB_BEGIN_NAMESPACE

/** ImpliedSample: Interface for sampling from an 
    Implied Distribution */
class MARKET_DLL ImpliedSample: virtual public IObject {
public:
    static CClassConstSP const TYPE;

    /** Produces a sample from a uniform random variable */
    virtual double sample(double uniform) const = 0;

    /** Returns whether this is an implied sample of a 
        return or of an absolute price */
    virtual bool isRelativeReturn() const = 0;

private:    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

typedef smartPtr<ImpliedSample> ImpliedSampleSP;
typedef smartConstPtr<ImpliedSample> ImpliedSampleConstSP;

typedef array<ImpliedSampleSP, ImpliedSample> ImpliedSampleArray;
typedef smartPtr<ImpliedSampleArray> ImpliedSampleArraySP;
typedef smartConstPtr<ImpliedSampleArray> ImpliedSampleArrayConstSP;


/** Wraps a LinearInterpolant on the Inverse of the Implied 
    Distribution with optimal lookup values for the Hunt algorithm. 
    Provides adjusted samples to match the forward */
class MARKET_DLL LinearImpliedSample: public CObject,
                           virtual public ImpliedSample {
public:
    static CClassConstSP const TYPE;
    friend class LinearImpliedSampleHelper;

    /** Full constructor */
    LinearImpliedSample(const CDoubleArray& strikes,             // tabulated strikes
                        const CDoubleArray& probs,               // tabulated probs
                        double              fwdFwd,              // the theoretical mean of the distribution
                        int                 nbDown,              // number of strikes used in the left tail
                        int                 nbMid,               // number of strikes used in the middle
                        int                 nbUp,                // number of strikes used in the right tail
                        bool                computeAdjustment,   // adjust the tabulated probs to match the Fwd perfectly
                        bool                isRelative = false); // relative return or absolute price

    /** Produces a sample from the implied distribution by 
        inversion of the cumulative probability by delegating
        it to sampleNonVirtual. Part of the ImpliedSample
        interface. */
    virtual double sample(double uniform) const;

    /** Returns whether this is an implied sample of a 
        return or of an absolute price */
    virtual bool isRelativeReturn() const;

    /** Produces a sample from the implied distribution by 
        inversion of the cumulative probability. Not a virtual call */
    double sampleNonVirtual(double uniform) const;

    /** Computes Probability(S > K) for a given strike K */
    double probability(double strike) const;
    
    /** Computes the density function at a given strike K */
    double density(double strike) const;
    
    /** Computes the expected value of a function under the 
        tabulated probability distribution */
    double expectedValue(const Function1DDouble& func,
                         const Function1DDouble& indefiniteIntegral) const;

    /** Gets the adjustment used to match perfectly the Fwd */
    double getAdjustment() const;

    /** Gets the probability mass in the left tail */
    double leftTailMass();
    
    /** Gets the probability mass in the right tail */
    double rightTailMass();

private:
    LinearImpliedSample();

    double              fwdFwd;                 //!< Asset's forward
    int                 nbDown;                 //!< Number of strikes on left tail
    int                 nbMid;                  //!< Number of strikes in middle
    int                 nbUp;                   //!< Number of strikes on right tail
    bool                computeAdjustment;      //!< Whether to adjust samples to match the fwd

    // Fast non virtual linear interpolant
    LinearInterpolantNonVirtualSP interp;       //!< Interpolant on (distribution, strike)
    
    double              adjustment;             //!< Multiplicative adjustment
    double              leftTailProb;           //!< Probability mass on left tail
    double              rightTailProb;          //!< Probability mass on right tail
    bool                isRelative;             //!< Whether it is a return or absolute price
    
    // Helper value for fast sample method: (nbMid - 1) / (leftTailProb - rightTailProb)
    double              invMidProbPerStrike; 

    // Required for computing the forward
    class MARKET_DLL UnitMap: public Function1DDouble {
    public:
        virtual double operator()(double  x) const {
            return x;   
        };
    };

    // Required for computing the forward
    class MARKET_DLL IntegralUnitMap: public Function1DDouble {
    public:
        virtual double operator()(double  x) const {
            return x * x / 2.0;
        };
    };
};

typedef smartPtr<LinearImpliedSample> LinearImpliedSampleSP;
typedef smartConstPtr<LinearImpliedSample> LinearImpliedSampleConstSP;

typedef array<LinearImpliedSampleSP, LinearImpliedSample> LinearImpliedSampleArray;
typedef smartPtr<LinearImpliedSampleArray> LinearImpliedSampleArraySP;
typedef smartConstPtr<LinearImpliedSampleArray> LinearImpliedSampleArrayConstSP;

//////////////////////////////////////////////////////////////
// DegenerateImpliedSample: a sampler that returns a fixed
// number. Used for situations where the user requires an
// ImpliedSample* point of view but the asset is deterministic.
// This is our view of a Dirac density function
//
// Useful for zero trading time, or for dropping the first AvgIn
// date in MCImplied
//////////////////////////////////////////////////////////////
class MARKET_DLL DegenerateImpliedSample: virtual public ImpliedSample,
                                       public CObject {
public:
    static CClassConstSP const TYPE;
    friend class DegenerateImpliedSampleHelper;

    DegenerateImpliedSample(double fixedSample, bool isRelative);

    /** Produces a sample from the implied distribution by 
        inversion of the cumulative probability. Here the fixed
        sample is returned */
    virtual double sample(double uniform) const;

    /** Returns whether this is an implied sample of a 
        return or of an absolute price */
    virtual bool isRelativeReturn() const;

private:
    DegenerateImpliedSample();    
    
    double fixedSample;     //!< Deterministic sample
    bool   isRelative;      //!< Whether it is a return or absolute price
};

typedef smartPtr<DegenerateImpliedSample> DegenerateImpliedSampleSP;
typedef smartConstPtr<DegenerateImpliedSample> DegenerateImpliedSampleConstSP;


///////////////////////////////////////////
// PDFPARAMS
///////////////////////////////////////////
// class that contains parameters required in computing a 2-stage strikes/distributions pair
class MARKET_DLL PDFParams: public CObject {
public:
    static CClassConstSP const TYPE;
    
    PDFParams();
    
    PDFParams(int    numMidStrikesPerLogStrike,
              int    numTailStrikesPerLogStrike,
              int    numMiddleStrikes,           
              int    numTailStrikes,
              double propMiddleStrikes,
              double numMiddleStdDevs,
              double numTailStdDevs,
              int    maxStdDevs,
              int    maxFineStdDevs,
              bool   failProbs);

    int    numMidStrikesPerLogStrike;   //!< Number of strikes per logstrike used in the middle of the distribution
    int    numTailStrikesPerLogStrike;  //!< Number of strikes per logstrike used in the tail of the distribution
    int    numMiddleStrikes;            //!< The minimum number of strikes that can be used in the middle
    int    numTailStrikes;              //!< The minimum number of strikes that can be used in the tail
    double propMiddleStrikes;           //!< The proportion of middle strikes that will be used in the first pass
    
    double numMiddleStdDevs;            //!< The number of standard deviations that describe the middle of the distribution
    double numTailStdDevs;              //!< The number of standard deviations that describe the tail of the distribution

    int    maxStdDevs;                  //!< Maximum number of standard deviations used for bracketing the desired stdDevs
    int    maxFineStdDevs;              //!< Maximum number of standard deviations used for finer bracketing of the desired stdDevs

    bool   failProbs;                   //!< Fail Vs Massaging ill-defined densities

    void validatePop2Object();

    // DEFAULT VALUES
    static int    default_numMidStrikesPerLogStrike;
    static int    default_numTailStrikesPerLogStrike;
    static int    default_numMiddleStrikes;
    static int    default_numTailStrikes;
    static double default_propMiddleStrikes;
    
    static double default_numMiddleStdDevs;
    static double default_numTailStdDevs;

    static int    default_maxStdDevs;
    static int    default_maxFineStdDevs;

    static bool   default_failProbs;

private:
    static IObject* defaultPDFParams();
    static void load(CClassSP& clazz);
};

typedef smartPtr<PDFParams> PDFParamsSP;
typedef smartConstPtr<PDFParams> PDFParamsConstSP;


///////////////////////////////
// BRACKETS
///////////////////////////////
class MARKET_DLL Brackets: public CObject {
public:
    static CClassConstSP const TYPE;
    
    Brackets(double                stdDevs,           // target stdDevs
             const PDFCalculator*  pdfCalc,           // pdfCalculator
             const DateTime&       volStartDate,      // the start date of the Vol
             const DateTimeArray&  datesTo,           // list of dates at which to bracket the stdDevs
             const DoubleArray&    fwds,              // fwds at datesTo
             double                fwdAtVolStart,     // the forward at the volStartDate
             const DoubleArray&    sqrtTotalVar,      // the ATM sqrt(variance) for datesTo
             int                   maxStdDevs,        // maximum number of standard deviations used for bracketing the desired stdDevs
             int                   maxFineStdDevs,    // maximum number of standard deviations used for finer bracketing of the desired stdDevs
             bool                  isVolFwdStarting); // is the Vol forward starting

    DoubleArray downStrikes;  //!< Strikes that bracket -stdDevs
    DoubleArray upStrikes;    //!< Strikes that bracket +stdDevs
    DoubleArray downStdDevs;  //!< The actual number of stdDevs used to bracket -stdDevs
    DoubleArray upStdDevs;    //!< The actual number of stdDevs used to bracket +stdDevs

private:
    void findBracket(const CLatticeDouble& probs,
                     const CLatticeDouble& strikes,
                     double                normalStdDev,
                     DoubleArray&          downStrike,
                     DoubleArray&          upStrike,
                     IntArray&             downVol,
                     IntArray&             upVol);
    
    Brackets();

    static void load(CClassSP& clazz);
    static IObject* defaultBrackets();
};

typedef smartPtr<Brackets> BracketsSP;
typedef smartConstPtr<Brackets> BracketsConstSP;


////////////////////////////////
// NUMBER OF STRIKES
////////////////////////////////
class MARKET_DLL StrikeNumbers: public CObject {
public:
    static CClassConstSP const TYPE;
    
    StrikeNumbers(BracketsSP midBrackets,                // the Brackets for the middle of the distribution
                  BracketsSP tailBrackets,               // the Brackets for the tail of the distribution
                  int        numMidStrikesPerLogStrike,  // number of strikes per logstrike used in the middle of the distribution
                  int        numMiddleStrikes,           // the minimum number of strikes that can be used in the middle
                  int        numTailStrikesPerLogStrike, // number of strikes per logstrike used in the tail of the distribution
                  int        numTailStrikes,             // the minimum number of strikes that can be used in the tail
                  double     propMiddleStrikes);         // the proportion of middle strikes that will be used in the first pass
    
    /** Returns the total number of strikes used */
    int totalNumStrikes(int iStep) const;
    
    /** Returns the total number of strikes used for the first pass */
    int totalNumStrikesFP(int iStep) const;
    
    IntArray nbDownTail;    //!< Number of strikes used in the left tail
    IntArray nbMid;         //!< Number of strikes used in the middle
    IntArray nbUpTail;      //!< Number of strikes used in the right tail
    
    IntArray nbMidFP;       //!< Number of strikes used in the middle for the first pass

private:
    StrikeNumbers();
    static IObject* defaultStrikeNumbers();
    static void load(CClassSP& clazz);
};

typedef smartPtr<StrikeNumbers> StrikeNumbersSP;
typedef smartConstPtr<StrikeNumbers> StrikeNumbersConstSP;


//////////////////////////////////
// STRIKES PARTITION
//////////////////////////////////
class MARKET_DLL StrikesPartition: public CObject {
public:
    static CClassConstSP const TYPE;
    
    StrikesPartition(const PDFCalculator*  pdfCalc,             // pdfCalculator
                     const DateTime&       volStartDate,        // the start date of the vol
                     const DateTimeArray&  datesTo,             // list of dates
                     const DoubleArray&    fwds,                // list of fwds for datesTo
                     double                fwdAtVolStart,       // the fwd at the vol start date
                     const DoubleArray&    sqrtTotalVar,        // the ATM sqrt(variance) for datesTo
                     const PDFParams&      params,              // PDFParams
                     bool                  isVolFwdStarting);   // is vol forward starting
    
    // transient
    BracketsSP        midBrackets;      //!< The brackets for the middle of the distribution
    BracketsSP        tailBrackets;     //!< The brackets for the tail of the distribution
    StrikeNumbersSP   nbStrikes;        //!< A breakdown of the number of strikes used in every part of the distribution

protected:
    StrikesPartition();
    static IObject* defaultStrikesPartition();
    static void load(CClassSP& clazz);
};

typedef smartPtr<StrikesPartition> StrikesPartitionSP;
typedef smartConstPtr<StrikesPartition> StrikesPartitionConstSP;


//////////////////////////////////////////////////////////////
// LinearImpliedSampler
//////////////////////////////////////////////////////////////
class MARKET_DLL LinearImpliedSampler: public CObject {
public:
    static CClassConstSP const TYPE;
    friend class LinearImpliedSamplerHelper;

    // Constructor with basic inputs
    LinearImpliedSampler(const CAssetConstSP&        asset,         // asset
                         const VolRequestLNStrikeSP& volRequest,    // volRequest
                         const DateTime&             today,         // valueDate
                         const DateTime&             startDate,     // startDate (different from today for FwdStarting)
                         const DateTimeArray&        datesTo,       // list of dates
                         const PDFParams&            params,        // PDFParams
                         const PDFRequestSP&         pdfRequest);   // pdfRequest

    // Constructor with additional inputs: fwds and vols
    LinearImpliedSampler(const CAssetConstSP&        asset,         // asset
                         const VolRequestLNStrikeSP& volRequest,    // volRequest
                         const DateTime&             today,         // valueDate
                         const DateTime&             startDate,     // startDate (different from today for FwdStarting)
                         const DateTimeArray&        datesTo,       // list of dates
                         const PDFParams&            params,        // PDFParams
                         const PDFRequestSP&         pdfRequest,    // pdfRequest
                         const DoubleArraySP&        fwds,          // fwds at datesTo
                         const DoubleArraySP&        sqrtTotalVar); // ATM sqrt(variance) for datesTo

    // Constructor that takes everything as input
    LinearImpliedSampler(const CAssetConstSP&        asset,         //asset
                         const VolRequestLNStrikeSP& volRequest,    // volRequest
                         const DateTime&             today,         // valueDate
                         const DateTime&             startDate,     // startDate (different from today for FwdStarting)
                         const DateTimeArray&        datesTo,       // list of dates
                         const PDFParams&            params,        // PDFParams
                         const PDFRequestSP&         pdfRequest,    // pdfRequest
                         const DoubleArraySP&        fwds,          // fwds at datesTo
                         const DoubleArraySP&        sqrtTotalVar,  // ATM sqrt(variance) for datesTo
                         const CLatticeDoubleSP&     strikes,       // tabulated strikes
                         const CLatticeDoubleSP&     probs,         // tabulated probs
                         const StrikesPartitionSP&   partition);    // strikes partition

    // Create a LinearImpliedSampler from an existing sampler, but with 
    // explicit partition & strikes
    LinearImpliedSampler(const LinearImpliedSampler* base,
                         const CLatticeDoubleSP&     strikes,       // tabulated strikes
                         const StrikesPartitionSP&   partition);    // strikes partition
                         
    /** Gets the strikes partition of the distribution */
    StrikesPartitionSP         getPartition();
    
    /** Gets the tabulated strikes */
    CLatticeDoubleSP           getStrikes();
    
    /** Gets the tabulated probs */
    CLatticeDoubleSP           getProbs();
    
    /** Gets an array of ImpliedSamples */
    LinearImpliedSampleArraySP getImpliedSamples();

    /** Gets array sqrtTotalVar */
    DoubleArraySP getSqrtTotalVar(double interpFwdSpot);

    /** Gets array fwds */
    DoubleArraySP getFwd();

    virtual void validatePop2Object();

protected:
    /** Computes implied probs using a PDFCalculator and smooths them */
    void computeImpliedProbs(const DateTimeArray&  datesTo,
                             const CLatticeDouble& strikes,
                             CLatticeDouble&       probs) const;

    /** Algorithm that smooths ill-defined distributions */
    void massageDistribution(const CSliceDouble& sliceStrikes,   // the strikes
                             CSliceDouble&       sliceDist,      // the distributions
                             bool                midPercentile   /* whether to start from the midStrike or
                                                                    the midProbability */ ) const;

private:
    // Basic constructor
    LinearImpliedSampler();
    
    // required from empty shell
    CAssetConstSP               asset;
    VolRequestLNStrikeSP        volRequest;
    DateTime                    today;
    DateTime                    startDate;
    DateTimeArray               datesTo;
    PDFParams                   params;
    PDFRequestSP                pdfRequest;

    // Optional
    DoubleArraySP               fwds;
    DoubleArraySP               sqrtTotalVar;
    CLatticeDoubleSP            strikes; // $unregistered
    CLatticeDoubleSP            probs; // $unregistered
    StrikesPartitionSP          partition;

    LinearImpliedSampleArraySP  samples;
    
    // Transient
    string                      assetName;
    int                         numDates;
    DateTime                    volStartDate;
    double                      fwdAtVolStart;
    bool                        isVolFwdStarting;

    PDFCalculatorSP             pdfCalc;
};

typedef smartPtr<LinearImpliedSampler> LinearImpliedSamplerSP;
typedef smartConstPtr<LinearImpliedSampler> LinearImpliedSamplerConstSP;


/** Implied Mapping Interface for mapping a spot to driver
    and vice versa */
class MARKET_DLL ImpliedMapping: virtual public IObject {
public:
    static CClassConstSP const TYPE;

    /** Maps the driver to spot */
    virtual double spotFromDriver(double driver) const = 0;

    /** Maps the spot to driver */
    virtual double driverFromSpot(double spot) const = 0;

    /** Returns whether this is an implied sample of a 
        return or of an absolute price */
    virtual bool isRelativeReturn() const = 0;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

typedef smartPtr<ImpliedMapping> ImpliedMappingSP;
typedef smartConstPtr<ImpliedMapping> ImpliedMappingConstSP;

typedef array<ImpliedMappingSP, ImpliedMapping> ImpliedMappingArray;
typedef smartPtr<ImpliedMappingArray> ImpliedMappingArraySP;
typedef smartConstPtr<ImpliedMappingArray> ImpliedMappingArrayConstSP;

/** Gaussian Implied Mapping i.e. driver is gaussian */
class MARKET_DLL GaussianImpliedMapping: public CObject,
                              virtual public ImpliedMapping {
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    GaussianImpliedMapping(LinearImpliedSampleSP marketSample,
                           double vol);

    /** Maps the driver to spot */
    virtual double spotFromDriver(double driver) const;

    /** Maps the spot to driver */
    virtual double driverFromSpot(double spot) const;

    /** Returns whether this is an implied sample of a 
        return or of an absolute price */
    virtual bool isRelativeReturn() const;

private:
    /** Default constructor */
    GaussianImpliedMapping();

    /** Basic validation on parameters */
    virtual void validatePop2Object();

    /** For reflection */
    static IObject* defaultGaussianImpliedMapping();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    LinearImpliedSampleSP marketSample;      //!< Mapping from uniform to spot
    double                vol;               //!< Driver's vol for that date
};

typedef smartPtr<GaussianImpliedMapping> GaussianImpliedMappingSP;
typedef smartConstPtr<GaussianImpliedMapping> GaussianImpliedMappingConstSP;


/** Gaussian Implied Mapping i.e. driver is gaussian but the 
    cumulative distribution is computed using an interpolant 
    for speed in MCImplied */
class MARKET_DLL GaussianInterpImpliedMapping: public CObject,
                                    virtual public ImpliedMapping {
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    GaussianInterpImpliedMapping(EquidistantLinearInterpolantNonVirtualSP tabulatedNorm,
                                 LinearImpliedSampleSP marketSample,
                                 double vol);

    /** Maps the driver to spot */
    virtual double spotFromDriver(double driver) const;

    /** Maps the spot to driver */
    virtual double driverFromSpot(double spot) const;

    /** Override close method to take shallow copy of tabulated norm */
    virtual IObject* clone() const;

    /** Returns whether this is an implied sample of a 
        return or of an absolute price */
    virtual bool isRelativeReturn() const;

private:
    /** Default constructor */
    GaussianInterpImpliedMapping();

    /** Basic validation on parameters */
    virtual void validatePop2Object();

    /** For reflection */
    static IObject* defaultGaussianInterpImpliedMapping();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    EquidistantLinearInterpolantNonVirtualSP tabulatedNorm;     //!< Tabulated standard normal distribution
    LinearImpliedSampleSP                    marketSample;      //!< Mapping from uniform to spot
    double                                   vol;               //!< Driver's vol for that date
};

typedef smartPtr<GaussianInterpImpliedMapping> GaussianInterpImpliedMappingSP;
typedef smartConstPtr<GaussianInterpImpliedMapping> GaussianInterpImpliedMappingConstSP;


/////////////////////////////////////////////
// Degenerate Implied Mapping
/////////////////////////////////////////////
class MARKET_DLL DegenerateImpliedMapping: public CObject,
                                virtual public ImpliedMapping {
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    DegenerateImpliedMapping(DegenerateImpliedSampleSP marketSample);

    /** Maps the driver to spot */
    virtual double spotFromDriver(double driver) const;

    /** Maps the spot to driver */
    virtual double driverFromSpot(double spot) const;

    /** Returns whether this is an implied sample of a 
        return or of an absolute price */
    virtual bool isRelativeReturn() const;

    /** Basic validation on parameters */
    virtual void validatePop2Object();

private:
    /** Default constructor */
    DegenerateImpliedMapping();

    /** For reflection */
    static IObject* defaultDegenerateImpliedMapping();
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    
    DegenerateImpliedSampleSP marketSample;     //!< Mapping from uniform to spot
    double                    spotFwd;          //!< Spot forward
};

typedef smartPtr<DegenerateImpliedSample> DegenerateImpliedSampleSP;
typedef smartConstPtr<DegenerateImpliedSample> DegenerateImpliedSampleConstSP;


DRLIB_END_NAMESPACE

#endif
