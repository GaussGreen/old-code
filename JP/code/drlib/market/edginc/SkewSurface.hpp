//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SkewSurface.hpp
//
//   Description : Surface of skews (fast and normal) indexed by strikes and maturities
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SKEW_SURFACE_HPP
#define EDR_SKEW_SURFACE_HPP

//#include "edginc/DateTime.hpp"
#include "edginc/BetaSkewParallel.hpp"
#include "edginc/BetaSkewPointwiseTweak.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/QuasiContractualBaseCorrelation.hpp"
#include "edginc/TimePoint2D.hpp"

DRLIB_BEGIN_NAMESPACE

// Forward declaration of SkewSurface
class SkewSurface;

typedef smartPtr<SkewSurface> SkewSurfaceSP;
typedef smartConstPtr<SkewSurface> SkewSurfaceConstSP;

typedef array<SkewSurfaceSP, SkewSurface> SkewSurfaceArray;
typedef smartPtr<SkewSurfaceArray> SkewSurfaceArraySP;

// Support for wrapper class
typedef MarketWrapper<SkewSurface> SkewSurfaceWrapper;

/** 
 * Surface of skews (fast and normal) indexed by strikes and maturities.
 * */
class MARKET_DLL SkewSurface : 
    public MarketObject,
    public virtual BetaSkewPointwiseTweak::IShift,
    virtual public QuasiContractualBaseCorrelation::IShift,
    public virtual Calibrator::IAdjustable
{
public:

    /** Type for 1 dimension linear interpolation */
    static const string LINEAR_1D_INTERPOLATION;

    /** Type for 1 dimension flat interpolation */
    static const string FLAT_1D_INTERPOLATION;
    
    /** Type for 1 dimension splice interpolation */
    static const string SPLICE_1D_INTERPOLATION;
    
    /** Type for 1 dimension spline interpolation */
    static const string SPLINE_1D_INTERPOLATION;

    /**
     * Types of skew :
     * - Normal skew (computed using usual convolution algorithm) : used
     * when number of names < threshold 
     * 
     * - Fast skew (computed using fast convolution algorithm) : used
     * when number of names >= threshold 
     * */
    enum SkewType {
        NORMAL,
        FAST
    };


    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    ~SkewSurface();
    
    /** Called immediately after object constructed */
    void validatePop2Object();

    /**
     * Returns interpolated skew of given type (normal, fast...)
     * at point (strike, date)
     * */
    double getSkew(double strike, const DateTime& date, SkewType skewType) const;

	/** get arrays of skews and strikes defined at particular date in skewSurface 
		date must lie on skew surface timeline otherwise an error is thrown */
	void getStrikesAndSkews(
		DoubleArray & strikes,	// (O) strikes to populate
		DoubleArray & skews,	// (O) skews to populate
		const DateTime & date,
		SkewType skewType = SkewSurface::NORMAL) const;

    /**
     * Returns "neighbours" of the given point P{strike, date}
     * i.e. points needed to interpolate at P
     * */
    BetaSkewGridPointSet getNeighbours(double strike, 
                                       const DateTime& date,
                                       const int numberOfNeighbours) const;

	/** get the maturities */
	DateTimeArrayConstSP getMaturities() const { return maturities; }

	/** get maturity interpolation type */
	string getMatInterpType() const { return matInterpolationType; }

    /** Set a name for this function (useful for tweaking) */
    void setName(string name);
    
    /**
     * Returns the name of this object.
     * [Implements Calibrator::IAdjustable]
     * */
    virtual string getName() const;

    /** Returns all points defining the surface */
    BetaSkewGridPointArrayConstSP getAllPoints() const;

    /** Implementation of BetaSkewPointwiseTweak::IShift */
    virtual string sensName(BetaSkewPointwiseTweak* shift) const;
    virtual bool sensShift(BetaSkewPointwiseTweak* shift);

    /** Implementation of QuasiContractualBaseCorrelation::IShift */
    virtual bool sensShift(QuasiContractualBaseCorrelation* shift);

    /** Tweaking methods */
    void sensShift(BetaSkewParallel* shift);
    void sensRestore(BetaSkewParallel* shift);
    
    /** give acces to bootstrpper */
    friend class SkewSurfaceBootstrapper;
    friend class BSLPBootstrapper;
    
private:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Only build instances of that class using reflection */
    SkewSurface();
    
    /** Default constructor */
    static IObject* defaultConstructor();
    
    /**
     * This method returns the strikes and skews for a given maturity index.
     * This method assumes that the list of points is sorted by increasing
     * strikes.
     */
    void findMaturitySkew(
        SkewType skewType,
        int idx,
        int* nbStrikeMat,
        double** strikeMat,
        double** skewMat) const;

    // ----------------
    // MANDATORY FIELDS
    // ----------------
    
    // The fields strikes, maturities, skews and fastSkews are
    // arrays of the same length N describing the list of
    // available points (strike, maturity, skew, skew fast).
    // We use a 'flat' structure here because the skew matrix can 
    // have undefined values for some (strike, maturity) points.

    /** Array of strikes, length N */
    CDoubleArraySP strikes;

    /** Array of maturities, length N */
    DateTimeArraySP maturities;
    
    /** Array of  skews (normal in beta skew case), length N */
    CDoubleArraySP skews;       
    
    // ---------------
    // OPTIONAL FIELDS
    // ---------------
	
	/** Array of fast beta skews, length N */       
	CDoubleArraySP fastSkews;

    /** 
     * Type of interpolation in the strike dimension 
     * (should be the same as the one used to compute skews)
     * */
    string strikeInterpolationType;
    
    /**
     * Type of interpolation in the time dimension 
     * (should be the same as the one used to compute skews)
     * */
    string matInterpolationType;
    
    /** Name (for tweaking)*/
    string name;



};

DRLIB_END_NAMESPACE

#endif
