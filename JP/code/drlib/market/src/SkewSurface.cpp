//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SkewSurface.cpp
//
//   Description : 
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SkewSurface.hpp"
#include "edginc/Spline.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

#define TINY 3e-15


// -------------------------------------
// Declarations of interpolation methods
// -------------------------------------

/** linear interpolation between 2 time points */
double LinearTimePointInterp(
    const DateTime& prevMat,
    double betaLowPrevMat,
    const DateTime& nextMat,
    double betaLowNextMat,
    const DateTime& t);

/** flat interpolation between 2 time points */
double FlatTimePointInterp(
    const DateTime& prevMat,
    double betaLowPrevMat,
    const DateTime& nextMat,
    double betaLowNextMat,
    const DateTime& t);
    
/**
 * This function interpolates in the strike dimension.
 * It uses splice interpolation.
 */
void SpliceStrikeInterp(
        int             n,      /* (I) nb of points in the input arrays */
        const double   *xData,  /* (I) x values                         */
        const double   *yData,  /* (I) y values                         */
        double          x,      /* (I) value of x that musst be interp  */
        double         *y);     /* (O) interpolated y value             */

/**
 * This function interpolates in the strike dimension.
 * It uses spline interpolation.
 */
void SplineStrikeInterp(
        int             n,      /* (I) nb of points in the input arrays */
        const double   *xData,  /* (I) x values                         */
        const double   *yData,  /* (I) y values                         */
        double          x,      /* (I) value of x that musst be interp   */
        double         *y);     /* (O) interpolated y value             */

/**
* This function interpolates in the strike dimension.
* It uses linear interpolation and flat extrapolation.
*/
void LinearStrikeInterp(
						int             n,      /* (I) nb of points in the input arrays */
						const double   *xData,  /* (I) x values (need to be increasing) */
						const double   *yData,  /* (I) y values                         */
						double          x,      /* (I) value of x that must be interp   */
						double         *y);      /* (O) interpolated y value             */

/**
* This function interpolates in the strike dimension.
* It uses flat interpolation and flat extrapolation.
*/
void FlatStrikeInterp(
					  int             n,      /* (I) nb of points in the input arrays */
					  const double   *xData,  /* (I) x values (need to be increasing) */
					  const double   *yData,  /* (I) y values                         */
					  double          x,      /* (I) value of x that must be interp   */
					  double         *y);      /* (O) interpolated y value             */


// --------------
// Implementation
// --------------

/** Only build instances of that class using reflection */
SkewSurface::SkewSurface() :
    MarketObject(TYPE),
    strikes(new DoubleArray(0)),
    maturities(new DateTimeArray(0)),
    skews(new DoubleArray(0)),
    fastSkews(NULL),
    strikeInterpolationType(SPLINE_1D_INTERPOLATION),
    matInterpolationType(LINEAR_1D_INTERPOLATION),
    name("")
	{}

/** Default constructor */
IObject* SkewSurface::defaultConstructor() {
    return new SkewSurface();
}

/** Destructor */
SkewSurface::~SkewSurface() {}

/**
 * Returns interpolated skew of given type (normal, fast...)
 * at point (strike, date)
 * */
double SkewSurface::getSkew(double strike, const DateTime& date, SkewType skewType) const {
    static const string method("SkewSurface::getSkew");
 
    int idxPrev = date.findLower(*maturities);
    DateTime prevMat;
    if (idxPrev != -1) {
        prevMat = (*maturities)[idxPrev];
    }

    int idxNext = date.findUpper(*maturities);
    DateTime nextMat;
    if (idxNext != maturities->size()) {
        nextMat = (*maturities)[idxNext];
    }

    int nbStrikePrevMat = 0;
    int nbStrikeNextMat = 0;
    double *strikePrevMat   = NULL;
    double *strikeNextMat   = NULL;
    double *skewPrevMat     = NULL;
    double *skewNextMat     = NULL;
    
    if (idxPrev != -1) {
        findMaturitySkew(
            skewType, idxPrev, &nbStrikePrevMat, &strikePrevMat, &skewPrevMat);
    }

    if (idxNext != maturities->size()) {
        findMaturitySkew(
            skewType, idxNext, &nbStrikeNextMat, &strikeNextMat, &skewNextMat);
    }

    if (nbStrikePrevMat == 0 && nbStrikeNextMat == 0) {
        throw ModelException(method,
                "Failed : no strikes to interpolate");
    }

    /* 
     * we extrapolate flatly before the 1st mat and after the last mat
     * so, if a previous maturity for the index does not exist,
     * we set the skew for previous existing maturity equal to the 
     * skew for the next existing maturity.
     * We do the similar thing if a next maturity does not exists.
     */
    if (nbStrikePrevMat == 0)
    {
        /* prev mat skews is equal to next mat skew */
        nbStrikePrevMat = nbStrikeNextMat;
        strikePrevMat   = strikeNextMat;
        skewPrevMat     = skewNextMat;
    }

    if (nbStrikeNextMat == 0)
    {
        /* prev mat skews is equal to next mat skew */
        nbStrikeNextMat = nbStrikePrevMat;
        strikeNextMat   = strikePrevMat;
        skewNextMat     = skewPrevMat;
    }

    double betaPrevMat = 0.;
    double betaNextMat = 0.;

    if (strikeInterpolationType == SPLICE_1D_INTERPOLATION) {
        /* splice interpolation for strike on both skews */
        SpliceStrikeInterp(
            nbStrikePrevMat, strikePrevMat, skewPrevMat, strike, &betaPrevMat);
        SpliceStrikeInterp(
            nbStrikeNextMat, strikeNextMat, skewNextMat, strike, &betaNextMat);
    } else if (strikeInterpolationType == SPLINE_1D_INTERPOLATION) {
        /* spline interpolation for strike on both skews */
        SplineStrikeInterp(
            nbStrikePrevMat, strikePrevMat, skewPrevMat, strike, &betaPrevMat);
        SplineStrikeInterp(
            nbStrikeNextMat, strikeNextMat, skewNextMat, strike, &betaNextMat);
    } 
	else if (strikeInterpolationType == LINEAR_1D_INTERPOLATION) {
		/* Linear interpolation for strike on both skews */
		LinearStrikeInterp(
			nbStrikePrevMat, strikePrevMat, skewPrevMat, strike, &betaPrevMat);
		LinearStrikeInterp(
			nbStrikeNextMat, strikeNextMat, skewNextMat, strike, &betaNextMat);
	}
	else if (strikeInterpolationType == FLAT_1D_INTERPOLATION) {
		/* Linear interpolation for strike on both skews */
		FlatStrikeInterp(
			nbStrikePrevMat, strikePrevMat, skewPrevMat, strike, &betaPrevMat);
		FlatStrikeInterp(
			nbStrikeNextMat, strikeNextMat, skewNextMat, strike, &betaNextMat);
	}
	else {
        throw ModelException(method,
                strikeInterpolationType +
                " : strike interpolation type not supported");
    }

    double result = 0.0;
    if (matInterpolationType == LINEAR_1D_INTERPOLATION) {
        /* linear interpolation in time */
        result = LinearTimePointInterp(
            prevMat, betaPrevMat, nextMat, betaNextMat, date);
    } else if (matInterpolationType == FLAT_1D_INTERPOLATION) {
        /* flat interpolation in time          */
        result = FlatTimePointInterp(
            prevMat, betaPrevMat, nextMat, betaNextMat, date);
    } else {
        throw ModelException(method,
                matInterpolationType +
                " : time interpolation type not supported");
    }

    return result;
}

/** get arrays of skews and strikes defined at particular date in skewSurface 
date must lie on skew surface timeline otherwise an error is thrown */
void SkewSurface::getStrikesAndSkews(
						DoubleArray & strikes,	// (O) strikes to populate
						DoubleArray & skews,	// (O) skews to populate
						const DateTime & date,
						SkewType skewType) const
{
	static const string method = "SkewSurface::getStrikesAndSkews";
	
	//index of first occurance of date in maturities. Fails if not found
	int dateIdx = date.find(*maturities); 

	int nbStrikesMat = 0;
	double *strikesMat   = NULL;
	double *skewsMat     = NULL;
	
	findMaturitySkew(
			skewType, dateIdx, &nbStrikesMat, &strikesMat, &skewsMat);
	

	// resize output arrays
	strikes.resize(nbStrikesMat);
	skews.resize(nbStrikesMat);

	// populate
	int i = 0;
	for(; i < nbStrikesMat; i++ )
	{
		strikes[i] = strikesMat[i];
		skews[i] = skewsMat[i];
	}
	

}
/**
 * Returns "neighbours" of the given point P{strike, date}
 * i.e. points needed to interpolate at P
 * */
BetaSkewGridPointSet SkewSurface::getNeighbours(
    double strike, 
    const DateTime& date,
    const int numberOfNeighbours) const 
{
    static const string method("SkewSurface::getNeighbours");
    
    int idxPrev = date.findLower(*maturities);
    int idxNext = date.findUpper(*maturities);

    DateTime prevMat;
    if (idxPrev != -1) {
        prevMat = (*maturities)[idxPrev];
    } else {
        prevMat = (*maturities)[idxNext];
    }

    DateTime nextMat;
    if (idxNext != maturities->size()) {
        nextMat = (*maturities)[idxNext];
    } else {
        nextMat = (*maturities)[idxPrev];
    }

    int nbStrikePrevMat = 0;
    int nbStrikeNextMat = 0;
    double *strikePrevMat   = NULL;
    double *strikeNextMat   = NULL;
    double *skewPrevMat     = NULL;
    double *skewNextMat     = NULL;
    
    if (idxPrev != -1) {
        findMaturitySkew(
            NORMAL, idxPrev, &nbStrikePrevMat, &strikePrevMat, &skewPrevMat);
    }

    if (idxNext != maturities->size()) {
        findMaturitySkew(
            NORMAL, idxNext, &nbStrikeNextMat, &strikeNextMat, &skewNextMat);
    }

    if (nbStrikePrevMat == 0 && nbStrikeNextMat == 0) {
        throw ModelException(method,
                "Failed : no strikes to interpolate");
    }

    /* 
     * we extrapolate flatly before the 1st mat and after the last mat
     * so, if a previous maturity for the index does not exist,
     * we set the skew for previous existing maturity equal to the 
     * skew for the next existing maturity.
     * We do the similar thing if a next maturity does not exists.
     */
    if (nbStrikePrevMat == 0) {
        /* prev mat skews is equal to next mat skew */
        nbStrikePrevMat = nbStrikeNextMat;
        strikePrevMat   = strikeNextMat;
        skewPrevMat     = skewNextMat;
    }

    if (nbStrikeNextMat == 0) {
        /* next mat skews is equal to prev mat skew */
        nbStrikeNextMat = nbStrikePrevMat;
        strikeNextMat   = strikePrevMat;
        skewNextMat     = skewPrevMat;
    }
    
    /* Finds immediate neighbours of "strike" in strikePrevMat*/
    bool found = false;
    int idx = 0;
    for(; idx < nbStrikePrevMat && !found;) {
        if (strike < strikePrevMat[idx]) {
            found = true;
        } else {
            idx++;    
        }
    }
    int idxStrikePrevMatInf, idxStrikePrevMatSup;
    if (!found) {
        idxStrikePrevMatInf = nbStrikePrevMat-1;
        idxStrikePrevMatSup = nbStrikePrevMat-1;
    } else {
        idxStrikePrevMatInf = (idx == 0) ? 0 : (idx-1);
        idxStrikePrevMatSup = idx;     
    }
    
    /* Finds immediate neighbours of "strike" in strikeNextMat*/
    found = false;
    idx = 0;
    for(; idx<nbStrikeNextMat && !found;) {
        if (strike < strikeNextMat[idx]) {
            found = true;
        } else {
            idx++;
        }    
    }
    int idxStrikeNextMatInf, idxStrikeNextMatSup;
    if (!found) {
        idxStrikeNextMatInf = nbStrikeNextMat-1;
        idxStrikeNextMatSup = nbStrikeNextMat-1;
    } else {
        idxStrikeNextMatInf = (idx == 0) ? 0 : (idx-1);
        idxStrikeNextMatSup = idx;     
    }
    
    // HACK:
    // if first strike is negative (say -99%), we will actually
    // use beta corresponding to -99% to 
    // linearly interpolate between 0% and [first positive strike]%
    // so in that case the nearest previous strike of a point between 0%
    // and [first positive strike]% is actually -99%
    if (strikePrevMat[0] < 0.0) {
        if (idxStrikePrevMatInf == 1) {
            idxStrikePrevMatInf = 0;
        }
    }
    if (strikeNextMat[0] < 0.0) {
        if (idxStrikeNextMatInf == 1) {
            idxStrikeNextMatInf = 0;
        }
    }

    /* Builds the result set */
    BetaSkewGridPointSet result;
    for (int i=0; i < numberOfNeighbours; ++i) {
        if (idxStrikePrevMatInf-i >= 0) {
            result.insert(BetaSkewGridPointSP(
                new BetaSkewGridPoint(prevMat, strikePrevMat[idxStrikePrevMatInf-i])));
        }
        if (idxStrikePrevMatSup+i < nbStrikePrevMat) {
            result.insert(BetaSkewGridPointSP(
                new BetaSkewGridPoint(prevMat, strikePrevMat[idxStrikePrevMatSup+i])));
        }


        if (idxStrikeNextMatInf-i >= 0) {
            result.insert(BetaSkewGridPointSP(
                new BetaSkewGridPoint(nextMat, strikeNextMat[idxStrikeNextMatInf-i])));
        }
        if (idxStrikeNextMatSup+i < nbStrikeNextMat) {
          result.insert(BetaSkewGridPointSP(
                new BetaSkewGridPoint(nextMat, strikeNextMat[idxStrikeNextMatSup+i])));
        }   
    }
    ASSERT((int)result.size() <= numberOfNeighbours*4);

    return result;
}


/** Called immediately after object constructed */
void SkewSurface::validatePop2Object() {
    static const string method("SkewSurface::validatePop2Object");        
    try{
        int nbMat = maturities->size();           
        
        // checks that maturities, strikes, skews and fastSkews have same length
		if(!fastSkews)
		{
			if (nbMat != strikes->size() ||
				nbMat != skews->size() ) {
				throw ModelException(method,
					"maturities, strikes, and skews "
					"don't have the same length");            
			} 
		}
		else
		{
			if (nbMat != strikes->size() ||
				nbMat != skews->size() ||
				nbMat != fastSkews->size()) {
					throw ModelException(method,
						"maturities, strikes, skews and fastSkews "
						"don't have the same length");            
				} 
		}
		
       
        // checks that maturities are increasing
        DateTime::ensureIncreasing(
            *maturities,
            "skew maturities",
            true);
        
        // checks that strikes are increasing for a given maturity
        for(int i=1; i<nbMat; i++) {
            if ((*strikes)[i-1] >= (*strikes)[i] && 
                (*maturities)[i-1] == (*maturities)[i]) {
                throw ModelException(method,
                    "Strikes not increasing for maturity " +
                    (*maturities)[i].toString());
            }
        }

		

    } catch (exception& e){
        throw ModelException(e, "SkewSurface::validatePop2Object");
    }
}

/** Tweaking methods */
void SkewSurface::sensShift(BetaSkewParallel* shift)
{
    for (int i=0; i<strikes->size(); i++)
    {
        //either skews or fastSkews are in use
        //so move both together
        (*skews)[i] += shift->getShiftSize();
        (*fastSkews)[i] += shift->getShiftSize();
    }
}

void SkewSurface::sensRestore(BetaSkewParallel* shift)
{
    for (int i=0; i<strikes->size(); i++)
    {
        //either skews or fastSkews are in use
        //so move both together
        (*skews)[i] -= shift->getShiftSize();
        if(!!fastSkews) (*fastSkews)[i] -= shift->getShiftSize();
    }
}

/** Invoked when Class is 'loaded' */
void SkewSurface::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SkewSurface, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(BetaSkewPointwiseTweak::IShift);
    IMPLEMENTS(QuasiContractualBaseCorrelation::IShift);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultConstructor);
    
    FIELD(strikes, "Array of strikes, length N");
    FIELD(maturities, "Array of maturities, length N");
    FIELD(skews, "Array of skews (normal), length N");
    FIELD(fastSkews, "Array of fast skews (used for beta skew), length N");
	FIELD_MAKE_OPTIONAL(fastSkews);

    FIELD(strikeInterpolationType,
        "Type of interpolation in the strike dimension "
        "(should be the same as the one used to compute skews)");
    FIELD_MAKE_OPTIONAL(strikeInterpolationType);
    FIELD(matInterpolationType,
        "Type of interpolation in the time dimension "
        "(should be the same as the one used to compute skews)");
    FIELD_MAKE_OPTIONAL(matInterpolationType);
    FIELD(name, "name");
    FIELD_MAKE_OPTIONAL(name);

    
    // Register "skew" as a field that can be calibrated
    Calibrator::IAdjustable::registerField(
        clazz,
        "skews",
       new Range(OpenBoundary(0.0),  OpenBoundary(1.0)));

    // Register "fastSkew" as a field that can be calibrated
    Calibrator::IAdjustable::registerField(
        clazz,
        "fastSkews",
        new Range(OpenBoundary(0.0),  OpenBoundary(1.0)));
}

/** Type for 1 dimension linear interpolation */
const string SkewSurface::LINEAR_1D_INTERPOLATION = "LINEAR";

/** Type for 1 dimension flat interpolation */
const string SkewSurface::FLAT_1D_INTERPOLATION = "FLAT";
    
/** Type for 1 dimension splice interpolation */
const string SkewSurface::SPLICE_1D_INTERPOLATION = "SPLICE";
    
/** Type for 1 dimension spline interpolation */
const string SkewSurface::SPLINE_1D_INTERPOLATION = "SPLINE";




/** TYPE (for reflection) */
CClassConstSP const SkewSurface::TYPE = CClass::registerClassLoadMethod(
    "SkewSurface", typeid(SkewSurface), SkewSurface::load);

/** TYPE for SkewSurfaceArray */
DEFINE_TEMPLATE_TYPE(SkewSurfaceArray);

/** TYPE for SkewSurfaceWrapper */
DEFINE_TEMPLATE_TYPE(SkewSurfaceWrapper);

/**
 * This method returns the strikes and skews for a given maturity index.
 * This method assumes that the list of points is sorted by increasing
 * strikes.
 */
void SkewSurface::findMaturitySkew(
    SkewType skewType,
    int idx,
    int* nbStrikeMat,
    double** strikeMat,
    double** skewMat) const
{
    static const string method("SkewSurface::findMaturitySkew");
    
    // Finds idxLow such that 
    // (*maturities)[idxLow+1] == *maturities)[idx]
    // and
    // (*maturities)[idxLow] < *maturities)[idx] (or idxLow == -1)
    int idxLow = idx;
    for(; (idxLow < 0) ? (false) : ((*maturities)[idxLow] == (*maturities)[idx]); idxLow--) {}
    
    int nbMat = maturities->size();
    // Finds idxHigh such that 
    // (*maturities)[idxHigh-1] == *maturities)[idx]
    // and
    // (*maturities)[idxHigh] > *maturities)[idx] (or idxHigh >= nbMat)
    int idxHigh = idx;
    for(; (idxHigh >= nbMat) ? (false) : ((*maturities)[idxHigh] == (*maturities)[idx]); idxHigh++) {}
    
    
    *nbStrikeMat = idxHigh - idxLow - 1;
    *strikeMat = &(*strikes)[idxLow+1];
    switch (skewType) {
        case NORMAL:
            *skewMat = &(*skews)[idxLow+1];
            break;
        case FAST:
			if(!fastSkews)
			{
				throw ModelException("skewType = FAST but fastSkews are NULL", method);
			}
            *skewMat = &(*fastSkews)[idxLow+1];
            break;
    }
}

/** Set a name for this function (useful for tweaking) */
void SkewSurface::setName(string name) {
    this->name = name;
}

/**
 * Returns the name of this object.
 * [Implements Calibrator::IAdjustable]
 * */
string SkewSurface::getName() const {
	return name;
}

/** Returns all points defining the surface */
BetaSkewGridPointArrayConstSP SkewSurface::getAllPoints() const {
    BetaSkewGridPointArraySP points(new BetaSkewGridPointArray(strikes->size()));
    for(int i=0; i<strikes->size();i++) {
        (*points)[i] = BetaSkewGridPointSP(
            new BetaSkewGridPoint((*maturities)[i], (*strikes)[i]));
    }
    return points;
}

/** Implementation of BetaSkewPointwiseTweak::IShift */
string SkewSurface::sensName(BetaSkewPointwiseTweak* shift) const {
    return name;
}

/** Implementation of BetaSkewPointwiseTweak::IShift */
bool SkewSurface::sensShift(BetaSkewPointwiseTweak* shift) {

    static const string method = "SkewSurface::sensShift(BetaSkewPointwiseTweak* shift)";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)){
            // find where to shift
            BetaSkewGridPoint pointToShift = shift->getCurrentPointToTweak();

            // naive search
            // (TODO: could be improved by using a stl set internally in SkewSurface)
            DateTime date = pointToShift.getMaturity();
            double strike = pointToShift.getStrike();
            bool found = false;
            for(int  i=0;i<strikes->size() && !found;i++) {
                if ((*strikes)[i] == strike && (*maturities)[i] == date) {
                    found = true;
                    (*skews)[i] += shiftSize;
					if(!!fastSkews) (*fastSkews)[i] += shiftSize;
                }
            }
            
            if (!found) {
                throw ModelException(method,
                     "couldn't find " + 
                    date.toString() + " maturity and" +
                    Format::toString(strike) + " strike in " + 
                    name);
            }
        }
        return true;//false; // none of our components has a skew type sensitivity
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Implementation of QuasiContractualBaseCorrelation::IShift */
bool SkewSurface::sensShift(QuasiContractualBaseCorrelation* shift) {

    static const string method = 
        "SkewSurface::sensShift(QuasiContractualBaseCorrelation* shift)";

    try {
        double skewLevel = shift->getSkewLevel();

        for (int i=0; i<strikes->size(); ++i) {
            (*skews)[i]     = skewLevel;
            if(!!fastSkews) (*fastSkews)[i] = skewLevel;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

    return false;
}


/*****  linterp  ************************************************************/
/*
*       Linear interpolation between 2 points.
*/
void linterp(
    double x, 
    double *y, 
    double x0, 
    double x1, 
    double y0, 
    double y1)
{
    double a,b,l0,l1;

    /* basic check */
    if (fabs(x0-x1)<TINY) {
        throw ModelException("linterp",
                "x inputs supplied for interpolation are the same.");
    }

    if ((x == x0) || (y0 == y1)) {
        *y = y0;
    } else if (x == x1) {
        *y = y1;
    } else {
        a=x-x1;
        b=x0-x1;
        l0=a/b;

        a=x-x0;
        b=x1-x0;
        l1=a/b;

        *y=l0*y0+l1*y1;
    }
}  /* linterp */

/** linear interpolation between 2 time points */
double LinearTimePointInterp(
    const DateTime& prevMat,
    double betaLowPrevMat,
    const DateTime& nextMat,
    double betaLowNextMat,
    const DateTime& t)
{
    long dcc1;
    long dcc;

    if (prevMat.empty()) {
        /* return the skew for the 1st maturity */
        return betaLowNextMat;
    }

    if (nextMat.empty()) {
        /* return the skew for the 1st maturity */
        return betaLowPrevMat;
    }
    
    dcc1 = t.daysDiff(prevMat);
    dcc  = nextMat.daysDiff(prevMat);

    /* check dcc to avoid division by zero */
    if (dcc < TINY) {
        return betaLowNextMat;
    }
    
    return betaLowPrevMat + (double) dcc1 / (double) dcc * (betaLowNextMat - betaLowPrevMat);
}

/** flat interpolation between 2 time points */
double FlatTimePointInterp(
    const DateTime& prevMat,
    double betaLowPrevMat,
    const DateTime& nextMat,
    double betaLowNextMat,
    const DateTime& t)
{
    if (prevMat.empty()) {
        /* return the skew for the 1st maturity */
        return betaLowNextMat;
    }

    if (nextMat.empty()) {
        /* return the skew for the 1st maturity */
        return betaLowPrevMat;
    }
    
    return betaLowNextMat;
}

/*

    NAME:             spliced.c 

    SUMMARY:        interpolation method using spliced parabolas. 

     AUTHOR:            G. A. Gatoff.

    HISTORY:        created on 23rd August 1995.
                    Added slope and curve outputs on March 4th, 1996.
                    May 2003: merged with fourWeights.c and added getSpliceWeights

    RETURNS:        void.

    CAUTION:        None.

    INPUTS:

            wantY        -     1 if you want y(x) interpolated, else 0.
            wantSlope    -     1 if you want y'(x) interpolated, else 0.
            wantCurve    -     1 if you want y''(x) interpolated, else 0.
            n            -    number of points in data
            xData          -    array of size n (x_1,... x_n)
            yData        -    array of size n (y_1,... y_n)
            gamma1        -    slope for x < x1
            gamman        -    slope for x >= xn
            npoints     -    number of interpolted points.
            x            -    x points for which y needs interpolating.
    OUTPUTS:
            y            -    interpolated y points
            slope         -    interploated slopes (deltas)
            curve         -   interpolated curves (gammas)

*/

void spliced(int wantY, int wantSlope, int wantCurve, int n, const double *xData, const double *yData,
             double gamma1, double gamman, int npoints, const double *x,
             double *y, double *slope, double *curve) {
    int i,j;
    double pj=0.0,pjPlus1=0.0,pjdash=0.0;
    double pjPlus1dash=0.0,pjdashdash=0.0,pjPlus1dashdash=0.0;
    double x1 = xData[0];
    double x2 = xData[1];
    double x3 = xData[2];
    double xnMin2 = xData[n-3];
    double xnMin1 = xData[n-2];
    double xn    = xData[n-1];
    double y1 = yData[0];
    double y2 = yData[1];
    double y3 = yData[2];
    double ynMin2 = yData[n-3];
    double ynMin1 = yData[n-2];
    double yn    = yData[n-1];
    double A1,B1,An,Bn,p2dash,pnMin1dash;        
    static char routine[] = "spliced";

    /* checking input data */
    if (n<=2)
    {
        throw ModelException(routine,
            "number of y points is lower than 2. Not supported for splice interpolation.");    
    }
    
    for (j=1 ; j < n ; j++)
    {
        if (xData[j] == xData[j-1])
        {
            throw ModelException(routine,
                "x points are duplicated");
        }
        if (xData[j] < xData[j-1])
        {
            throw ModelException(routine,
                "x points not increasing");
        
        }
    }

    for (i=0; i<npoints ; i++)
    {

        if      (x[i]  < x1)
        {
            if (wantY)             y[i] =    y1-gamma1*(x1-x[i]);
            if (wantSlope)    slope[i] =     gamma1;
            if (wantCurve)     curve[i] =     0.0;
        }
        else if (xn <= x[i])
        {
            if (wantY)          y[i] =    yn+gamman*(x[i]- xn);
            if (wantSlope)     slope[i] =    gamman;
            if (wantCurve)     curve[i] =    0.0;
        }
        else if ((x1 <= x[i]) && (x[i] < x2))
        {
            p2dash = (x2-x3)*y1/((x1-x2)*(x1-x3)) + y2/(x2-x1)
                    + y2/(x2-x3) + (x2-x1)*y3/ ((x3-x1)*(x3-x2));
            A1    =    ((x2-x1)*(p2dash+gamma1)-2*(y2-y1))/pow((x2-x1),3.0);
            B1    =    (3*(y2-y1)-(x2-x1)*
                (p2dash+2*gamma1))/((x2-x1)*(x2-x1));
            if (wantY)
                y[i]    =    (x[i]-x1)*(x[i]-x1)*(x[i]- x1)*A1+
                    (x[i]-x1)*(x[i]-x1)*B1+(x[i]-x1)*gamma1+y1;
            if (wantSlope)
                slope[i] = 3*(x[i]-x1)*(x[i]-x1)*A1+2*(x[i]-x1)*B1+gamma1;
            if (wantCurve)
                curve[i] = 6*(x[i]-x1)*A1+2*B1;
        
        }
        else if ((xnMin1 <= x[i]) && (x[i] < xn))
        {
            pnMin1dash = 
                (xnMin1-xn)*ynMin2/((xnMin2-xnMin1)*(xnMin2-xn))
                    + ynMin1/(xnMin1-xnMin2) + ynMin1/(xnMin1-xn) 
                    + (xnMin1-xnMin2)*yn /((xn-xnMin2)*(xn-xnMin1));
            An    =    ((xn-xnMin1)*(pnMin1dash+gamman)-2*(yn-ynMin1))/
                        pow((xn-xnMin1),3.0);
            Bn    =    ((xn-xnMin1)*(pnMin1dash+2*gamman)-3*(yn-ynMin1))/
                    ((xn-xnMin1)*(xn-xnMin1));
            if (wantY)
                y[i]    =    (x[i]-xn)*(x[i]-xn)*(x[i]-xn)*An+
                    (x[i]-xn)*(x[i]-xn)*Bn+(x[i]-xn)*gamman+yn;
            if (wantSlope)
                slope[i] = 3*(x[i]-xn)*(x[i]-xn)*An+2*(x[i]-xn)*Bn+gamman;
            if (wantCurve)
                curve[i] = 6*(x[i]-xn)*An+2*Bn;

        }
        else 
        {
            for (j=1; j<(n-2); j++)
            {
                if ((xData[j] <= x[i])    && (x[i] < xData[j+1])) 
                {
                    break;
                }
            }

            if ((wantY) || (wantSlope))
            {
                pj     =      ((x[i]- xData[j])*(x[i]- xData[j+1]))*yData[j-1]/
                        ((xData[j-1]-xData[j])*(xData[j-1]-xData[j+1]))    
                    +((x[i]- xData[j-1])*(x[i]- xData[j+1]))*yData[j]/
                        ((xData[j]-xData[j-1])*(xData[j]-xData[j+1]))
                    +((x[i]- xData[j-1])*(x[i]- xData[j]))*yData[j+1]/
                        ((xData[j+1]-xData[j-1])*(xData[j+1]-xData[j]));
                pjPlus1 = ((x[i]- xData[j+1])*(x[i]- xData[j+2]))*yData[j]/
                        ((xData[j]-xData[j+1])*(xData[j]-xData[j+2]))    
                    +((x[i]- xData[j])*(x[i]- xData[j+2]))*yData[j+1]/
                        ((xData[j+1]-xData[j])*(xData[j+1]-xData[j+2]))
                    +((x[i]- xData[j])*(x[i]- xData[j+1]))*yData[j+2]/
                        ((xData[j+2]-xData[j])*(xData[j+2]-xData[j+1]));
            }
            if ((wantSlope) || (wantCurve))
            {
                pjdash         =    
                    (2*x[i]-xData[j]-xData[j+1])*yData[j-1]/
                        ((xData[j-1]-xData[j])*(xData[j-1]-xData[j+1]))    
                    +(2*x[i]-xData[j-1]-xData[j+1])*yData[j]/
                        ((xData[j]-xData[j-1])*(xData[j]-xData[j+1]))
                    +(2*x[i]-xData[j-1]-xData[j])*yData[j+1]/
                        ((xData[j+1]-xData[j-1])*(xData[j+1]-xData[j]));
                pjPlus1dash  =
                        (2*x[i]- xData[j+1]- xData[j+2])*yData[j]/
                        ((xData[j]-xData[j+1])*(xData[j]-xData[j+2]))    
                    +(2*x[i]- xData[j]-xData[j+2])*yData[j+1]/
                        ((xData[j+1]-xData[j])*(xData[j+1]-xData[j+2]))
                    +(2*x[i]- xData[j]-xData[j+1])*yData[j+2]/
                        ((xData[j+2]-xData[j])*(xData[j+2]-xData[j+1]));
            }

            if (wantCurve)
            {
                pjdashdash         =
                        2*yData[j-1]/
                        ((xData[j-1]-xData[j])*(xData[j-1]-xData[j+1]))    
                    +2*yData[j]/
                        ((xData[j]-xData[j-1])*(xData[j]-xData[j+1]))
                    +2*yData[j+1]/
                        ((xData[j+1]-xData[j-1])*(xData[j+1]-xData[j]));
                pjPlus1dashdash    =
                        2*yData[j]/
                        ((xData[j]-xData[j+1])*(xData[j]-xData[j+2]))    
                    +2*yData[j+1]/
                        ((xData[j+1]-xData[j])*(xData[j+1]-xData[j+2]))
                    +2*yData[j+2]/
                        ((xData[j+2]-xData[j])*(xData[j+2]-xData[j+1]));

            }
            if (wantY)     
                    y[i] =      (xData[j+1]-x[i])*pj/(xData[j+1]-xData[j])
                        + (x[i]- xData[j])*pjPlus1/(xData[j+1]-xData[j]);
            if (wantSlope)
                slope[i] = -pj/(xData[j+1]-xData[j])
                        + pjPlus1/(xData[j+1]-xData[j])
                        + (xData[j+1]-x[i])*pjdash/(xData[j+1]-xData[j])
                        +(x[i]- xData[j])*pjPlus1dash/(xData[j+1]-xData[j]);
            if (wantCurve)
                curve[i] = -2*pjdash/(xData[j+1]-xData[j])
                        + 2*pjPlus1dash/(xData[j+1]-xData[j])
                        + (xData[j+1]-x[i])*pjdashdash/(xData[j+1]-xData[j])
                        + (x[i]- xData[j])*pjPlus1dashdash/
                        (xData[j+1]-xData[j]);
        }
    }
    
    return;

}


/**
 * This function interpolates in the strike dimension.
 * It uses splice interpolation.
 */
// HACK: uses xData[0] as a flag
// * if xData[0] < 0:
//   -> LINEAR interpolation before xData[2], using (0, yData[0]) and (xData[2], yData[2])
//   -> SPLICE interpolation after xData[2], using {(xData[i], yData[i])} points (i between 1 and n-1)
// * if xData[0] >= 0:
//   -> SPLICE interpolation using {(xData[i], yData[i])} points (i between 0 and n-1)
void SpliceStrikeInterp(
        int             n,      /* (I) nb of points in the input arrays */
        const double   *xData,  /* (I) x values                         */
        const double   *yData,  /* (I) y values                         */
        double          x,      /* (I) value of x that musst be interp   */
        double         *y)      /* (O) interpolated y value             */
{
    static const string routine("SpliceStrikeInterp");
    double gamma1 =  0.;
    double gamman =  0.;

    if (n == 0) {
      throw ModelException(routine,
        "number of points is 0 in strike data");
    } else if (n == 1) {
        gamma1 = gamman = 0.;
    } else {
        if (xData[0] < 0) {
            // relevant points start at (xdata[1], ydata[1])
            gamma1 = (yData[2] - yData[1]) / (xData[2] - xData[1]);
        } else {
            // relevant points start at (xdata[0], ydata[0])
            gamma1 = (yData[1] - yData[0]) / (xData[1] - xData[0]);
        }
        gamman = (yData[n-1] - yData[n-2]) / (xData[n-1] - xData[n-2]);
    }
    
    if (xData[0] < 0) {
        if (x <= xData[2]) {
            // LINEAR interpolation before xData[2], using (0, yData[0]) and (xData[2], yData[2])
            linterp(x, y, 0, xData[2], yData[0], yData[2]);
       } else {
            // SPLICE interpolation after xData[2], using {(xData[i], yData[i])} points (i between 1 and n-1)
            spliced(1, 0, 0, n-1, xData+1, yData+1, gamma1, gamman, 1, &x, y , NULL, NULL);
      }
    } else {
        // SPLICE interpolation using {(xData[i], yData[i])} points (i between 0 and n-1)
        spliced(1, 0, 0, n, xData, yData, gamma1, gamman, 1, &x, y , NULL, NULL);
    }
}


/**
 * This function interpolates in the strike dimension.
 * It uses spline interpolation.
 */
// HACK: uses xData[0] as a flag
// * if xData[0] < 0:
//   -> LINEAR interpolation before xData[2], using (0, yData[0]) and (xData[2], yData[2])
//   -> LINEAR interpolation after xData[n-1], using (xData[n-2], yData[n-2]) and (xData[n-1], yData[n-1])
//   -> SPLINE interpolation between [xData[2]-xData[n-1]], using {(xData[i], yData[i])} points (i between 1 and n-1)
// * if xData[0] >= 0:
//   -> LINEAR interpolation before xData[0], using (xData[0], yData[0]) and (xData[1], yData[1])
//   -> LINEAR interpolation after xData[n-1], using (xData[n-2], yData[n-2]) and (xData[n-1], yData[n-1])
//   -> SPLINE interpolation between [xData[0]-xData[n-1]], using {(xData[i], yData[i])} points (i between 0 and n-1)
void SplineStrikeInterp(
        int             n,      /* (I) nb of points in the input arrays */
        const double   *xData,  /* (I) x values                         */
        const double   *yData,  /* (I) y values                         */
        double          x,      /* (I) value of x that musst be interp   */
        double         *y)      /* (O) interpolated y value             */
{
    static const string routine("SplineStrikeInterp");

    if (n == 0) {
      throw ModelException(routine,
        "number of points is 0 in strike data");
    } 

    if (xData[0] < 0) {
        if (x <= xData[2]) {
            // LINEAR interpolation before xData[2], using (0, yData[0]) and (xData[2], yData[2])
            linterp(x, y, 0.0, xData[2], yData[0], yData[2]);
        } else if (x >= xData[n-1]) {
            // LINEAR interpolation after xData[n-1], using (xData[n-2], yData[n-2]) and (xData[n-1], yData[n-1])
            linterp(x, y, xData[n-2], xData[n-1], yData[n-2], yData[n-1]);
        } else {
            // builds NRSpline with left and right 2nd derivatives = 0
            NRSpline spline(2, 0, 2, 0);
            // SPLINE interpolation between [xData[2]-xData[n-1]], using {(xData[i], yData[i])} points (i between 1 and n-1)
            Interpolator::InterpolantConstSP interpolator = spline.computeInterp(xData + 1, yData + 1, n - 1);
            *y = interpolator->value(x);    
        }
    } else {
        if (x <= xData[0]) {
            // LINEAR interpolation before xData[0], using (xData[0], yData[0]) and (xData[1], yData[1])
            linterp(x, y, xData[0], xData[1], yData[0], yData[1]);
        } else if (x >= xData[n-1]) {
            // LINEAR interpolation after xData[n-1], using (xData[n-2], yData[n-2]) and (xData[n-1], yData[n-1])
            linterp(x, y, xData[n-2], xData[n-1], yData[n-2], yData[n-1]);
        } else {
            // builds NRSpline with left and right 2nd derivatives = 0
            NRSpline spline(2, 0, 2, 0);
            // SPLINE interpolation between [xData[0]-xData[n-1]], using {(xData[i], yData[i])} points (i between 0 and n-1)
            Interpolator::InterpolantConstSP interpolator = spline.computeInterp(xData, yData, n);
            *y = interpolator->value(x);    
        }
    }
}
/**
* This function interpolates in the strike dimension.
* It uses linear interpolation and flat extrapolation.
*/
void LinearStrikeInterp(
						int             n,      /* (I) nb of points in the input arrays */
						const double   *xData,  /* (I) x values (need to be increasing) */
						const double   *yData,  /* (I) y values                         */
						double          x,      /* (I) value of x that must be interp   */
						double         *y)      /* (O) interpolated y value             */
{
	static const string routine("LinearStrikeInterp");

	if(xData == NULL)
	{
		throw ModelException("xData is NULL", routine);
	}
	if(yData == NULL)
	{
		throw ModelException("yData is NULL", routine);
	}
	if (n == 0) {
		throw ModelException(routine,
			"number of points is 0 in strike data");
	} 
	// check if x is outside range in which case we use flat extrapolation
	if( x <= xData[0])
	{
		*y = yData[0];
	}
	else if (x >= xData[n-1])
	{
		*y = yData[n-1];
	}
	else 
	{
		// if here means that we do linear interpolation
		/** index in xData such that xData[highIdx] >= x */
		int highIdx = 0;
		while ( highIdx < n && xData[highIdx] < x ) highIdx++;
		
		if(highIdx < 1 || highIdx >= n)
		{
			// this can't happen
			throw ModelException("highIdx out of bounds", routine);
		}
		// call linterp routine to calculate y 
		linterp(x, y, xData[highIdx-1], xData[highIdx], yData[highIdx-1], yData[highIdx]);
	}


} /** LinearStrikeInterp */

/**
* This function interpolates in the strike dimension.
* It uses flat interpolation and flat extrapolation.
*/
void FlatStrikeInterp(
						int             n,      /* (I) nb of points in the input arrays */
						const double   *xData,  /* (I) x values (need to be increasing) */
						const double   *yData,  /* (I) y values                         */
						double          x,      /* (I) value of x that must be interp   */
						double         *y)      /* (O) interpolated y value             */
{
	static const string routine("FlatStrikeInterp");

	if(xData == NULL)
	{
		throw ModelException("xData is NULL", routine);
	}
	if(yData == NULL)
	{
		throw ModelException("yData is NULL", routine);
	} 

	if (x >= xData[n-1])
	{
		// if x is past range then assign last point
		*y = yData[n-1];
	}
	else 
	{	
		/** index in xData such that xData[highIdx] >= x */
		int highIdx = 0;
		// note that because of previous check, highIdx always < n
		while ( highIdx < n-1 && xData[highIdx] < x ) highIdx++;
		

		// assign y value for xData point just after x
		*y = yData[highIdx];
	}


} /** FlatStrikeInterp */




DRLIB_END_NAMESPACE
