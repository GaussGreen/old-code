//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolHyperTrig.cpp
//
//   Description : Parameterized Volatility
//                 The parameterization is based on a specification of the local
//                 volatilities in terms of hyperbolic trigonometric functions,
//                 and is a variation to the local volatility parameterization
//                 developped by Henrik Rasmussen for long-term FX volatilities.
//                 The local volatility parameterization is made of two components:
//                      - the skew (in the form of 'the tanh function')
//                      - the convexity (in the form of '1.0 minus the sech function')
//                 This implementation uses an approximation (Gatheral's formula)
//                 to derive the corresponding implied volatilities. This means the
//                 functional is close to arbitrage-free both in the time dimension 
//                 and in the strike dimension.
//
//
//   Author      : Regis Guichard
//
//   Date        : 15 Oct 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolBaseParam.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Addin.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/VHTSkewFactorsPointwise.hpp"
#include "edginc/VHTSkewCelerityFactorsPointwise.hpp"
#include "edginc/VHTConvexityFactorsPointwise.hpp"
#include "edginc/VHTConvexityCelerityFactorsPointwise.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/VolAtmParallelConstConvx.hpp"
#include "edginc/VolAtmPointwiseConstConvx.hpp"
// for additional addin
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VarSwapBasis.hpp"
#include "edginc/DeltaToStrike.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/Delta.hpp"
#include "edginc/Results.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/BenchmarkDate.hpp"

// #define CREATE_MODIFIED_XML_FILES

DRLIB_BEGIN_NAMESPACE

/** Pade approximation for log(1+exp(-x)) */
class PadeApproxLogCosh: public CObject {
public:
    static CClassConstSP const TYPE;

    PadeApproxLogCosh(int nbPadeCoeffs): CObject(TYPE), nbPadeCoeffs(nbPadeCoeffs) {
        coeffs = getLogCoshPadeCoeffs(nbPadeCoeffs);
    }

    PadeApproxLogCosh(): CObject(TYPE), nbPadeCoeffs(0) { }

    ~PadeApproxLogCosh() {}

    /** log(ch(x)) = x - ln(2) + ln(1+exp(-2x))
    The Pade approximant is used on ln(1+exp(-x)) */
    double operator() (double x) const {
        double y = fabs(x);
        double result = y - logTwo + ratval((2.0 * y), &coeffs[0], nbPadeCoeffs, nbPadeCoeffs);
        return result;
    }

private:
    static void load(CClassSP& clazz){
        REGISTER(PadeApproxLogCosh, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(nbPadeCoeffs, "Number of Pade coefficients");
        FIELD(coeffs,       "Pade coeffiecients");
    }

    static IObject* defaultCtor(){
        return new PadeApproxLogCosh();
    }

    /** Pade coefficients of ln(1+exp(-x)) */
    static DoubleArray getLogCoshPadeCoeffs(int nbPadeCoeffs) {
        if(nbPadeCoeffs > 10) {
            throw ModelException("getLogCoshPadeCoeffs",
                                 "The number of Pade coefficients has to be <= 10");
        }

        if(nbPadeCoeffs < 0) {
            throw ModelException("getLogCoshPadeCoeffs",
                                 "The number of Pade coefficients has to be > 0");
        }

        DoubleArray coeffs(21);
        // Taylor coefficients
        coeffs[0]  = logTwo;
        coeffs[1]  = -1.0/2.0;
        coeffs[2]  = +1.0/8.0;
        coeffs[3]  =  0.0;
        coeffs[4]  = -1.0/192.0;
        coeffs[5]  =  0.0;
        coeffs[6]  = +1.0/2880.0;
        coeffs[7]  =  0.0;
        coeffs[8]  = -17.0/645120.0;
        coeffs[9]  =  0.0;
        coeffs[10] = +31.0/14515200.0;
        coeffs[11] =  0.0;
        coeffs[12] = -691.0/3832012800.0;
        coeffs[13] =  0.0;
        coeffs[14] = +5461.0/348713164800.0;
        coeffs[15] =  0.0;
        coeffs[16] = -929569.0/669529276416000.0;
        coeffs[17] =  0.0;
        coeffs[18] = +3202291.0/25609494822912000.0;
        coeffs[19] =  0.0;
        coeffs[20] = -221930581.0/19463216065413120000.0;

        for(int i = 20; i > 2 * nbPadeCoeffs; i--) {
            // Drop the last coefficients to be on the safe side
            coeffs.pop_back();
        }

        double resid = 0.0;
        // calculate Pade coefficients
        pade(&coeffs[0], nbPadeCoeffs, &resid);

        return coeffs;
    }

    static const double logTwo;
    int                 nbPadeCoeffs;
    mutable DoubleArray coeffs;         // mutable because ratval expects non-const double*
};

const double PadeApproxLogCosh::logTwo = log(2.0);

CClassConstSP const PadeApproxLogCosh::TYPE =
CClass::registerClassLoadMethod("PadeApproxLogCosh", typeid(PadeApproxLogCosh), PadeApproxLogCosh::load);

/** Pade approximation for -2 * (arctan(exp(x) - Pi/4) */
class PadeApproxAtanExp: public CObject {
public:
    static CClassConstSP const TYPE;

    PadeApproxAtanExp(): CObject(TYPE), nbPadeCoeffs(0) { }

    PadeApproxAtanExp(int nbPadeCoeffs): CObject(TYPE), nbPadeCoeffs(nbPadeCoeffs) {
        coeffs = getAtanExpPadeCoeffs(nbPadeCoeffs);
    }

    ~PadeApproxAtanExp() {}

    /** x - 2 * (arctan(exp(x)) - Pi/4) = (x - Pi/2) + (-2 * (arctan(th(x/2)) - Pi/4))
        The Pade approximant is used on (-2 * (arctan(th(x)) - Pi/4)) */
    double operator() (double x) const {
        double y = fabs(x);
        double result = y - halfPi + ratval((0.5 * y), &coeffs[0], nbPadeCoeffs, nbPadeCoeffs);
        return (x >= 0) ? result : -result;
    }

private:
    static void load(CClassSP& clazz){
        REGISTER(PadeApproxAtanExp, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(nbPadeCoeffs, "Number of Pade coefficients");
        FIELD(coeffs,       "Pade coeffiecients");
    }

    static IObject* defaultCtor(){
        return new PadeApproxAtanExp();
    }

    /** Pade coefficients of -2 * (arctan(tanh(x)) - Pi/4) */
    static DoubleArray getAtanExpPadeCoeffs(int nbPadeCoeffs) {
        if(nbPadeCoeffs > 10) {
            throw ModelException("getAtanExpPadeCoeffs",
                                 "The number of Pade coefficients has to be <= 10");
        }

        if(nbPadeCoeffs < 0) {
            throw ModelException("getAtanExpPadeCoeffs",
                                 "The number of Pade coefficients has to be > 0");
        }

        DoubleArray coeffs(21);
        // Taylor coefficients
        coeffs[0]  = halfPi;
        coeffs[1]  = -2.0;
        coeffs[2]  =  0.0;
        coeffs[3]  = +4.0/3.0;
        coeffs[4]  =  0.0;
        coeffs[5]  = -4.0/3.0;
        coeffs[6]  =  0.0;
        coeffs[7]  = +488.0/315.0;
        coeffs[8]  =  0.0;
        coeffs[9]  = -1108.0/567.0;
        coeffs[10] =  0.0;
        coeffs[11] = +404168.0/155925.0;
        coeffs[12] =  0.0;
        coeffs[13] = -332648.0/93555.0;
        coeffs[14] =  0.0;
        coeffs[15] = +3189775696.0/638512875.0;
        coeffs[16] =  0.0;
        coeffs[17] = -912541748.0/127702575.0;
        coeffs[18] =  0.0;
        coeffs[19] = +19239037403528.0/1856156927625.0;
        coeffs[20] =  0.0;

        for(int i = 20; i > 2 * nbPadeCoeffs; i--) {
            // Drop the last coefficients to be on the safe side
            coeffs.pop_back();
        }

        double resid = 0.0;
        // calculate Pade coefficients
        pade(&coeffs[0], nbPadeCoeffs, &resid);

        return coeffs;
    }

    static const double halfPi;
    int                 nbPadeCoeffs;
    mutable DoubleArray coeffs;         // mutable because ratval expects non-const double*
};

const double PadeApproxAtanExp::halfPi = 0.5 * Maths::PI;

CClassConstSP const PadeApproxAtanExp::TYPE =
CClass::registerClassLoadMethod("PadeApproxAtanExp", typeid(PadeApproxAtanExp), PadeApproxAtanExp::load);

/** Approximation for log(cosh(x1)/cosh(x2)) for large x1, x2 */
static double truncatedLogCoshCosh(double x1, double x2) {
    double result = 0.0;
    // note: cosh(x) = maxDouble = 1.7E308 <=> x = 308 *ln(10) + ln(3.4) = 710.41
    if ((fabs(x1) > 700.0) || (fabs(x2) > 700.0)) {
        result = x1 - x2;
    } else {
        result = log (cosh(x1) / cosh(x2));
    }
    return result;
}

class VolHyperTrig;
class ParameterOutputAddin;

class VolHyperTrigPrimer: public CObject {
public:
    // VolGrid class to precompute and cache implied vols at benchmark dates
    class VolGrid: public CObject{
    public:
        static CClassConstSP const TYPE;

        static const string GRID_TYPE_NONE;
        static const string GRID_TYPE_LINEAR;
        static const string GRID_TYPE_SPLINE;
        static const string SPEEDUP_METHOD_GRID;

        typedef enum{GT_NONE, GT_LINEAR, GT_SPLINE} GRID_Type;

        // constructor
        VolGrid():
            CObject(TYPE),
            gridType(GRID_TYPE_NONE),
            gridEnumType(GT_NONE),
            numGrid(1000),
            TruncationStd(5.0){}

        VolGrid(string  gridType,
                int     numGrid,
                double  TruncationStd):
            CObject(TYPE),
            gridType(gridType),
            numGrid(numGrid),
            TruncationStd(TruncationStd){
                if( CString::equalsIgnoreCase(gridType, GRID_TYPE_LINEAR) )
                    gridEnumType = GT_LINEAR;
                else if( CString::equalsIgnoreCase(gridType, GRID_TYPE_SPLINE) ) 
                    gridEnumType = GT_SPLINE;
                else 
                    gridEnumType = GT_NONE;            
            }

        void validatePop2Object(){
            static const string method("VolHyperTrigPrimer::VolGrid::ValidatePop2Object");
            try{
                if( !CString::equalsIgnoreCase(gridType, GRID_TYPE_NONE) )
                {
                    if( !CString::equalsIgnoreCase(gridType, GRID_TYPE_LINEAR) &&
                        !CString::equalsIgnoreCase(gridType, GRID_TYPE_SPLINE) )
                        throw ModelException(method, "Only N, L, S allowed for grid type");
                    if( !Maths::isPositive(TruncationStd) )
                        throw ModelException(method, "TruncationStd must be positive");
                    // make sure numGrid is odd to have a grid point at middle
                    if( numGrid < 2 )
                        throw ModelException(method, "num grid at least 2");
                    numGrid = (numGrid/2)*2 + 1;
                }
            } catch(exception& e) {
                throw ModelException(e,method);
            }
        }

        string getGridType(){
            return gridType;
        }

        void setupGrid(VolHyperTrigPrimer* vhtPrimer){
            const double y1 = 2e30, yn = 2e30; // default to natural spline
            int iGrid;

            // create cache of vols for benchmark dates
            xGrid.resize(numGrid);
            for(iGrid=0; iGrid<numGrid; iGrid++)
                xGrid[iGrid] = -TruncationStd + 2.0 * iGrid * TruncationStd / (numGrid-1);

            int numTradYears = vhtPrimer->tradYears.size();
            varGrid.resize(numTradYears);
            varGrid2.resize(numTradYears);

            double lastAtmVar = 0;
            double totalAtmVar = 0;
            int iTradYear = 0;
            for (; iTradYear < numTradYears; ++iTradYear){
                varGrid[iTradYear].resize(numGrid);
                varGrid2[iTradYear].resize(numGrid);

                totalAtmVar = vhtPrimer->atmVars[iTradYear];
                double A = vhtPrimer->backboneSkewParams[iTradYear];
                double B = vhtPrimer->backboneConvexityParams[iTradYear];
                double C = vhtPrimer->backboneSkewCelerities[iTradYear];
                double D = vhtPrimer->backboneConvexityCelerities[iTradYear];
                double skewComp, convComp;
                for(iGrid=0; iGrid<numGrid; iGrid++) {
                    // special case for ATM
                    if(Maths::isZero(xGrid[iGrid])) {
                        varGrid[iTradYear][iGrid] = totalAtmVar;
                        continue;
                    }

                    double& var = varGrid[iTradYear][iGrid];
					double thisx = xGrid[iGrid];
                    double lastx = thisx * lastAtmVar / totalAtmVar;
                    double deltax = thisx - lastx;

                    // interpolate vol at current grid for next step
                    if(iTradYear==0)
                        var = 0.0;
                    else
						splint(&*xGrid.begin()-1, &*varGrid[iTradYear-1].begin()-1, &*varGrid2[iTradYear-1].begin()-1,
                                 xGrid.size(), lastx, &var);

                    /**********************************************************************************
                        * Some background on the equivalence of the two calls:
                        * d = x1 - x2
                        * 1 - 2 * (atan(exp(x1) - atan(exp(x2)) / d =
                        * ( x1 - 2 * atan(exp(x1) + pi/2.0) - (x2 - 2 * atan(exp(x2)) + pi/2.0) ) / d =
                        * ( fastAtanExp(x1) - fastAtanExp(x2) ) / d
                    **********************************************************************************/
                    if(vhtPrimer->usePade) {
                        // use the Pade approximation
                        skewComp = A * (vhtPrimer->fastLogCosh(C * thisx) - vhtPrimer->fastLogCosh(C * lastx))
                            / (C * deltax);
                        // fudge a bit
                        convComp = Maths::isZero(D) ? 0.0
                            : B * (vhtPrimer->fastAtanExp(D * thisx) - vhtPrimer->fastAtanExp(D * lastx))
                            / (D * deltax);
                    } else {
                        // use the closed-formula
                        skewComp = A * truncatedLogCoshCosh(C * thisx, C * lastx) / (C * deltax);
                        //skewComp = A * log(cosh(C * thisx) / cosh(C * lastx)) / (C * deltax);
                        // fudge a bit
                        convComp = Maths::isZero(D) ? 0.0
                            : B * (1.0 - 2.0 * (atan(exp(D * thisx)) - atan(exp(D * lastx)))
                            / (D * deltax));
                    }

                    double normalizedVarIncrement = 1.0 + skewComp + convComp;
                    // for extreme values of parameters the 'squared' local vol becomes
                    // negative !!! Therefore, we need to fudge things here...
                    normalizedVarIncrement = Maths::max(0.0, normalizedVarIncrement);
                    var += (totalAtmVar - lastAtmVar) * normalizedVarIncrement;
                }

                spline(&*xGrid.begin()-1, &*varGrid[iTradYear].begin()-1, numGrid, y1, yn, &*varGrid2[iTradYear].begin()-1);

                lastAtmVar = totalAtmVar;
            }

            if( CString::equalsIgnoreCase(gridType, GRID_TYPE_LINEAR) )
            {
                varGrid2.clear();
            }
        }

        // if outside grid, do flat extrapolation. otherwise, linear/spline interp
        double interpImpVar(double      xnorm,
                            int         gridIdx){
            static const string method("VolHyperTrigPrimer::VolGrid::interpImpVar");

            double var = 0.0;
            if( xnorm <= xGrid.front() )
                var = varGrid[gridIdx].front();
            else if( xnorm >= xGrid.back() )
                var = varGrid[gridIdx].back();
            else if( gridEnumType == GT_LINEAR ) {

                double lo = (xnorm - xGrid[0])/(xGrid[1]-xGrid[0]);
                int klo = (int)lo;
                double *vlo = &varGrid[gridIdx][klo];
                var = *vlo + (lo-klo) * (*(vlo+1) - *(vlo));

            } else if( gridEnumType == GT_SPLINE ) {

                splint(&*xGrid.begin()-1, &*varGrid[gridIdx].begin()-1, &*varGrid2[gridIdx].begin()-1,
                    xGrid.size(), xnorm, &var);

            } else {
                throw ModelException(method, "no interpolation of implied var if gridType is "+ GRID_TYPE_NONE);
            }

            return var;
        }

        ~VolGrid(){};

    protected:
        string      gridType;
        int         numGrid;
        double      TruncationStd;

        //transient fields
        mutable DoubleArray         xGrid;
        mutable DoubleArrayArray    varGrid;
        mutable DoubleArrayArray    varGrid2;
        int     gridEnumType;           // to use internally to reduce calc time.

    private:
        static void load(CClassSP& clazz){
            REGISTER(VolGrid, clazz);
            SUPERCLASS(CObject);
            EMPTY_SHELL_METHOD(defaultVolGrid);
            FIELD(gridType, "type of vol interpolation: NONE, LINEAR or SPLINE");
            FIELD_MAKE_OPTIONAL(gridType);
            FIELD(numGrid, "number of grid points");
            FIELD_MAKE_OPTIONAL(numGrid);
            FIELD(TruncationStd, "trancation of standard normal variable");
            FIELD_MAKE_OPTIONAL(TruncationStd);
            FIELD(xGrid, "");
            FIELD_MAKE_TRANSIENT(xGrid);
            FIELD(varGrid, "");
            FIELD_MAKE_TRANSIENT(varGrid);
            FIELD(varGrid2, "");
            FIELD_MAKE_TRANSIENT(varGrid2);
            FIELD(gridEnumType, "");
            FIELD_MAKE_TRANSIENT(gridEnumType);
        }

        static IObject* defaultVolGrid(){
            return new VolGrid();
        }
    };
    typedef smartPtr<VolGrid> VolGridSP;


    static CClassConstSP const TYPE;
    friend class ParameterOutputAddin;
    friend class VolHyperTrigPrimer::VolGrid;

    VolHyperTrigPrimer(): CObject(TYPE) {}
    VolHyperTrigPrimer(CClassConstSP clazz): CObject(clazz) {}

    ~VolHyperTrigPrimer() {}

    // Reduced constructor with 2 term structures
    VolHyperTrigPrimer(CClassConstSP        clazz,
                       const string&        approxPade,
                       const string&        tenorRef,
                       const StringArray&   strikeRefRatioBMs,
                       const DoubleArray&   strikeRefRatios,
                       double               strikeRef,
                       const StringArray&   skewFactorBMs,
                       const DoubleArray&   skewFactors,
                       double               skew,
                       const StringArray&   convexityFactorBMs,
                       const DoubleArray&   convexityFactors,
                       double               convexity,
                       const string&        tenorTime,
                       const bool           useTimeDependentCelerities,
                       double               convexityPower,
                       double               skewBackboneScaling,
                       double               convexityBackboneScaling,
                       VolSurfaceConstSP    backbone,
                       VarSwapBasisSP       basis,
                       DateTime             markedDate):
    CObject(clazz), approxPade(approxPade),
    tenorRef(tenorRef), strikeRefRatioBMs(strikeRefRatioBMs), strikeRefRatios(strikeRefRatios), strikeRef(strikeRef),
    skewFactorBMs(skewFactorBMs), skewFactors(skewFactors), skew(skew),
    convexityFactorBMs(convexityFactorBMs), convexityFactors(convexityFactors), convexity(convexity),
    skewCelerity(0.0), convexityCelerity(0.0),
    tenorTime(tenorTime), useTimeDependentCelerities(useTimeDependentCelerities), convexityPower(convexityPower),
    skewBackboneScaling(skewBackboneScaling), convexityBackboneScaling(convexityBackboneScaling), backbone(backbone),
    basis(basis), markedDate(markedDate) {

        static const string method("VolHyperTrigPrimer::VolHyperTrigPrimer");
        validatePop2Object();
    }

    // Reduced constructor with 2 term structures + VolGrid
    VolHyperTrigPrimer(CClassConstSP        clazz,
                       const string&        approxPade,
                       const string&        tenorRef,
                       const StringArray&   strikeRefRatioBMs,
                       const DoubleArray&   strikeRefRatios,
                       double               strikeRef,
                       const StringArray&   skewFactorBMs,
                       const DoubleArray&   skewFactors,
                       double               skew,
                       const StringArray&   convexityFactorBMs,
                       const DoubleArray&   convexityFactors,
                       double               convexity,
                       const string&        tenorTime,
                       const bool           useTimeDependentCelerities,
                       double               convexityPower,
                       double               skewBackboneScaling,
                       double               convexityBackboneScaling,
                       VolSurfaceConstSP    backbone,
                       string               gridType,
                       int                  numGrid,
                       double               TruncationStd,
                       bool                 useGrid,
                       VarSwapBasisSP       basis,
                       DateTime             markedDate):
    CObject(clazz), approxPade(approxPade),
    tenorRef(tenorRef), strikeRefRatioBMs(strikeRefRatioBMs), strikeRefRatios(strikeRefRatios), strikeRef(strikeRef),
    skewFactorBMs(skewFactorBMs), skewFactors(skewFactors), skew(skew),
    convexityFactorBMs(convexityFactorBMs), convexityFactors(convexityFactors), convexity(convexity),
    skewCelerity(0.0), convexityCelerity(0.0),
    tenorTime(tenorTime), useTimeDependentCelerities(useTimeDependentCelerities), convexityPower(convexityPower),
    skewBackboneScaling(skewBackboneScaling), convexityBackboneScaling(convexityBackboneScaling), backbone(backbone),
    useGrid(useGrid),volGrid(VolGridSP(new VolGrid(gridType,numGrid,TruncationStd))),basis(basis), markedDate(markedDate){

        static const string method("VolHyperTrigPrimer::VolHyperTrigPrimer");
        // Dummy try/catch to avoid vc6.opt crash
		try {
            validatePop2Object();
        } catch (...) { throw; }
    }


    // Full constructor with 4 term structures
    VolHyperTrigPrimer(CClassConstSP        clazz,
                       const string&        approxPade,
                       const string&        tenorRef,
                       const StringArray&   strikeRefRatioBMs,
                       const DoubleArray&   strikeRefRatios,
                       double               strikeRef,
                       const StringArray&   skewFactorBMs,
                       const DoubleArray&   skewFactors,
                       double               skew,
                       const StringArray&   convexityFactorBMs,
                       const DoubleArray&   convexityFactors,
                       double               convexity,
                       const StringArray&   skewCelerityFactorBMs,
                       const DoubleArray&   skewCelerityFactors,
                       double               skewCelerity,
                       const StringArray&   convexityCelerityFactorBMs,
                       const DoubleArray&   convexityCelerityFactors,
                       double               convexityCelerity,
                       const string&        tenorTime,
                       const bool           useTimeDependentCelerities,
                       double               convexityPower,
                       double               skewBackboneScaling,
                       double               convexityBackboneScaling,
                       VolSurfaceConstSP    backbone,
                       string               gridType,
                       int                  numGrid,
                       double               TruncationStd,
                       bool                 useGrid,
                       VarSwapBasisSP       basis,
                       DateTime             markedDate):
    CObject(clazz), approxPade(approxPade),
    tenorRef(tenorRef), strikeRefRatioBMs(strikeRefRatioBMs), strikeRefRatios(strikeRefRatios), strikeRef(strikeRef),
    skewFactorBMs(skewFactorBMs), skewFactors(skewFactors), skew(skew),
    convexityFactorBMs(convexityFactorBMs), convexityFactors(convexityFactors), convexity(convexity),
    skewCelerityFactorBMs(skewCelerityFactorBMs), skewCelerityFactors(skewCelerityFactors), skewCelerity(skewCelerity),
    convexityCelerityFactorBMs(convexityCelerityFactorBMs), convexityCelerityFactors(convexityCelerityFactors), convexityCelerity(convexityCelerity),
    tenorTime(tenorTime), useTimeDependentCelerities(useTimeDependentCelerities), convexityPower(convexityPower),
    skewBackboneScaling(skewBackboneScaling), convexityBackboneScaling(convexityBackboneScaling), backbone(backbone),
    useGrid(useGrid), volGrid(VolGridSP(new VolGrid(gridType,numGrid,TruncationStd))),basis(basis),markedDate(markedDate){

        static const string method("VolHyperTrigPrimer::VolHyperTrigPrimer");
        validatePop2Object();
    }

    VarSwapBasis::VarSwapBasisProc* setup(const CAsset* asset,
                                          const VolSurface* backbone,
                                          const VarSwapBasis::VarSwapBasisRequest* req,
                                          const string& tenorTime) {
        return basis->setup(asset, backbone, req, tenorTime);
    }

    // the update method is helper class specific and thus must be implemented there
    virtual void update(double             skewMod,
                        double             convexityMod,
                        DoubleArray        skewFactorsMod,
                        DoubleArray        convexityFactorsMod,
                        double             skewCelerity,
                        double             convexityCelerity) {
        static const string method = "VolHyperTrigPrimer::update";
        throw ModelException(method, "Not implemented");
    }

    // the update method is helper class specific and thus must be implemented there
    virtual void update(double skew,
                        double convexity,
                        double skewCelerity,
                        double convexityCelerity,
                        const DoubleArray& skewFactors,
                        const DoubleArray& convexityFactors,
                        const DoubleArray& skewCelerityFactors,
                        const DoubleArray& convexityCelerityFactors) {
        static const string method = "VolHyperTrigPrimer::update";
        throw ModelException(method, "Not implemented");
    }


private:
    void calcFactorTerms(const StringArray&  BMs,
                         const DoubleArray&  factors,
                         DateTimeArray&      dates,
                         DoubleArray&        tradYears,
                         string              what){
        static const string method = "VolHyperTrigPrimer::calcFactorTerms";
        int date1YFound = false;
        int iBM = 0;
        for (; iBM < BMs.size(); ++iBM){
            MaturityPeriod period(BMs[iBM]);
            dates[iBM] = DateTime(period.toDate(baseDate).getDate(), tenorTimeInt);
            if (iBM > 0
                && dates[iBM - 1] >= dates[iBM]){
                throw ModelException(method,
                                     "The "
                                     + what
                                     + " benchmarks must be increasing. Compare the "
                                     + Format::toString(iBM)
                                     + "-th benchmark date ("
                                     + dates[iBM - 1].toString()
                                     + ") with the "
                                     + Format::toString(iBM + 1)
                                     + "-th benchmark date ("
                                     + dates[iBM].toString()
                                     + ")");
            }

            // Dummy try/catch to ""solve"" vc6.opt crashing

            try {
                if (dates[iBM] == dateRef){
                    if (!Maths::isZero(factors[iBM] - 1.0)){
                        throw ModelException(method,
                                             "the "
                                             + Format::toString(iBM + 1)
                                             + "-th "
                                             + what
                                             + " ("
                                             + Format::toString(factors[iBM])
                                             + ") should be equal to 1.0,\nas the "
                                             + Format::toString(iBM + 1)
                                             + "-th "
                                             + what
                                             + " benchmark date is equal to the reference date ("
                                             + dates[iBM].toString()
                                             + ")");
                    }
                    date1YFound = true;
                }
            }
            catch (...) {
                throw;
            }

            tradYears[iBM]
                = timeMetric->yearFrac(baseDate, dates[iBM]);
        }
        // if the 1y ref date hasn't been found, need to insert it
        if (!date1YFound){
            DateTimeArray datesNew(BMs.size() + 1);
            DoubleArray tradYearsNew(BMs.size() + 1);
            iBM = 0;
            int iNewBM = 0;
            while (iBM < BMs.size()
                   && dates[iBM] < dateRef){
                datesNew[iNewBM] = dates[iBM];
                tradYearsNew[iNewBM] = tradYears[iBM];
                ++iBM;
                ++iNewBM;
            }
            datesNew[iNewBM] = dateRef;
            tradYearsNew[iNewBM] = tradYearRef;
            ++iNewBM;
            while (iBM < BMs.size()){
                datesNew[iNewBM] = dates[iBM];
                tradYearsNew[iNewBM] = tradYears[iBM];
                ++iBM;
                ++iNewBM;
            }
            dates = datesNew;
            tradYears = tradYearsNew;
        }
    }

    double calcVol(double   strike,
                   double   tradYear) const {
        static const string method("VolHyperTrigPrimer:calcVol");
        try{
			// Determine tradYearsIdx such that: tradYears[tradYearsIdx-1] <= tradYear < tradYears[tradYearsIdx]
			int tradYearsIdx = 0;
			while (tradYearsIdx < tradYears.size() && tradYears[tradYearsIdx] <= tradYear){
                tradYearsIdx++;
            }
            
			//skew cutoff is relevant here, not price cutoff. Get skew cutoff strike
            if(basis.get()) {
                double cutoffStrike = basis->interpSkewCutoff(tradYear);
                strike = Maths::max(strike, cutoffStrike);
            }

            if (Maths::isZero(tradYear)){
                return sqrt(atmVars[0] / tradYears[0]);
            }
            double k = strike / getStrikeRef(tradYear);

            double totalAtmVar;
            if (tradYearsIdx <= 0){
                totalAtmVar = atmVars[0] / tradYears[0] * tradYear;
            }
            else if(tradYearsIdx > tradYears.size() - 1){
                totalAtmVar = atmVars.back() / tradYears.back() * tradYear;
            }
            else{
                totalAtmVar = atmVars[tradYearsIdx - 1]
                              + (atmVars[tradYearsIdx] - atmVars[tradYearsIdx - 1])
                                / (tradYears[tradYearsIdx] - tradYears[tradYearsIdx - 1])
                                * (tradYear - tradYears[tradYearsIdx - 1]);
            }

            if(Maths::isZero(k)) {
                // Work out the limit for x = -inf

                // Integrate var(0) * (1 - A(t) + B(t)) in [0, T]
                double var = 0.0;
                double lastAtmVar = 0.0;
                // Integrate over benchmarks
                int iTradYear;
                for (iTradYear = 0; iTradYear < tradYearsIdx && iTradYear < tradYears.size(); ++iTradYear){
                    double At = backboneSkewParams[iTradYear];
                    double Bt = backboneConvexityParams[iTradYear];
                    double varIncrement = Maths::max(1.0 - At + Bt, 0.0);
                    double thisAtmVar = atmVars[iTradYear];
                    double atmVarDiff = thisAtmVar - lastAtmVar;
                    var += atmVarDiff * varIncrement;
                    lastAtmVar = thisAtmVar;
                }

                // Add contribution of last point
                if(totalAtmVar > lastAtmVar) {
                    iTradYear = Maths::min(Maths::max(tradYearsIdx, 0), tradYears.size() - 1);
                    double At = backboneSkewParams[iTradYear];
                    double Bt = backboneConvexityParams[iTradYear];
                    double varIncrement = Maths::max(1.0 - At + Bt, 0.0);
                    double thisAtmVar = atmVars[iTradYear];
                    double atmVarDiff = thisAtmVar - lastAtmVar;
                    var += atmVarDiff * varIncrement;
                }

                // Return vol
                return sqrt(var / tradYear);
            }

            double x = log(k);
            double var = 0.0;
            double lastx = 0.0;
            double lastAtmVar = 0.0;
            int iTradYear = 0;

			double very_small_number = 0.00000001;
            if (tradYearsIdx > 0){
                if(!useGrid){

                    for (; iTradYear < tradYearsIdx && iTradYear < tradYears.size(); ++iTradYear){

                        double normalizedVarIncrement = 1.0;
                        double thisAtmVar = atmVars[iTradYear];
                        double thisx = thisAtmVar / totalAtmVar * x;
                        double deltax = thisx - lastx;

                        if (!Maths::isZero(deltax)){
                            double A = backboneSkewParams[iTradYear];
                            double B = backboneConvexityParams[iTradYear];
                            double C = backboneSkewCelerities[iTradYear];
                            double D = backboneConvexityCelerities[iTradYear];
                            double skewComp, convComp;

                            /**********************************************************************************
                            * Some background on the equivalence of the two calls:
                            * d = x1 - x2
                            * 1 - 2 * (atan(exp(x1) - atan(exp(x2)) / d =
                            * ( x1 - 2 * atan(exp(x1) + pi/2.0) - (x2 - 2 * atan(exp(x2)) + pi/2.0) ) / d =
                            * ( fastAtanExp(x1) - fastAtanExp(x2) ) / d
                            **********************************************************************************/

							if( abs(x) < very_small_number ){
								// use taylor expansion ...
								skewComp = A * C * (lastx + thisx)/2.;
								convComp = B * D * D * (thisx*thisx + thisx*lastx + lastx*lastx)/6.;
							} else if (usePade) {
                                // use the Pade approximation
                                skewComp = A * (fastLogCosh(C * thisx) - fastLogCosh(C * lastx))
                                    / (C * deltax);
                                // fudge a bit
                                convComp = Maths::isZero(D) ? 0.0
                                    : B * (fastAtanExp(D * thisx) - fastAtanExp(D * lastx))
                                        / (D * deltax);
                            } else {
                                // use the closed-formula
                                skewComp = A * truncatedLogCoshCosh(C * thisx, C * lastx) / (C * deltax);
                                //skewComp = A * log(cosh(C * thisx) / cosh(C * lastx)) / (C * deltax);
                                // fudge a bit
                                convComp = Maths::isZero(D) ? 0.0
                                    : B * (1.0 - 2.0 * (atan(exp(D * thisx)) - atan(exp(D * lastx)))
                                        / (D * deltax));
                            }

                            normalizedVarIncrement += skewComp + convComp;
                            // for extreme values of parameters the 'squared' local vol becomes
                            // negative !!! Therefore, we need to fudge things here...
                            normalizedVarIncrement = Maths::max(0.0, normalizedVarIncrement);
						}

                        lastx = thisx;
                        double atmVarDiff = thisAtmVar - lastAtmVar;
                        var += atmVarDiff * normalizedVarIncrement;
                        lastAtmVar = thisAtmVar;
                    }
                } else {
                    int gridIdx = Maths::min(tradYearsIdx, tradYears.size())-1;
                    lastAtmVar = atmVars[gridIdx];
					double xnorm = x*lastAtmVar/totalAtmVar;
                    var = volGrid.get()->interpImpVar(xnorm,gridIdx);
                }
            }

            iTradYear = Maths::min(Maths::max(tradYearsIdx, 0), tradYears.size() - 1);
            lastx = x*lastAtmVar/totalAtmVar;
            double normalizedVarIncrement = 1.0;
            double deltax = x - lastx;
            if (!Maths::isZero(deltax)){
                double A = backboneSkewParams[iTradYear];
                double B = backboneConvexityParams[iTradYear];
                double C = backboneSkewCelerities[iTradYear];
                double D = backboneConvexityCelerities[iTradYear];
                double skewComp, convComp;

				if( abs(x) < very_small_number ){
					// use taylor expansion ...
					skewComp = A * C * (lastx + x)/2.;
					convComp = B * D * D * (x*x + x*lastx + lastx*lastx)/6.;
				} else if(usePade) {
                    // use the Pade approximation
                    skewComp = A * (fastLogCosh(C * x) - fastLogCosh(C * lastx))
                        / (C * deltax);
                    convComp = B * (fastAtanExp(D * x) - fastAtanExp(D * lastx))
                        / (D * deltax);
                } else {
                    // use the closed-formula
                    skewComp = A * truncatedLogCoshCosh(C * x, C * lastx) / (C * deltax);
                    //skewComp = A * log(cosh(C * x) / cosh(C * lastx)) / (C * deltax);
                    convComp = B * (1.0 - 2.0 * (atan(exp(D * x)) - atan(exp(D * lastx)))
                        / (D * deltax));
                }

                normalizedVarIncrement += skewComp + convComp;
                // for extreme values of parameters the 'squared' local vol becomes
                // negative !!! Therefore, we need to fudge things here...
                normalizedVarIncrement = Maths::max(0.0, normalizedVarIncrement);
            }

            double atmVarDiff = totalAtmVar - lastAtmVar;
            var += atmVarDiff * normalizedVarIncrement;
            return sqrt(var / tradYear);
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    // calculate strike ref
    double getStrikeRef(double tradYear) const {
        static int guessIdx = -1;
        if (strikeRefRatios.size() == 0){
            return strikeRef;
        }
        double ratio = strikeRefRatioInterp->valueWithGuess(tradYear, guessIdx);
        return (ratio * strikeRef);
    }

    CDoubleMatrix spotVolMatrixFromStrikes(const CDoubleArray& strikes) const {
        static const string method("VolHyperTrigPrimer::spotVolMatrixFromStrikes");
        try {
            DateTimeArray bbDates = backbone->getDates();
            CDoubleMatrix matrix(strikes.size(), bbDates.size());
            // loop over maturities
            for (int iMat = 0; iMat < bbDates.size(); iMat++) {
                double tradYear = timeMetric->yearFrac(baseDate, bbDates[iMat]);
                // loop over strikes
                for (int iStrike = 0; iStrike < strikes.size(); iStrike ++) {
                    // do somthing here
                    matrix[iStrike][iMat] = calcVol(strikes[iStrike],
                                                    tradYear);
                }
            }
            return matrix;
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }

protected:
    friend class VolHyperTrig;
    friend class VolHyperTrigDual;

    // everything is done, that is common to VolHyperTrig and VolHyperTrigDual
    // specific update methods are in the helper classes of VolHyperTrig and VolHyperTrigDual
    void updateParent() {
        static const string method("VolHyperTrigPrimer::updateParent");
        static const double ln100To90 = -log(0.9);
        try {
            // SKEW FACTORS
            // Interpolation of Skew Factors for Skew Term Structure
            // if the 1y ref date hasn't been found, need to insert it
            DoubleArray skewFactorsUsed;
            if (skewFactorTradYears.size() != skewFactors.size()){
                skewFactorsUsed.resize(skewFactorTradYears.size());
                int iSkewFactorBM = 0;
                int iNewSkewFactorBM = 0;
                while (iSkewFactorBM < skewFactorBMs.size()
                       && skewFactorDates[iSkewFactorBM] < dateRef){
                    skewFactorsUsed[iNewSkewFactorBM] = skewFactors[iSkewFactorBM];
                    ++iSkewFactorBM;
                    ++iNewSkewFactorBM;
                }
                skewFactorsUsed[iNewSkewFactorBM] = 1.0;
                ++iNewSkewFactorBM;
                while (iSkewFactorBM < skewFactorBMs.size()){
                    skewFactorsUsed[iNewSkewFactorBM] = skewFactors[iSkewFactorBM];
                    ++iSkewFactorBM;
                    ++iNewSkewFactorBM;
                }
            }
            else{
                skewFactorsUsed = skewFactors;
            }
            // interpolate the skew factors linearly
            LinearInterpolator interpolator;
            skewFactorInterp = Interpolator::InterpolantSP::constCast(
                interpolator.computeInterp(skewFactorTradYears, skewFactorsUsed));


            // CONVEXITY FACTORS
            // Interpolation of Convexity Factors for Convexity Term Structure
            // if the 1y ref date hasn't been found, need to insert it
            DoubleArray convexityFactorsUsed;
            if (convexityFactorTradYears.size() != convexityFactors.size()){
                convexityFactorsUsed.resize(convexityFactorTradYears.size());
                int iConvexityFactorBM = 0;
                int iNewConvexityFactorBM = 0;
                while (iConvexityFactorBM < convexityFactorBMs.size()
                       && convexityFactorDates[iConvexityFactorBM] < dateRef){
                    convexityFactorsUsed[iNewConvexityFactorBM] = convexityFactors[iConvexityFactorBM];
                    ++iConvexityFactorBM;
                    ++iNewConvexityFactorBM;
                }
                convexityFactorsUsed[iNewConvexityFactorBM] = 1.0;
                ++iNewConvexityFactorBM;
                while (iConvexityFactorBM < convexityFactorBMs.size()){
                    convexityFactorsUsed[iNewConvexityFactorBM] = convexityFactors[iConvexityFactorBM];
                    ++iConvexityFactorBM;
                    ++iNewConvexityFactorBM;
                }
            }
            else{
                convexityFactorsUsed = convexityFactors;
            }
            // interpolate the convexity factors linearly
            convexityFactorInterp = Interpolator::InterpolantSP::constCast(
                interpolator.computeInterp(convexityFactorTradYears, convexityFactorsUsed));


            // SKEW CELERITY FACTORS
            // Interpolation of Skew Celerity Factors for Skew Celerity Term Structure
            // if the 1y ref date hasn't been found, need to insert it
            DoubleArray skewCelerityFactorsUsed;
            if (skewCelerityFactorTradYears.size() != skewCelerityFactors.size()){
                skewCelerityFactorsUsed.resize(skewCelerityFactorTradYears.size());
                int iSkewCelerityFactorBM = 0;
                int iNewSkewCelerityFactorBM = 0;
                while (iSkewCelerityFactorBM < skewCelerityFactorBMs.size()
                       && skewCelerityFactorDates[iSkewCelerityFactorBM] < dateRef){
                    skewCelerityFactorsUsed[iNewSkewCelerityFactorBM] = skewCelerityFactors[iSkewCelerityFactorBM];
                    ++iSkewCelerityFactorBM;
                    ++iNewSkewCelerityFactorBM;
                }
                skewCelerityFactorsUsed[iNewSkewCelerityFactorBM] = 1.0;
                ++iNewSkewCelerityFactorBM;
                while (iSkewCelerityFactorBM < skewCelerityFactorBMs.size()){
                    skewCelerityFactorsUsed[iNewSkewCelerityFactorBM] = skewCelerityFactors[iSkewCelerityFactorBM];
                    ++iSkewCelerityFactorBM;
                    ++iNewSkewCelerityFactorBM;
                }
            }
            else{
                skewCelerityFactorsUsed = skewCelerityFactors;
            }
            // interpolate the skew Celerity factors linearly
            skewCelerityFactorInterp = Interpolator::InterpolantSP::constCast(
                interpolator.computeInterp(skewCelerityFactorTradYears, skewCelerityFactorsUsed));


            // CONVEXITY CELERITY FACTORS
            // Interpolation of Convexity Factors for Convexity Term Structure
            // if the 1y ref date hasn't been found, need to insert it
            DoubleArray convexityCelerityFactorsUsed;
            if (convexityCelerityFactorTradYears.size() != convexityCelerityFactors.size()){
                convexityCelerityFactorsUsed.resize(convexityCelerityFactorTradYears.size());
                int iConvexityCelerityFactorBM = 0;
                int iNewConvexityCelerityFactorBM = 0;
                while (iConvexityCelerityFactorBM < convexityCelerityFactorBMs.size()
                       && convexityCelerityFactorDates[iConvexityCelerityFactorBM] < dateRef){
                    convexityCelerityFactorsUsed[iNewConvexityCelerityFactorBM] = convexityCelerityFactors[iConvexityCelerityFactorBM];
                    ++iConvexityCelerityFactorBM;
                    ++iNewConvexityCelerityFactorBM;
                }
                convexityCelerityFactorsUsed[iNewConvexityCelerityFactorBM] = 1.0;
                ++iNewConvexityCelerityFactorBM;
                while (iConvexityCelerityFactorBM < convexityCelerityFactorBMs.size()){
                    convexityCelerityFactorsUsed[iNewConvexityCelerityFactorBM] = convexityCelerityFactors[iConvexityCelerityFactorBM];
                    ++iConvexityCelerityFactorBM;
                    ++iNewConvexityCelerityFactorBM;
                }
            }
            else{
                convexityCelerityFactorsUsed = convexityCelerityFactors;
            }
            // interpolate the convexity factors linearly
            convexityCelerityFactorInterp = Interpolator::InterpolantSP::constCast(
                interpolator.computeInterp(convexityCelerityFactorTradYears, convexityCelerityFactorsUsed));


            // preparation in order to back out parameters A,B,C,D
            // however, these parameters are then only backed out in the update methods of the helper functions
            double dAtmVol_dx = skew / ln100To90; // scale skew
            double d2AtmVol_dx2 = convexity / (ln100To90 * ln100To90); // scale convexity

            for (int iTradYear = 0; iTradYear < tradYears.size(); ++iTradYear){
                double thisTradYear = tradYears[iTradYear];
                // skew(T) = skew(1yr) * SCF(T) / sqrt(T/1yr)
                skewVector[iTradYear] = dAtmVol_dx * skewFactorInterp->value(thisTradYear)
                                    / sqrt(thisTradYear / tradYearRef);
                // convexity(T) = convexity(1yr) * CCF(T) / sqrt(T/1yr)
                convexityVector[iTradYear] = d2AtmVol_dx2 * convexityFactorInterp->value(thisTradYear)
                                        / pow(thisTradYear / tradYearRef, convexityPower);

                skewCelerityVector[iTradYear] = skewCelerityFactorInterp->value(thisTradYear) * skewCelerity /
                    sqrt(thisTradYear / tradYearRef);

                convexityCelerityVector[iTradYear] = convexityCelerityFactorInterp->value(thisTradYear) * convexityCelerity /
                    sqrt(thisTradYear / tradYearRef);

                double atmVar = atmVars[iTradYear];
                double dAtmVar_dx = 2.0 * thisTradYear * sqrt(atmVar/thisTradYear) * skewVector[iTradYear];
                double d2AtmVar_dx2 = 2.0  * thisTradYear *
                    (skewVector[iTradYear] * skewVector[iTradYear] +
                        sqrt(atmVar / thisTradYear) * convexityVector[iTradYear]);
                vectorLinEq1[iTradYear] = dAtmVar_dx * atmVar;
                vectorLinEq2[iTradYear] = d2AtmVar_dx2 * atmVar * atmVar;
                // another fudge in order to make the calibrator work ...
                if ((iTradYear > 0) &&
                        (Maths::isNegative(vectorLinEq2[iTradYear] - vectorLinEq2[iTradYear-1]))) {
                    vectorLinEq2[iTradYear] = vectorLinEq2[iTradYear-1];
                }
            }
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    // validate & initialize all the fields ...
    void validatePop2Object(){
        static const string method("VolHyperTrigPrimer::validatePop2Object");
        try{
            // pade coefficients
            if (approxPade == "EXACT") {
                usePade = false;
            } else {
                // Figure out the number of pade coefficients
                char* endPtr;
                int nbPadeCoeffs = (int)strtol(approxPade.c_str(), &endPtr, 10 /* base 10 */);
                if (*endPtr != '\0') {
                    throw ModelException(method, "Failed to convert approxPade string ("+
                        string(approxPade)+ ") to an int");
                }
                // Create approximants
                usePade = true;
                fastLogCosh = PadeApproxLogCosh(nbPadeCoeffs);
                fastAtanExp = PadeApproxAtanExp(nbPadeCoeffs);
            }

            // Check strike Refs
            Maths::checkPositive(strikeRef, "strikeRef");

            int nbStrikeRefRatios = strikeRefRatios.size();
            if (nbStrikeRefRatios != strikeRefRatioBMs.size()){
                throw ModelException(method,
                                     "The number of reference strike ratios ("
                                     + Format::toString(nbStrikeRefRatios)
                                     + ") must be equal to the number of reference strike ratio benchmarks ("
                                     + Format::toString(strikeRefRatioBMs.size())
                                     + ")");
            }
            Maths::checkPositive(strikeRefRatios, "strikeRefRatios");
            strikeRefRatioDates.resize(nbStrikeRefRatios);
            strikeRefRatioTradYears.resize(nbStrikeRefRatios);

            // Check skew factors
            int nbSkewFactors = skewFactors.size();
            if (nbSkewFactors != skewFactorBMs.size()){
                throw ModelException(method,
                                     "The number of skew factors ("
                                     + Format::toString(skewFactors.size())
                                     + ") must be equal to the number of skew factor benchmarks ("
                                     + Format::toString(skewFactorBMs.size())
                                     + ")");
            }
            skewFactorDates.resize(nbSkewFactors);
            skewFactorTradYears.resize(nbSkewFactors);

            // Check convexity factors
            int nbConvexityFactors = convexityFactors.size();
            if (nbConvexityFactors != convexityFactorBMs.size()){
                throw ModelException(method,
                                     "The number of convexity factors ("
                                     + Format::toString(convexityFactors.size())
                                     + ") must be equal to the number of convexity factor benchmarks ("
                                     + Format::toString(convexityFactorBMs.size())
                                     + ")");
            }
            convexityFactorDates.resize(nbConvexityFactors);
            convexityFactorTradYears.resize(nbConvexityFactors);

            // Check skew celerity factors
            int nbSkewCelerityFactors = skewCelerityFactors.size();
            if (nbSkewCelerityFactors != skewCelerityFactorBMs.size()){
                throw ModelException(method,
                                     "The number of skew celerity factors ("
                                     + Format::toString(skewCelerityFactors.size())
                                     + ") must be equal to the number of skew celerity factor benchmarks ("
                                     + Format::toString(skewCelerityFactorBMs.size())
                                     + ")");
            }
            skewCelerityFactorDates.resize(nbSkewCelerityFactors);
            skewCelerityFactorTradYears.resize(nbSkewCelerityFactors);

            // Check convexity factors
            int nbConvexityCelerityFactors = convexityCelerityFactors.size();
            if (nbConvexityCelerityFactors != convexityCelerityFactorBMs.size()){
                throw ModelException(method,
                                     "The number of convexity elerity factors ("
                                     + Format::toString(convexityCelerityFactors.size())
                                     + ") must be equal to the number of convexity celerity factor benchmarks ("
                                     + Format::toString(convexityCelerityFactorBMs.size())
                                     + ")");
            }
            convexityCelerityFactorDates.resize(nbConvexityCelerityFactors);
            convexityCelerityFactorTradYears.resize(nbConvexityCelerityFactors);

            // extract information from backbone
            ExpiryArrayConstSP expiries = backbone->getExpiries();
            timeMetric = backbone->getTimeMetric();

            // Set useMarkedDate=true if markedDate present and backbone consists entirely of fixed dates
            bool useMarkedDate = true;
            if (markedDate.empty()) {
                useMarkedDate = false;
            }

            BenchmarkDate* date;
            for (int i=0; i<expiries->size(); i++) {
                date = dynamic_cast<BenchmarkDate *>((*expiries)[i].get());
                if (!date) {
                    useMarkedDate = false;
                }
            }

            VolSurfaceSP rollBackbone(VolSurfaceSP(dynamic_cast<VolSurface *>(backbone->clone())));
            if (useMarkedDate && markedDate.getDate() < backbone->getBaseDate().getDate()) {
                // Roll the backbone so that it's baseDate matches last marked date
                // Needed to compute the atmVars[] array
                int interval = markedDate.daysDiff(rollBackbone->getBaseDate());
                ThetaSP rollToMark(new Theta(interval, HolidaySP(Holiday::noHolidays())));
                rollBackbone->sensShift(rollToMark.get());

                // temporary fix for pricing at EOD - roll time of markedDate to pricing time 
                // Strategic fix is to allow partial day shifts in Theta
                markedDate = DateTime(markedDate.getDate(), rollBackbone->getBaseDate().getTime());

                baseDate = markedDate;

                //DateTime::Interval offset = markedDate.subtract(rollBackbone->getBaseDate());
                //ThetaSP rollToMark(new Theta(offset, HolidaySP(Holiday::noHolidays())));
            } else {
                // markedDate is not used if:
                // 1. markedDate is not present or 
                // 2. backbone contains at least one rolling period or
                // 3. markedDate is in the future
                baseDate = backbone->getBaseDate();
            }

            // Set dates to rolled backbone dates
            dates = rollBackbone->getDates();
            tenorTimeInt = DateTime::timeConvert(tenorTime);

            int iStrikeRefRatio = 0;
            for (; iStrikeRefRatio < strikeRefRatioBMs.size(); ++iStrikeRefRatio){
                MaturityPeriod period(strikeRefRatioBMs[iStrikeRefRatio]);
                strikeRefRatioDates[iStrikeRefRatio] = DateTime(period.toDate(baseDate).getDate(), tenorTimeInt);
                if (iStrikeRefRatio > 0
                    && strikeRefRatioDates[iStrikeRefRatio - 1] >= strikeRefRatioDates[iStrikeRefRatio]){
                    throw ModelException(method,
                                         "The strike ref ratio benchmarks must be increasing.\nCompare the "
                                         + Format::toString(iStrikeRefRatio)
                                         + "-th benchmark date ("
                                         + strikeRefRatioDates[iStrikeRefRatio - 1].toString()
                                         + ") with the "
                                         + Format::toString(iStrikeRefRatio + 1)
                                         + "-th benchmark date ("
                                         + strikeRefRatioDates[iStrikeRefRatio].toString()
                                         + ")");
                }
                strikeRefRatioTradYears[iStrikeRefRatio]
                    = timeMetric->yearFrac(baseDate, strikeRefRatioDates[iStrikeRefRatio]);
            }
            // create a linear interpolator for strike refs
            if (strikeRefRatios.size() > 0){
                LinearInterpolator inter;
                strikeRefRatioInterp
                    = LinearInterpolantConstSP::dynamicCast(
                            inter.computeInterp(strikeRefRatioTradYears,
                                                strikeRefRatios));
            }

            int iDate = 0;
            MaturityPeriod periodRef(tenorRef);
            dateRef = DateTime(periodRef.toDate(baseDate).getDate(), tenorTimeInt);

            int nbDates = dates.size();

            dateRefIdx=0;
            bool tenorRefFound = false;
            for (iDate = 0; iDate < dates.size(); ++iDate){
                if (dateRef == dates[iDate]){
                    tenorRefFound = true;
                    break;
                }
                if (dateRef < dates[iDate]) {
                    // dateRef not on backbone, gone past it
                    break;
                }
                dateRefIdx++;
            }

            // tenorRefFound = true: dateRefIdx is the index of the dates array corresponding to dateRef
            // tenorRefFound = false: dateRefIdx is the required insertion point

            if (!tenorRefFound) {
                nbDates ++;
            }

            atmVars.resize(nbDates);
            tradYears.resize(nbDates);

            // Calculate tradYears and atmVars.
            // Insert interpolated vol at dateRef if not already on backbone
            DateTimeArray newDates(nbDates);
            DateTime curDate;
            for (iDate = 0; iDate < nbDates; ++iDate){
                if (tenorRefFound || (iDate<dateRefIdx)) {
                    curDate = dates[iDate];
                } else if (iDate==dateRefIdx) {
                    curDate = dateRef;
                } else {
                    curDate = dates[iDate-1];
                }
                newDates[iDate] = curDate;
                tradYears[iDate] = timeMetric->yearFrac(baseDate, curDate);
                CVolProcessedBSSP atmVolCurve =
                    CVolProcessedBSSP(rollBackbone->getProcessedVol(getStrikeRef(tradYears[iDate])));
                atmVars[iDate] = atmVolCurve->CalcVar(baseDate, curDate);
            }
            dates = newDates;

             // resize A,B,C,D
            backboneSkewCelerities.resize(nbDates);
            backboneSkewParams.resize(nbDates);
            backboneConvexityCelerities.resize(nbDates);
            backboneConvexityParams.resize(nbDates);

            // check tradYears and vars are increasing
            if (!Maths::isPositive(tradYears[0])){
                throw ModelException(method,
                                     "there is zero trading time between today ("
                                     + baseDate.toString()
                                     + ") and the 1st benchmark date ("
                                     + dates[0].toString()
                                     + ")");
            }
            if (!Maths::isPositive(atmVars[0])){
                throw ModelException(method,
                                     "there is zero atm variance between today ("
                                     + baseDate.toString()
                                     + ") and the 1st benchmark date ("
                                     + dates[0].toString()
                                     + ")");
            }
            for (iDate = 1; iDate < nbDates; ++iDate){
                // trad years
                if (!Maths::isPositive(tradYears[iDate] - tradYears[iDate - 1])){
                    throw ModelException(method,
                                         "there is zero trading time between the "
                                         + Format::toString(iDate)
                                         + "-th benchmark date ("
                                         + dates[iDate].toString()
                                         + ") and the "
                                         + Format::toString(iDate + 1)
                                         + "-th benchmark date ("
                                         + dates[iDate].toString()
                                         + ")");
                }
                // variances
                double fwdVar = atmVars[iDate] - atmVars[iDate - 1];
                if (!Maths::isPositive(fwdVar)){
                    throw ModelException(method,
                                         "the forward variance ("
                                         + Format::toString(fwdVar)
                                         + ") between the "
                                         + Format::toString(iDate)
                                         + "-th benchmark date ("
                                         + dates[iDate].toString()
                                         + ") and the "
                                         + Format::toString(iDate + 1)
                                         + "-th benchmark date ("
                                         + dates[iDate].toString()
                                         + ") is non-positive");
                }
            }

            // calc vols at ref point
            atmVarRef = atmVars[dateRefIdx];
            atmVolRef = sqrt(atmVars[dateRefIdx] / tradYears[dateRefIdx]);
            tradYearRef = tradYears[dateRefIdx];

            // calculates skew factor dates and years
            // checks that the skew factor BMs are increasing
            // checks that if the 1y ref is among the ske factor BMs then
            // the corresponding skew factor == 1
            // Dummy try/catch to avoid vc6.opt crash
			try {
                calcFactorTerms(skewFactorBMs,
                                skewFactors,
                                skewFactorDates,
                                skewFactorTradYears,
                                "skew factor");
			} catch (...) { throw; }
// do the same for the convexity factors
            calcFactorTerms(convexityFactorBMs,
                            convexityFactors,
                            convexityFactorDates,
                            convexityFactorTradYears,
                            "convexity factor");
            // do the same for the skew celerity factors
            calcFactorTerms(skewCelerityFactorBMs,
                            skewCelerityFactors,
                            skewCelerityFactorDates,
                            skewCelerityFactorTradYears,
                            "skew celerity factor");
            // do the same for the convexity celerity factors
            calcFactorTerms(convexityCelerityFactorBMs,
                            convexityCelerityFactors,
                            convexityCelerityFactorDates,
                            convexityCelerityFactorTradYears,
                            "convexity celerity factor");

            // resize some more things ...
            skewVector.resize(nbDates);
            convexityVector.resize(nbDates);
            vectorLinEq1.resize(nbDates);
            vectorLinEq2.resize(nbDates);
            skewCelerityVector.resize(nbDates);
            convexityCelerityVector.resize(nbDates);

            // validate volgrid
            volGrid.get()->validatePop2Object();

        } catch(exception& e) {
           throw ModelException(e, method);
        }
    }

    // MANDATORY / OPTINAL FIELDS

    // registered fields from VHT resp VHTdual passed to the constructor of VolHyperTrigPrimer
    string                      approxPade;

    // Input strikeRef data
    string                      tenorRef;
    StringArray                 strikeRefRatioBMs;
    DoubleArray                 strikeRefRatios;
    double                      strikeRef;

    // Input skew data
    StringArray                 skewFactorBMs;
    DoubleArray                 skewFactors;
    double                      skew;

    // Input convexity data
    StringArray                 convexityFactorBMs;
    DoubleArray                 convexityFactors;
    double                      convexity;

    // Input skew celerity data
    StringArray                 skewCelerityFactorBMs;
    DoubleArray                 skewCelerityFactors;
    double                      skewCelerity;

    // Input convexity celerity data
    StringArray                 convexityCelerityFactorBMs;
    DoubleArray                 convexityCelerityFactors;
    double                      convexityCelerity;

    string                      tenorTime;
    bool                        useTimeDependentCelerities;
    double                      convexityPower;
    double                      skewBackboneScaling;
    double                      convexityBackboneScaling;

    // TRANSIENT FIELDS
    VolSurfaceConstSP           backbone;

    // additional fields
    DateTime                    baseDate;
    TimeMetricConstSP           timeMetric;
    DateTimeArray               dates;
    DoubleArray                 atmVars;

    // Internal strikeRef term structures
    DateTimeArray               strikeRefRatioDates;
    DoubleArray                 strikeRefRatioTradYears;
    LinearInterpolantConstSP    strikeRefRatioInterp;

    // Internal skew term structures
    DateTimeArray               skewFactorDates;
    DoubleArray                 skewFactorTradYears;
    DoubleArray                 skewVector;

    // Internal convexity term structures
    DateTimeArray               convexityFactorDates;
    DoubleArray                 convexityFactorTradYears;
    DoubleArray                 convexityVector;

    // Internal skew celerity term structures
    DateTimeArray               skewCelerityFactorDates;
    DoubleArray                 skewCelerityFactorTradYears;
    DoubleArray                 skewCelerityVector;

    // Internal convexity celerity term structures
    DateTimeArray               convexityCelerityFactorDates;
    DoubleArray                 convexityCelerityFactorTradYears;
    DoubleArray                 convexityCelerityVector;

    DoubleArray                 vectorLinEq1;
    DoubleArray                 vectorLinEq2;

    DoubleArray                 tradYears;
    int                         tenorTimeInt;

    // pade approximation
    bool                        usePade;
    PadeApproxLogCosh           fastLogCosh;
    PadeApproxAtanExp           fastAtanExp;

    // Skew & convexity factor interpolants
    Interpolator::InterpolantSP skewFactorInterp;           // skew factor term structure interpolating polynomial
    Interpolator::InterpolantSP convexityFactorInterp;      // convexity factor term structure interpolating polynomial

    // Skew & convexity celerity factor interpolants
    Interpolator::InterpolantSP skewCelerityFactorInterp;
    Interpolator::InterpolantSP convexityCelerityFactorInterp;


    // ref point specifics
    double            atmVolRef;            // atm vol
    double            atmVarRef;            // atm var
    DateTime          dateRef;              // date at ref tenor
    int               dateRefIdx;           // new field
    double            tradYearRef;

    // A,B,C,D -- resize them here and assign values in constructor of Helper Classes
    DoubleArray             backboneSkewCelerities;
    DoubleArray             backboneSkewParams;
    DoubleArray             backboneConvexityCelerities;
    DoubleArray             backboneConvexityParams;

    // for VHT speed up
    bool                    useGrid;
    VolGridSP               volGrid;

    // for cutoff
    VarSwapBasisSP          basis;
    DateTime                markedDate;

        static void load(CClassSP& clazz){
        REGISTER(VolHyperTrigPrimer, clazz);
        SUPERCLASS(CObject);
        FIELD(approxPade, "");
        FIELD(tenorRef, "");

        FIELD(strikeRefRatioBMs, "");
        FIELD(strikeRefRatios, "");
        FIELD(strikeRef, "");

        FIELD(skewFactorBMs, "");
        FIELD(skewFactors, "");
        FIELD(skew, "");

        FIELD(convexityFactorBMs, "");
        FIELD(convexityFactors, "");
        FIELD(convexity, "");

        FIELD(skewCelerityFactorBMs, "");
        FIELD(skewCelerityFactors, "");
        FIELD(skewCelerity, "");

        FIELD(convexityCelerityFactorBMs, "");
        FIELD(convexityCelerityFactors, "");
        FIELD(convexityCelerity, "");

        FIELD(tenorTime, "");
        FIELD(useTimeDependentCelerities, "");
        FIELD_MAKE_OPTIONAL(useTimeDependentCelerities);
        FIELD(convexityPower, "");
        FIELD_MAKE_OPTIONAL(convexityPower);
        FIELD(skewBackboneScaling, "");
        FIELD_MAKE_OPTIONAL(skewBackboneScaling);
        FIELD(convexityBackboneScaling, "");
        FIELD_MAKE_OPTIONAL(convexityBackboneScaling);

        // transient fields
        FIELD(backbone, "");
        FIELD_MAKE_TRANSIENT(backbone);

        FIELD(baseDate, "");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(timeMetric, "");
        FIELD_MAKE_TRANSIENT(timeMetric);
        FIELD(dates, "");
        FIELD_MAKE_TRANSIENT(dates);
        FIELD(atmVars, "");
        FIELD_MAKE_TRANSIENT(atmVars);

        FIELD(strikeRefRatioDates, "");
        FIELD_MAKE_TRANSIENT(strikeRefRatioDates);
        FIELD(strikeRefRatioTradYears, "");
        FIELD_MAKE_TRANSIENT(strikeRefRatioTradYears);
        FIELD(strikeRefRatioInterp, "");
        FIELD_MAKE_TRANSIENT(strikeRefRatioInterp);

        FIELD(skewFactorDates, "");
        FIELD_MAKE_TRANSIENT(skewFactorDates);
        FIELD(skewFactorTradYears, "");
        FIELD_MAKE_TRANSIENT(skewFactorTradYears);
        FIELD(skewVector, "");
        FIELD_MAKE_TRANSIENT(skewVector);

        FIELD(convexityFactorDates, "");
        FIELD_MAKE_TRANSIENT(convexityFactorDates);
        FIELD(convexityFactorTradYears, "");
        FIELD_MAKE_TRANSIENT(convexityFactorTradYears);
        FIELD(convexityVector, "");
        FIELD_MAKE_TRANSIENT(convexityVector);

        FIELD(skewCelerityFactorDates, "");
        FIELD_MAKE_TRANSIENT(skewCelerityFactorDates);
        FIELD(skewCelerityFactorTradYears, "");
        FIELD_MAKE_TRANSIENT(skewCelerityFactorTradYears);
        FIELD(skewCelerityVector, "");
        FIELD_MAKE_TRANSIENT(skewCelerityVector);

        FIELD(convexityCelerityFactorDates, "");
        FIELD_MAKE_TRANSIENT(convexityCelerityFactorDates);
        FIELD(convexityCelerityFactorTradYears, "");
        FIELD_MAKE_TRANSIENT(convexityCelerityFactorTradYears);
        FIELD(convexityCelerityVector, "");
        FIELD_MAKE_TRANSIENT(convexityCelerityVector);

        FIELD(vectorLinEq1, "");
        FIELD_MAKE_TRANSIENT(vectorLinEq1);
        FIELD(vectorLinEq2, "");
        FIELD_MAKE_TRANSIENT(vectorLinEq2);

        FIELD(tradYears, "");
        FIELD_MAKE_TRANSIENT(tradYears);
        FIELD(tenorTimeInt, "");
        FIELD_MAKE_TRANSIENT(tenorTimeInt);

        FIELD(usePade, "");
        FIELD_MAKE_TRANSIENT(usePade);
        FIELD(fastLogCosh, "");
        FIELD_MAKE_TRANSIENT(fastLogCosh);
        FIELD(fastAtanExp, "");
        FIELD_MAKE_TRANSIENT(fastAtanExp);

        FIELD(skewFactorInterp, "");
        FIELD_MAKE_TRANSIENT(skewFactorInterp);

        FIELD(convexityFactorInterp, "");
        FIELD_MAKE_TRANSIENT(convexityFactorInterp);

        FIELD(skewCelerityFactorInterp, "");
        FIELD_MAKE_TRANSIENT(skewCelerityFactorInterp);

        FIELD(convexityCelerityFactorInterp, "");
        FIELD_MAKE_TRANSIENT(convexityCelerityFactorInterp);

        FIELD(atmVolRef, "");
        FIELD_MAKE_TRANSIENT(atmVolRef);
        FIELD(atmVarRef, "");
        FIELD_MAKE_TRANSIENT(atmVarRef);
        FIELD(dateRef, "");
        FIELD_MAKE_TRANSIENT(dateRef);
        FIELD(dateRefIdx, "");
        FIELD_MAKE_TRANSIENT(dateRefIdx);
        FIELD(markedDate, "");
        FIELD_MAKE_TRANSIENT(markedDate);
        FIELD(tradYearRef, "");
        FIELD_MAKE_TRANSIENT(tradYearRef);

        FIELD(backboneSkewCelerities, "");
        FIELD_MAKE_TRANSIENT(backboneSkewCelerities);
        FIELD(backboneSkewParams, "");
        FIELD_MAKE_TRANSIENT(backboneSkewParams);
        FIELD(backboneConvexityCelerities, "");
        FIELD_MAKE_TRANSIENT(backboneConvexityCelerities);
        FIELD(backboneConvexityParams, "");
        FIELD_MAKE_TRANSIENT(backboneConvexityParams);

        FIELD(useGrid, "");
        FIELD_MAKE_TRANSIENT(useGrid);
        FIELD(volGrid, "");
        FIELD_MAKE_TRANSIENT(volGrid);
        FIELD(basis, "");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(basis);
    }
};

typedef smartPtr<VolHyperTrigPrimer> VolHyperTrigPrimerSP;

CClassConstSP const VolHyperTrigPrimer::TYPE =
CClass::registerClassLoadMethod("VolHyperTrigPrimer", typeid(VolHyperTrigPrimer), VolHyperTrigPrimer::load);

const string VolHyperTrigPrimer::VolGrid::GRID_TYPE_NONE    =   "none";
const string VolHyperTrigPrimer::VolGrid::GRID_TYPE_LINEAR  =   "linear";
const string VolHyperTrigPrimer::VolGrid::GRID_TYPE_SPLINE  =   "spline";
const string VolHyperTrigPrimer::VolGrid::SPEEDUP_METHOD_GRID  =   "grid";

CClassConstSP const VolHyperTrigPrimer::VolGrid::TYPE =
CClass::registerClassLoadMethod("VolHyperTrigPrimer::VolGrid", typeid(VolHyperTrigPrimer::VolGrid), VolHyperTrigPrimer::VolGrid::load);

/** skewParam parameterised CVolBase that supports BS and DVF views of the world */
class VolHyperTrig: public VolBaseParam,
                    virtual public IVolProcessed,
                    virtual public IVolatilityBS,
                    virtual public IVolatilityDVF,
                    virtual public DeltaSurface::IShift,
                    public virtual Calibrator::IAdjustable,
                    virtual public VHTSkewFactorsPointwise::IShift,
                    virtual public VHTSkewCelerityFactorsPointwise::IShift,
                    virtual public VHTConvexityFactorsPointwise::IShift,
                    virtual public VHTConvexityCelerityFactorsPointwise::IShift,
					virtual public ITweakableWithRespectTo<VolAtmParallelConstConvx>,
					virtual public ITweakableWithRespectTo<VolAtmPointwiseConstConvx> {
private:
    // helper class for correct interpretation of input parameters
    class HelperFull: public VolHyperTrigPrimer{
    public:
        static CClassConstSP const TYPE;

        // constructor
        HelperFull(): VolHyperTrigPrimer(TYPE) {}
        HelperFull(const VolHyperTrig&  vht, VolSurfaceConstSP backbone):
        VolHyperTrigPrimer(TYPE,
                          // registered fields
                          vht.approxPade,
                          vht.tenorRef,
                          vht.strikeRefRatioBMs,
                          vht.strikeRefRatios,
                          vht.strikeRef,
                          vht.skewFactorBMs,
                          vht.skewFactors,
                          vht.skew,
                          vht.convexityFactorBMs,
                          vht.convexityFactors,
                          vht.convexity,
                          vht.skewCelerityFactorBMs,
                          vht.skewCelerityFactors,
                          vht.skewCelerity,
                          vht.convexityCelerityFactorBMs,
                          vht.convexityCelerityFactors,
                          vht.convexityCelerity,
                          vht.tenorTime,
                          vht.useTimeDependentCelerities,
                          vht.convexityPower,
                          vht.skewBackboneScaling,
                          vht.convexityBackboneScaling,
                          // additional fields
                          backbone,
                          vht.gridType,
                          vht.numGrid,
                          vht.TruncationStd,
                          vht.useGrid,
                          vht.basis,
                          vht.markedDate) {

            static const string method("VolHyperTrig::HelperFull::HelperFull");
            // step I:   general update already in constructor of VolHyperTrigPrimer via call to validatePop2Object
            // step II: specific update function of Helper class (diff for VHT and VHTDual)
            update(skew,
                   convexity,
                   skewCelerity,
                   convexityCelerity,
                   skewFactors,
                   convexityFactors,
                   skewCelerityFactors,
                   convexityCelerityFactors);
        }

        /** Update 4 term structures */
        void update(double skewMod,
                    double convexityMod,
                    double skewCelerityMod,
                    double convexityCelerityMod,
                    const DoubleArray& skewFactorsMod,
                    const DoubleArray& convexityFactorsMod,
                    const DoubleArray& skewCelerityFactorsMod,
                    const DoubleArray& convexityCelerityFactorsMod) {

            static const string method("VolHyperTrig::HelperFull::update");

            try{
                // step I:  general update -> delegated to VolHyperTrigPrimer (the same for VHT and VHTdual)
                skew = skewMod;
                convexity = convexityMod;
                skewCelerity = skewCelerityMod;
                convexityCelerity = convexityCelerityMod;

                skewFactors = skewFactorsMod;
                convexityFactors = convexityFactorsMod;
                skewCelerityFactors = skewCelerityFactorsMod;
                convexityCelerityFactors = convexityCelerityFactorsMod;

                updateParent();

                // step II: specific update for VHT
                // note C(1yr) = skewCelerity and D(1yr) = convexityCelerity
                backboneSkewCelerities[dateRefIdx] = skewCelerity;
                backboneConvexityCelerities[dateRefIdx] = convexityCelerity;
                double atmVolRef = sqrt (atmVars[0] / tradYears[0]);

                backboneSkewCelerities[0] = skewCelerityVector[0];
                backboneConvexityCelerities[0] = convexityCelerityVector[0];
                backboneSkewParams[0] = Maths::isZero(backboneSkewCelerities[0]) ? 1.0
                    : 2.0 * vectorLinEq1[0] / (backboneSkewCelerities[0] * atmVars[0]*atmVars[0]);

                backboneConvexityParams[0] = Maths::isZero(backboneConvexityCelerities[0]) ? 1.0
                    : 3.0 * vectorLinEq2[0] /
                        (backboneConvexityCelerities[0]*backboneConvexityCelerities[0] * atmVars[0]*atmVars[0]*atmVars[0]);

                for (int i = 1; i < tradYears.size(); i++) {
                    double thisAtmVar = atmVars[i];
                    double prevAtmVar = atmVars[i-1];
                    atmVolRef = sqrt (atmVars[i] / tradYears[i]);

                    backboneSkewCelerities[i] = skewCelerityVector[i];
                    backboneConvexityCelerities[i] = convexityCelerityVector[i];
                    backboneSkewParams[i] = Maths::isZero(backboneSkewCelerities[i]) ? 1.0
                        : 2.0 * (vectorLinEq1[i] - vectorLinEq1[i-1])
                            / (backboneSkewCelerities[i] * (thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar));
                    backboneConvexityParams[i] = Maths::isZero(backboneConvexityCelerities[i]) ? 1.0
                        : 3.0 * (vectorLinEq2[i] - vectorLinEq2[i-1])
                            / (backboneConvexityCelerities[i]*backboneConvexityCelerities[i]*
                            (thisAtmVar*thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar*prevAtmVar));
                }

            } catch(exception& e) {
                    ModelException(e, "VolHyperTrig::update");
            }
        }

        void update(double             skewMod,
                    double             convexityMod,
                    DoubleArray        skewFactorsMod,
                    DoubleArray        convexityFactorsMod,
                    double             skewCelerity,
                    double             convexityCelerity) {
            static const string method("VolHyperTrig::HelperFull::updateFull");
            try{
                // step I:  general update -> delegated to VolHyperTrigPrimer (the same for VHT and VHTdual)
                skew = skewMod;
                convexity = convexityMod;
                skewFactors = skewFactorsMod;
                convexityFactors = convexityFactorsMod;
                updateParent();

                // step II: specific update for VHT
                // note C(1yr) = skewCelerity and D(1yr) = convexityCelerity
                backboneSkewCelerities[dateRefIdx] = skewCelerity;
                backboneConvexityCelerities[dateRefIdx] = convexityCelerity;

                if(useTimeDependentCelerities) {
                    // C and D are time dependent and normalized by ATM vols
                    // A and B are constant and not normalized by ATM vols

                    // note A(1yr) = skewParam and B(1yr) = convexityParam
                    double atmVarRef = atmVars[dateRefIdx];
                    double atmVarRefMinus1 = atmVars[dateRefIdx-1];
                    double atmVolRef = sqrt (atmVars[dateRefIdx] / tradYears[dateRefIdx]);

                    // normalization has to be done only once, for skewCelerity at dateRefIdx
                    backboneSkewCelerities[dateRefIdx] = skewCelerity * pow(atmVolRef, skewBackboneScaling);
                    backboneConvexityCelerities[dateRefIdx] = convexityCelerity * pow(atmVolRef, convexityBackboneScaling / 2.0);

                    // back out skewParam A and convexityParam B
                    double skewParam = 2.0 * (vectorLinEq1[dateRefIdx] - vectorLinEq1[dateRefIdx-1])
                    / (backboneSkewCelerities[dateRefIdx] * (atmVarRef*atmVarRef - atmVarRefMinus1*atmVarRefMinus1));
                    double convexityParam = 3.0 * (vectorLinEq2[dateRefIdx] - vectorLinEq2[dateRefIdx-1])
                    / (backboneConvexityCelerities[dateRefIdx] * backboneConvexityCelerities[dateRefIdx]
                    * (atmVarRef*atmVarRef*atmVarRef - atmVarRefMinus1*atmVarRefMinus1*atmVarRefMinus1));

                    backboneSkewCelerities[0] = Maths::isZero(skewParam) ? 1.0
                        : 2.0 * vectorLinEq1[0] / (skewParam * atmVars[0]*atmVars[0]) *
                          pow(sqrt (atmVars[0] / tradYears[0]), skewBackboneScaling);
                    backboneConvexityCelerities[0] = Maths::isZero(convexityParam) ? 1.0
                        : sqrt(3.0 * vectorLinEq2[0] / (convexityParam * atmVars[0]*atmVars[0]*atmVars[0])) *
                          pow(sqrt (atmVars[0] / tradYears[0]), convexityBackboneScaling / 2.0);
                    backboneSkewParams[0] = skewParam;
                    backboneConvexityParams[0] = convexityParam;
                    for (int i = 1; i < tradYears.size(); i++) {
                        backboneSkewParams[i] = skewParam;
                        backboneConvexityParams[i] = convexityParam;
                        double thisAtmVar = atmVars[i];
                        double prevAtmVar = atmVars[i-1];
                        backboneSkewCelerities[i] = Maths::isZero(skewParam) ? 1.0
                            : 2.0 * (vectorLinEq1[i] - vectorLinEq1[i-1])
                                / (skewParam * (thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar))
                                * pow(sqrt (atmVars[i] / tradYears[i]), skewBackboneScaling);
                        backboneConvexityCelerities[i] = Maths::isZero(convexityParam) ? 1.0
                            : sqrt(3.0 * (vectorLinEq2[i] - vectorLinEq2[i-1])
                                / (convexityParam *
                                (thisAtmVar*thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar*prevAtmVar)))
                                * pow(sqrt (atmVars[i] / tradYears[i]), convexityBackboneScaling / 2.0);
                    }
                } else {
                    // A and B are time dependent and normalized by ATM vols
                    // C and D are constant and not normalized by ATM vols

                    double atmVolRef = sqrt (atmVars[0] / tradYears[0]);
                    backboneSkewParams[0] = Maths::isZero(skewCelerity) ? 1.0
                        : 2.0 * vectorLinEq1[0] / (skewCelerity * atmVars[0]*atmVars[0])
                            * pow(atmVolRef, skewBackboneScaling);
                    backboneConvexityParams[0] = Maths::isZero(convexityCelerity) ? 1.0
                        : 3.0 * vectorLinEq2[0] /
                            (convexityCelerity*convexityCelerity * atmVars[0]*atmVars[0]*atmVars[0])
                            * pow(atmVolRef, convexityBackboneScaling);

                    backboneSkewCelerities[0] = skewCelerity;
                    backboneConvexityCelerities[0] = convexityCelerity;

                    for (int i = 1; i < tradYears.size(); i++) {
                        backboneSkewCelerities[i] = skewCelerity;
                        backboneConvexityCelerities[i] = convexityCelerity;
                        double thisAtmVar = atmVars[i];
                        double prevAtmVar = atmVars[i-1];
                        atmVolRef = sqrt (atmVars[i] / tradYears[i]);
                        backboneSkewParams[i] = Maths::isZero(skewCelerity) ? 1.0
                            : 2.0 * (vectorLinEq1[i] - vectorLinEq1[i-1])
                                / (skewCelerity * (thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar))
                                * pow(atmVolRef, skewBackboneScaling);
                        backboneConvexityParams[i] = Maths::isZero(convexityCelerity) ? 1.0
                            : 3.0 * (vectorLinEq2[i] - vectorLinEq2[i-1])
                                / (convexityCelerity*convexityCelerity *
                                (thisAtmVar*thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar*prevAtmVar))
                                * pow(atmVolRef, convexityBackboneScaling);
                    }
                }

            } catch(exception& e) {
                    ModelException(e, "VolHyperTrig::update");
            }
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(HelperFull, clazz);
            SUPERCLASS(VolHyperTrigPrimer);
            EMPTY_SHELL_METHOD(defaultHelperFull);
        }
        static IObject* defaultHelperFull(){
            return new HelperFull();
        }
    };
    typedef smartPtr<HelperFull> HelperFullSP;

    // helper class for correct interpretation of input parameters
    class Helper: public VolHyperTrigPrimer{
    public:
        static CClassConstSP const TYPE;

        // constructor
        Helper(): VolHyperTrigPrimer(TYPE) {}
        Helper(const VolHyperTrig&  vht,
               VolSurfaceConstSP    backbone): VolHyperTrigPrimer(TYPE,
                                                                  // registered fields
                                                                  vht.approxPade,
                                                                  vht.tenorRef,
                                                                  vht.strikeRefRatioBMs,
                                                                  vht.strikeRefRatios,
                                                                  vht.strikeRef,
                                                                  vht.skewFactorBMs,
                                                                  vht.skewFactors,
                                                                  vht.skew,
                                                                  vht.convexityFactorBMs,
                                                                  vht.convexityFactors,
                                                                  vht.convexity,
                                                                  vht.tenorTime,
                                                                  vht.useTimeDependentCelerities,
                                                                  vht.convexityPower,
                                                                  vht.skewBackboneScaling,
                                                                  vht.convexityBackboneScaling,
                                                                  // additional fields
                                                                  backbone,
                                                                  vht.gridType,
                                                                  vht.numGrid,
                                                                  vht.TruncationStd,
                                                                  vht.useGrid,
                                                                  vht.basis,
                                                                  vht.markedDate) {
            static const string method("VolHyperTrig::Helper::Helper");
            // step I:   general update already in constructor of VolHyperTrigPrimer via call to validatePop2Object
            // step II: specific update function of Helper class (diff for VHT and VHTDual)
            update(skew,
                   convexity,
                   skewFactors,
                   convexityFactors,
                   vht.skewCelerity,
                   vht.convexityCelerity);
        }

        void update(double             skewMod,
                    double             convexityMod,
                    DoubleArray        skewFactorsMod,
                    DoubleArray        convexityFactorsMod,
                    double             skewCelerity,
                    double             convexityCelerity) {
            static const string method("VolHyperTrig::Helper::update");
            try{
                // step I:  general update -> delegated to VolHyperTrigPrimer (the same for VHT and VHTdual)
                skew = skewMod;
                convexity = convexityMod;
                skewFactors = skewFactorsMod;
                convexityFactors = convexityFactorsMod;
                updateParent();

                // step II: specific update for VHT
                // note C(1yr) = skewCelerity and D(1yr) = convexityCelerity
                backboneSkewCelerities[dateRefIdx] = skewCelerity;
                backboneConvexityCelerities[dateRefIdx] = convexityCelerity;

                if(useTimeDependentCelerities) {
                    // C and D are time dependent and normalized by ATM vols
                    // A and B are constant and not normalized by ATM vols

                    // note A(1yr) = skewParam and B(1yr) = convexityParam
                    double atmVarRef = atmVars[dateRefIdx];
                    double atmVarRefMinus1 = atmVars[dateRefIdx-1];
                    double atmVolRef = sqrt (atmVars[dateRefIdx] / tradYears[dateRefIdx]);

                    // normalization has to be done only once, for skewCelerity at dateRefIdx
                    backboneSkewCelerities[dateRefIdx] = skewCelerity * pow(atmVolRef, skewBackboneScaling);
                    backboneConvexityCelerities[dateRefIdx] = convexityCelerity * pow(atmVolRef, convexityBackboneScaling / 2.0);

                    // back out skewParam A and convexityParam B
                    double skewParam = 2.0 * (vectorLinEq1[dateRefIdx] - vectorLinEq1[dateRefIdx-1])
                    / (backboneSkewCelerities[dateRefIdx] * (atmVarRef*atmVarRef - atmVarRefMinus1*atmVarRefMinus1));
                    double convexityParam = 3.0 * (vectorLinEq2[dateRefIdx] - vectorLinEq2[dateRefIdx-1])
                    / (backboneConvexityCelerities[dateRefIdx] * backboneConvexityCelerities[dateRefIdx]
                    * (atmVarRef*atmVarRef*atmVarRef - atmVarRefMinus1*atmVarRefMinus1*atmVarRefMinus1));

                    backboneSkewCelerities[0] = Maths::isZero(skewParam) ? 1.0
                        : 2.0 * vectorLinEq1[0] / (skewParam * atmVars[0]*atmVars[0]) *
                          pow(sqrt (atmVars[0] / tradYears[0]), skewBackboneScaling);
                    backboneConvexityCelerities[0] = Maths::isZero(convexityParam) ? 1.0
                        : sqrt(3.0 * vectorLinEq2[0] / (convexityParam * atmVars[0]*atmVars[0]*atmVars[0])) *
                          pow(sqrt (atmVars[0] / tradYears[0]), convexityBackboneScaling / 2.0);
                    backboneSkewParams[0] = skewParam;
                    backboneConvexityParams[0] = convexityParam;
                    for (int i = 1; i < tradYears.size(); i++) {
                        backboneSkewParams[i] = skewParam;
                        backboneConvexityParams[i] = convexityParam;
                        double thisAtmVar = atmVars[i];
                        double prevAtmVar = atmVars[i-1];
                        backboneSkewCelerities[i] = Maths::isZero(skewParam) ? 1.0
                            : 2.0 * (vectorLinEq1[i] - vectorLinEq1[i-1])
                                / (skewParam * (thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar))
                                * pow(sqrt (atmVars[i] / tradYears[i]), skewBackboneScaling);
                        backboneConvexityCelerities[i] = Maths::isZero(convexityParam) ? 1.0
                            : sqrt(3.0 * (vectorLinEq2[i] - vectorLinEq2[i-1])
                                / (convexityParam *
                                (thisAtmVar*thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar*prevAtmVar)))
                                * pow(sqrt (atmVars[i] / tradYears[i]), convexityBackboneScaling / 2.0);
                    }
                } else {
                    // A and B are time dependent and normalized by ATM vols
                    // C and D are constant and not normalized by ATM vols

                    double atmVolRef = sqrt (atmVars[0] / tradYears[0]);
                    backboneSkewParams[0] = Maths::isZero(skewCelerity) ? 1.0
                        : 2.0 * vectorLinEq1[0] / (skewCelerity * atmVars[0]*atmVars[0])
                            * pow(atmVolRef, skewBackboneScaling);
                    backboneConvexityParams[0] = Maths::isZero(convexityCelerity) ? 1.0
                        : 3.0 * vectorLinEq2[0] /
                            (convexityCelerity*convexityCelerity * atmVars[0]*atmVars[0]*atmVars[0])
                            * pow(atmVolRef, convexityBackboneScaling);

                    backboneSkewCelerities[0] = skewCelerity;
                    backboneConvexityCelerities[0] = convexityCelerity;

                    for (int i = 1; i < tradYears.size(); i++) {
                        backboneSkewCelerities[i] = skewCelerity;
                        backboneConvexityCelerities[i] = convexityCelerity;
                        double thisAtmVar = atmVars[i];
                        double prevAtmVar = atmVars[i-1];
                        atmVolRef = sqrt (atmVars[i] / tradYears[i]);
                        backboneSkewParams[i] = Maths::isZero(skewCelerity) ? 1.0
                            : 2.0 * (vectorLinEq1[i] - vectorLinEq1[i-1])
                                / (skewCelerity * (thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar))
                                * pow(atmVolRef, skewBackboneScaling);
                        backboneConvexityParams[i] = Maths::isZero(convexityCelerity) ? 1.0
                            : 3.0 * (vectorLinEq2[i] - vectorLinEq2[i-1])
                                / (convexityCelerity*convexityCelerity *
                                (thisAtmVar*thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar*prevAtmVar))
                                * pow(atmVolRef, convexityBackboneScaling);
                    }
                }

            } catch(exception& e) {
                    ModelException(e, "VolHyperTrig::update");
            }
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(Helper, clazz);
            SUPERCLASS(VolHyperTrigPrimer);
            EMPTY_SHELL_METHOD(defaultHelper);
        }
        static IObject* defaultHelper(){
            return new Helper();
        }
    };
    typedef smartPtr<Helper> HelperSP;

    // helper class for original interpretation of input parameters
    class HelperOriginal: public VolHyperTrigPrimer{
    public:
        static CClassConstSP const TYPE;

        // constructor
        HelperOriginal(): VolHyperTrigPrimer(TYPE) {}
        HelperOriginal(const VolHyperTrig&  vht,
                       VolSurfaceConstSP    backbone): VolHyperTrigPrimer(TYPE,
                                                                          // registered fields
                                                                          vht.approxPade,
                                                                          vht.tenorRef,
                                                                          vht.strikeRefRatioBMs,
                                                                          vht.strikeRefRatios,
                                                                          vht.strikeRef,
                                                                          vht.skewFactorBMs,
                                                                          vht.skewFactors,
                                                                          vht.skew,
                                                                          vht.convexityFactorBMs,
                                                                          vht.convexityFactors,
                                                                          vht.convexity,
                                                                          vht.tenorTime,
                                                                          vht.useTimeDependentCelerities,
                                                                          vht.convexityPower,
                                                                          vht.skewBackboneScaling,
                                                                          vht.convexityBackboneScaling,
                                                                          // additional fields
                                                                          backbone,
                                                                          vht.gridType,
                                                                          vht.numGrid,
                                                                          vht.TruncationStd,
                                                                          vht.useGrid,
                                                                          vht.basis,
                                                                          vht.markedDate) {
            static const string method("VolHyperTrig::HelperOriginal::HelperOriginal");
            // step I:   general update already in constructor of VolHyperTrigPrimer via call to validatePop2Object
            // step II: specific update function of HelperOriginal class (diff for VHT and VHTDual)
            update(skew,
                   convexity,
                   skewFactors,
                   convexityFactors,
                   vht.skewCelerity,
                   vht.convexityCelerity);
        }

        void update(double             skewMod,
                    double             convexityMod,
                    DoubleArray        skewFactorsMod,
                    DoubleArray        convexityFactorsMod,
                    double             skewCelerity,
                    double             convexityCelerity) {
            static const string method("VolHyperTrig::HelperOriginal::update");
            static const double ln100To90 = -log(0.9);
            try{
                // step I:  general update -> delegated to VolHyperTrigPrimer (the same for VHT and VHTdual)
                skew = skewMod;
                convexity = convexityMod;
                skewFactors = skewFactorsMod;
                convexityFactors = convexityFactorsMod;
                updateParent();

                // step II: specific update for VHT

                // calculate skew and convexity celerities at each backbone date
                int iTradYear = 0;
                for (iTradYear = 0; iTradYear < tradYears.size(); ++iTradYear){
                    double thisTradYear = tradYears[iTradYear];
                    backboneSkewCelerities[iTradYear] =
                        skewCelerity * skewFactorInterp->value(thisTradYear)
                            / sqrt(thisTradYear / tradYearRef);
                    backboneConvexityCelerities[iTradYear] =
                        convexityCelerity * convexityFactorInterp->value(thisTradYear)
                            / sqrt(thisTradYear / tradYearRef);
                }
                // skew + convexity
                double skewParam = 1.0;
                double convexityParam = 1.0;
                double lastVar = 0.0;
                double dAtmVar_dx = 0.0, d2AtmVar_dx2 = 0.0;
                for (int iDate = 0; iDate < dates.size() && dates[iDate] <= dateRef; ++iDate){
                    double thisVar = atmVars[iDate];
                    double thisSqVar = thisVar * thisVar;
                    double lastSqVar = lastVar * lastVar;
                    dAtmVar_dx += 0.5 * (thisSqVar - lastSqVar) * backboneSkewCelerities[iDate];
                    d2AtmVar_dx2 += (thisSqVar * thisVar - lastSqVar * lastVar) / 3.0
                                    * Maths::square(backboneConvexityCelerities[iDate]);
                    lastVar = thisVar;
                }
                dAtmVar_dx *= skewParam / atmVarRef;
                d2AtmVar_dx2 *= convexityParam / (atmVarRef * atmVarRef);
                double dAtmVol_dx = skew / ln100To90;
                double d2AtmVol_dx2 = convexity / (ln100To90 * ln100To90);
                skewParam = 2.0 * atmVolRef * dAtmVol_dx * tradYearRef / dAtmVar_dx;
                convexityParam = 2.0 * (dAtmVol_dx * dAtmVol_dx + atmVolRef * d2AtmVol_dx2)
                                 * tradYearRef / d2AtmVar_dx2;
                for (iTradYear = 0; iTradYear < tradYears.size(); ++iTradYear){
                    backboneSkewParams[iTradYear] = skewParam;
                    backboneConvexityParams[iTradYear] = convexityParam;
                }

            } catch(exception& e) {
                    ModelException(e, "VolHyperTrig::update");
            }
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(HelperOriginal, clazz);
            SUPERCLASS(VolHyperTrigPrimer);
            EMPTY_SHELL_METHOD(defaultHelperOriginal);
        }
        static IObject* defaultHelperOriginal(){
            return new HelperOriginal();
        }
    };
    typedef smartPtr<HelperOriginal> HelperOriginalSP;

    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class VolHyperTrigVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        VolHyperTrigVolParam(): CVolParam(TYPE){}

        ~VolHyperTrigVolParam(){}

        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const{
            // turn the vol into what we must have
            const VolHyperTrig* myVol = static_cast<const VolHyperTrig *>(vol);
            // then just pass through the parameterised vol
            myVol->computeImpVol(strikes, maturities, impV);
        }

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface* spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const{
            // turn the vol into what we must have
            const VolHyperTrig* myVol = static_cast<const VolHyperTrig *>(vol);
            // then just pass through the parameterised vol
            return myVol->spotVolSurfaceFromStrikes(strikes);
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(VolHyperTrigVolParam, clazz);
            SUPERCLASS(CVolParam);
            EMPTY_SHELL_METHOD(defaultCtor);
        }
        static IObject* defaultCtor(){
            return new VolHyperTrigVolParam();
        }
    };
public:
    friend class VolHyperTrigVolParam;
    friend class HelperFull;
    friend class Helper;
    friend class HelperOriginal;
    friend class ParameterOutputAddin;

    static const string GRID_TYPE_NONE;
    static const string GRID_TYPE_LINEAR;
    static const string GRID_TYPE_SPLINE;
    static const string SPEEDUP_METHOD_GRID;

    static CClassConstSP const TYPE;

    ~VolHyperTrig(){}

    void validatePop2Object(){
        static const string method("VolHyperTrig::validatePop2Object");
        try{
            // the main part is in validatePop2Object of VolHyperTrigPrimer
            Calibrator::IAdjustable::checkRange(this);

            // Cache of expiry equivalents to the BMs
            // Here to facilitate the pointwise sens of skewFactors etc
            // Note the Primer class does some work with the "BMs" which
            // should be equivalent to this...
            int tenorTimeInt = DateTime::timeConvert(tenorTime);
            int i;
            skewFactorExpiries = ExpiryArraySP(new ExpiryArray(0));
            for (i=0; i<skewFactorBMs.size(); i++) {
                skewFactorExpiries->push_back(ExpirySP(new MaturityTimePeriod(skewFactorBMs[i], tenorTimeInt)));
            }
            convexityFactorExpiries = ExpiryArraySP(new ExpiryArray(0));
            for (i=0; i<convexityFactorBMs.size(); i++) {
                convexityFactorExpiries->push_back(ExpirySP(new MaturityTimePeriod(convexityFactorBMs[i], tenorTimeInt)));
            }
            skewCelerityFactorExpiries = ExpiryArraySP(new ExpiryArray(0));
            for (i=0; i<skewCelerityFactorBMs.size(); i++) {
                skewCelerityFactorExpiries->push_back(ExpirySP(new MaturityTimePeriod(skewCelerityFactorBMs[i], tenorTimeInt)));
            }
            convexityCelerityFactorExpiries = ExpiryArraySP(new ExpiryArray(0));
            for (i=0; i<convexityCelerityFactorBMs.size(); i++) {
                convexityCelerityFactorExpiries->push_back(ExpirySP(new MaturityTimePeriod(convexityCelerityFactorBMs[i], tenorTimeInt)));
            }

        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
        return new VolHyperTrigVolParam();
    }

    string sensName(DeltaSurface* shift) const{
        return getName();
    }

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const{
        return CVolBaseParamSurface::getName();
    }

    bool sensShift(DeltaSurface* shift){
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)){
            // only bother if non zero
            // what to do - spot moves to spot + shift
            // simply move strike ref k -> k (S+dS)/S
            double spot    = shift->getSpot();
            double newSpot = spot * (1.0 + shiftSize);
            strikeRef *= newSpot/spot;
            // get new backbone
            // Dummy try/catch to avoid vc6.opt crash
            try { buildCache(); } catch (...) { throw; }
        }
        // DON'T tweak the backbone !!!
        return false; // no more tweaking required here
    }

    // ----- VHTSkewFactorsPointwise -----
    string sensName(VHTSkewFactorsPointwise* shift) const{
        return getName();
    }
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    ExpiryArrayConstSP sensExpiries(VHTSkewFactorsPointwise* shift) const{
        return skewFactorExpiries;
    }
    bool sensShift(VHTSkewFactorsPointwise* shift){
        double shiftSize = shift->getShiftSize();
        // only bother if non zero
        if (!Maths::isZero(shiftSize)){
            int expiryIdx = shift->getExpiry()->search(skewFactorExpiries.get());
            skewFactors[expiryIdx] += shiftSize;
            // get new backbone
            // Dummy try/catch to avoid vc6.opt crash
            try { buildCache(); } catch (...) { throw; }
        }
        // DON'T tweak the backbone !!!
        return false; // no more tweaking required here
    }

    // ----- VHTSkewCelerityFactorsPointwise -----
    string sensName(VHTSkewCelerityFactorsPointwise* shift) const{
        // This sens is only meaningful for later versions of VHT
        // By returning an empty name we stop the sens for appropriate cases
        if (version<2 || skewCelerityFactorBMs.size()<1) {
            return "";
        }
        return getName();
    }
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    ExpiryArrayConstSP sensExpiries(VHTSkewCelerityFactorsPointwise* shift) const{
        return skewCelerityFactorExpiries;
    }
    bool sensShift(VHTSkewCelerityFactorsPointwise* shift){
        double shiftSize = shift->getShiftSize();
        // only bother if non zero
        if (!Maths::isZero(shiftSize)){
            int expiryIdx = shift->getExpiry()->search(skewCelerityFactorExpiries.get());
            skewCelerityFactors[expiryIdx] += shiftSize;
            // get new backbone
            // Dummy try/catch to avoid vc6.opt crash
            try { buildCache(); } catch (...) { throw;  }
        }
        // DON'T tweak the backbone !!!
        return false; // no more tweaking required here
    }

    // ----- VHTConvexityFactorsPointwise -----
    string sensName(VHTConvexityFactorsPointwise* shift) const{
        return getName();
    }
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    ExpiryArrayConstSP sensExpiries(VHTConvexityFactorsPointwise* shift) const{
        return convexityFactorExpiries;
    }
    bool sensShift(VHTConvexityFactorsPointwise* shift){
        double shiftSize = shift->getShiftSize();
        // only bother if non zero
        if (!Maths::isZero(shiftSize)){
            int expiryIdx = shift->getExpiry()->search(convexityFactorExpiries.get());
            convexityFactors[expiryIdx] += shiftSize;
            // get new backbone
            // Dummy try/catch to avoid vc6.opt crash
            try { buildCache(); } catch (...) { throw; }
        }
        // DON'T tweak the backbone !!!
        return false; // no more tweaking required here
    }

    // ----- VHTConvexityCelerityFactorsPointwise -----
    string sensName(VHTConvexityCelerityFactorsPointwise* shift) const{
        // This sens is only meaningful for later versions of VHT
        // By returning an empty name we stop the sens for appropriate cases
        if (version<2 || convexityCelerityFactorBMs.size()<1) {
            return "";
        }
        return getName();
    }
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    ExpiryArrayConstSP sensExpiries(VHTConvexityCelerityFactorsPointwise* shift) const{
        return convexityCelerityFactorExpiries;
    }
    bool sensShift(VHTConvexityCelerityFactorsPointwise* shift){
        double shiftSize = shift->getShiftSize();
        // only bother if non zero
        if (!Maths::isZero(shiftSize)){
            int expiryIdx = shift->getExpiry()->search(convexityCelerityFactorExpiries.get());
            convexityCelerityFactors[expiryIdx] += shiftSize;
            // get new backbone
            // Dummy try/catch to avoid vc6.opt crash
            try { buildCache(); } catch (...) { throw; }
        }
        // DON'T tweak the backbone !!!
        return false; // no more tweaking required here
    }

	//------ VegaAtmParallelConstConvx ----- 
	/** Returns name identifying VolHyperTrig for VEGA_PARALLEL_CONST_CONVX */
	string sensName(const VolAtmParallelConstConvx*) const{
		return getName();
	}

	/** Shifts the object using given shift */
	TweakOutcome sensShift(const PropertyTweak<VolAtmParallelConstConvx>& tweak){
	static const string method = "VolHyperTrig::sensShift";
	try {
		if (!Maths::isZero(tweak.coefficient)){
			// tweak convexity and convexityFactor
			double atmVolRef = calcVol(strikeRef, vhtHelper->tradYearRef);
			convexity *= Maths::square(atmVolRef/(atmVolRef + tweak.coefficient));

			int iFactor = 0;
			int tenorTimeInt = DateTime::timeConvert(tenorTime);
			for (; iFactor < convexityFactors.size(); ++iFactor){
				MaturityPeriod	period(convexityFactorBMs[iFactor]);
				DateTime		date(period.toDate(vhtHelper->baseDate).getDate(),tenorTimeInt);
				double			tradYear = vhtHelper->timeMetric->yearFrac(vhtHelper->baseDate, date);
				double			strike = strikeRef;
				if (vhtHelper->strikeRefRatioInterp.get()){
					strike *= vhtHelper->strikeRefRatioInterp->value(tradYear);
				} 
				double atmVol = calcVol(strikeRef, tradYear);
				convexityFactors[iFactor] *= Maths::square(atmVol*(atmVolRef + tweak.coefficient)/(atmVol + tweak.coefficient)/atmVolRef);
			}

			// tweak backbone
			VolSurfaceSP backbone(const_cast<VolSurface*>(getBackboneSurface()));

			backbone->sensShift(PropertyTweak<VolParallel>(tweak.coefficient,
                                                           tweak.marketDataName));

			// update others ...
			buildCache(backbone);
		}
		return TweakOutcome(tweak.coefficient, false); // none of our components has a rho type sensitivity
	}
	catch (exception &e) {
		throw ModelException(e, method,
							"VolParallel tweak failed for " + getName());
		}
	}

	//------ VegaAtmPointwiseConstConvx ----- 
	/** Returns name identifying vol for vega pointwise */
	string sensName(const VolAtmPointwiseConstConvx*) const{
		return getName();
	}

	/** Returns the array of expiries (ie maturities/benchmark dates) that
		need to be tweaked for this vol */
	ExpiryWindowArrayConstSP sensQualifiers(const VolAtmPointwiseConstConvx*) const{
		return ExpiryWindow::series(getBackboneSurface()->getExpiries());
	}

	/** Shifts the object using given shift */
	TweakOutcome sensShift(const PropertyTweak<VolAtmPointwiseConstConvx>& tweak){
		static const string routine = "VolHyperTrig::sensShift<VolAtmPointwiseConstConvx>";
		try{
			if (!Maths::isZero(tweak.coefficient)){
				// tweak convexity or convexityFactor
				// check whether the maturity is in convexityFactorBMs. if not, insert it.
				int			tenorTimeInt = DateTime::timeConvert(tenorTime);
				DateTime	dateFrom(vhtHelper->baseDate.getDate(),tenorTimeInt);
				DateTime	dateTo(tweak.qualifier->expiry->toDate(vhtHelper->baseDate).getDate(), tenorTimeInt);
				MaturityPeriodSP period(MaturityPeriod::dateSubtract(dateTo, dateFrom));
				string		bm = period->toString();
				double		tradYear = vhtHelper->timeMetric->yearFrac(vhtHelper->baseDate, dateTo);
				double		strike = strikeRef;
				if (vhtHelper->strikeRefRatioInterp.get()){
					strike *= vhtHelper->strikeRefRatioInterp->value(tradYear);
				} 
				double		atmVol = calcVol(strike, tradYear);

				if(bm == tenorRef){
					convexity *= Maths::square(atmVol/(atmVol + tweak.coefficient));
				} else {
					int			iFactor = -1;
					DateTime	iDate; 
					do{
						MaturityPeriod	iPeriod(convexityFactorBMs[++iFactor]);
						iDate = DateTime(iPeriod.toDate(vhtHelper->baseDate).getDate(),tenorTimeInt);
					} while ( iDate < dateTo && (iFactor < convexityFactors.size()-1));

					if (iDate == dateTo){
						convexityFactors[iFactor] *= Maths::square(atmVol/(atmVol + tweak.coefficient));
					} else {
						StringArray tmpBMs(convexityFactorBMs);
						DoubleArray tmpFactors(convexityFactors);
						int size = convexityFactorBMs.size();
						convexityFactorBMs.resize(size + 1);
						convexityFactors.resize(size + 1);

						if(iDate < dateTo) ++iFactor;

						int iNewFactor = 0;
						while (iNewFactor < iFactor){
							convexityFactorBMs[iNewFactor] = tmpBMs[iNewFactor];
							convexityFactors[iNewFactor] = tmpFactors[iNewFactor];
							++iNewFactor;
						}
						convexityFactorBMs[iNewFactor] = bm;
						convexityFactors[iNewFactor] = vhtHelper->convexityFactorInterp->value(tradYear)
							* Maths::square(atmVol/(atmVol + tweak.coefficient));
						++iNewFactor;
						while (iNewFactor < convexityFactorBMs.size()){
							convexityFactorBMs[iNewFactor] = tmpBMs[iNewFactor-1];
							convexityFactors[iNewFactor] = tmpFactors[iNewFactor-1];
							++iNewFactor;
						}
					}
				}

				// tweak backbone
				VolSurfaceSP backbone(const_cast<VolSurface*>(getBackboneSurface()));

				backbone->sensShift(PropertyTweak<VolPointwise>(tweak.coefficient,
                                                                tweak.marketDataName,
                                                                tweak.qualifier));

				// update others ...
				buildCache(backbone);
			}
			return TweakOutcome(tweak.coefficient, false); 
		} catch (exception& e){
			throw ModelException(e, routine);
		}
	}

    IVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                   const CAsset*      asset) const {

        static const string routine = "VolHyperTrig::getProcessedVol";

        try{
            // Determine if it is a VarSwapBasis request
            const VarSwapBasis::VarSwapBasisRequest* basisRequest = dynamic_cast<const VarSwapBasis::VarSwapBasisRequest*>(volRequest);
            if(basisRequest) {
                if(basis.get()){
                    VolSurfaceConstSP   backbone(getBackboneSurface());
                    try {
                        return vhtHelper->setup(asset, backbone.get(), basisRequest, tenorTime);
                    } catch(exception& e) {
                        // New communication language - better than QLIB exception mechanism
                        return new VarSwapBasis::VarSwapBasisProcError(ModelException(e, routine));
                    }
                } else {
                    return 0;
                }
            }

            // determine useGrid here
            const LocVolRequest*  request = dynamic_cast<const LocVolRequest*>(volRequest);
            if (request) {
                bool gridAllowed = !CString::equalsIgnoreCase(gridType, GRID_TYPE_NONE);
                bool gridRequested = CString::equalsIgnoreCase(request->getSpeedup(), SPEEDUP_METHOD_GRID);
                bool useGrid = gridAllowed && gridRequested;
                VolHyperTrig*   vht = const_cast<VolHyperTrig*>(this);
                vht->useGrid = useGrid;
                (vht->vhtHelper.get())->useGrid = useGrid;
                if (useGrid)
                    vht->vhtHelper.get()->volGrid.get()->setupGrid(vht->vhtHelper.get());
			}else{
				bool useGrid = false;
                VolHyperTrig*   vht = const_cast<VolHyperTrig*>(this);
                vht->useGrid = useGrid;
                (vht->vhtHelper.get())->useGrid = useGrid;
			}
            return VolBaseParam::getProcessedVol(volRequest, asset);
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    void getMarket(const IModel* model,
                   const MarketData* market){


        // Fetch the VarSwapBasis
        //MarketDataFetcherSP fetcher(new MarketDataFetcher());
        //fetcher->setVarSwapBasis(true);

        basis = VarSwapBasis::fetchBasis(model, market, this->getName());
        if(basis.get()) {
            basis->getMarket(model, market);
        }

        /** Retrieves time metric from parent's vol surface */
        VolBaseParam::getMarket(model, market);
    }

    StringArray getSkewFactorBMs() const {
        return skewFactorBMs;
    }

    StringArray getConvexityFactorBMs() const {
        return convexityFactorBMs;
    }

private:

    void computeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const {
        static const string method("VolHyperTrig::computeImpVol");
        if ((maturities.size() != strikes.size()) ||
            (maturities.size() != impV.size())) {
            throw ModelException(method, "Size mismatch between strikes ("+
                                 Format::toString(strikes.size()) +
                                 "), maturities ("+
                                 Format::toString(maturities.size())+
                                 ") and impV ("+
                                 Format::toString(impV.size())+ ")");
        }
        for (int iMat = 0; iMat < maturities.size(); iMat++) {
            int nbStrikes = strikes[iMat].size();
            if (nbStrikes != impV[iMat].size()){
                throw ModelException(method, "Size mismatch between strikes"
                                     " & maturities for Mat " +
                                     maturities[iMat].toString() +
                                     " (n "+ Format::toString(iMat) + ")");
            }

            DateTimeArray dates = vhtHelper->dates;
            double tradYear = vhtHelper->timeMetric->yearFrac(vhtHelper->baseDate, maturities[iMat]);

            for (int iStrike = 0; iStrike < nbStrikes; iStrike++) {
                impV[iMat][iStrike] = calcVol(strikes[iMat][iStrike],
                                              tradYear);
            }
        }
    }

    // calculate strike ref
    double getStrikeRef(double tradYear) const {
        return vhtHelper->getStrikeRef(tradYear);
    }

    // given a strike, a trading time, calculate implied volatility
    double calcVol(double strike,
                   double tradYear) const{
        return vhtHelper->calcVol(strike,tradYear);
    }

    /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
    VolSurface* spotVolSurfaceFromStrikes(
        const CDoubleArray&   strikes) const{
        static const string method("VolHyperTrig::spotVolSurfaceFromStrikes");
        try{
            const VolSurface* backbone = getBackboneSurface();
            CDoubleMatrix matrix = vhtHelper->spotVolMatrixFromStrikes(strikes);
            /** for performance use special constructors */
            return (new VolSurface(backbone, strikes, matrix));
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }
    // reference to VolHyperTrigPrimer class and thus to helper class
    VolHyperTrigPrimerSP    vhtHelper;

    // REGISTERED FIELDS
    string            approxPade;           // optional: 'exact' if closed-form, nber of Pade coeffs if approx

    DateTime          markedDate;           // The last marked date
    double            strikeRef;
    StringArray       strikeRefRatioBMs;    // polymorphic types not supported by IMS
    DoubleArray       strikeRefRatios;

    // ref point specifics
    string            tenorRef;             // polymorphic types not supported by IMS
    string            tenorTime;            // optional; applies to factorSkewBMs too if any

    double            skew;
    double            convexity;
    double            skewCelerity;
    double            convexityCelerity;

    // term structure
    StringArray       skewFactorBMs;                // polymorphic types not supported by IMS
    DoubleArray       skewFactors;
    StringArray       convexityFactorBMs;           // polymorphic types not supported by IMS
    DoubleArray       convexityFactors;
    StringArray       skewCelerityFactorBMs;        // polymorphic types not supported by IMS
    DoubleArray       skewCelerityFactors;
    StringArray       convexityCelerityFactorBMs;   // polymorphic types not supported by IMS
    DoubleArray       convexityCelerityFactors;

    bool              useTimeDependentCelerities;   // true, if C and D are time dependent
                                                    // false, if A and B are time dependent

    int               version;              // which version of vht to choose
                                            // version 0 = original interpretation of skewFactors and convexityFactors
                                            // version 1 = correct interpretation of skewFactors and convexityFactors

    double            convexityPower;

    double            skewBackboneScaling;
    double            convexityBackboneScaling;

    bool              done;                 // needed for creation of xml files when mapping input parameters of
                                            // version 0 to input parameters of version 1

    // vo grid
    string            gridType;
    int               numGrid;
    double            TruncationStd;

    // Transient
    bool              useGrid;              // determined by request from pricing engine and gridType together
    VarSwapBasisSP    basis;
    ExpiryArraySP     skewFactorExpiries;
    ExpiryArraySP     convexityFactorExpiries;
    ExpiryArraySP     skewCelerityFactorExpiries;
    ExpiryArraySP     convexityCelerityFactorExpiries;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolHyperTrig, clazz);
        SUPERCLASS(VolBaseParam);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(VHTSkewFactorsPointwise::IShift);
        IMPLEMENTS(VHTSkewCelerityFactorsPointwise::IShift);
        IMPLEMENTS(VHTConvexityFactorsPointwise::IShift);
        IMPLEMENTS(VHTConvexityCelerityFactorsPointwise::IShift);
        IMPLEMENTS(Calibrator::IAdjustable);
		IMPLEMENTS(ITweakableWithRespectTo<VolAtmParallelConstConvx>);
		IMPLEMENTS(ITweakableWithRespectTo<VolAtmPointwiseConstConvx>);
        EMPTY_SHELL_METHOD(defaultCtor);

        FIELD(approxPade, "Closed Form or Nber Pade Coefficients");
        FIELD_MAKE_OPTIONAL(approxPade);
        FIELD(markedDate, "Last marked date");
        FIELD_MAKE_OPTIONAL(markedDate);
        FIELD(strikeRef, "Reference strike");
        FIELD(strikeRefRatioBMs, "Reference strike ratio bencharks");
        FIELD_MAKE_OPTIONAL(strikeRefRatioBMs);
        FIELD(strikeRefRatios, "Reference strike ratios");
        FIELD_MAKE_OPTIONAL(strikeRefRatios);
        FIELD(tenorRef, "Reference tenor");
        FIELD(skew, "Reference skew");
        FIELD(convexity, "Reference convexity");
        FIELD(skewCelerity, "Skew Celerity");
        FIELD(convexityCelerity, "Convexity Celerity");
        FIELD(skewFactorBMs, "Skew Factor Bencharks");
        FIELD(skewFactors, "Skew Factors");
        FIELD(convexityFactorBMs, "Convexity Factor Bencharks");
        FIELD(convexityFactors, "Convexity Factors");

        FIELD(skewCelerityFactorBMs, "Skew Celerity Factor Bencharks");
        FIELD_MAKE_OPTIONAL(skewCelerityFactorBMs);
        FIELD(skewCelerityFactors, "Skew Celerity Factors");
        FIELD_MAKE_OPTIONAL(skewCelerityFactors);
        FIELD(convexityCelerityFactorBMs, "Convexity Celerity Factor Bencharks");
        FIELD_MAKE_OPTIONAL(convexityCelerityFactorBMs);
        FIELD(convexityCelerityFactors, "Convexity Celerity Factors");
        FIELD_MAKE_OPTIONAL(convexityCelerityFactors);

        FIELD(tenorTime, "Time of Day for Tenors");
        FIELD_MAKE_OPTIONAL(tenorTime);
        FIELD(useTimeDependentCelerities, "use time dependent celerities");
        FIELD_MAKE_OPTIONAL(useTimeDependentCelerities);
        FIELD(version, "which version of vht to choose");
        FIELD_MAKE_OPTIONAL(version);
        FIELD(convexityPower, "Propagation power for convexity factors");
        FIELD_MAKE_OPTIONAL(convexityPower);
        FIELD(skewBackboneScaling, "Scaling power for skew parameter or skew celerity");
        FIELD_MAKE_OPTIONAL(skewBackboneScaling);
        FIELD(convexityBackboneScaling, "Scaling power for convexity parameter or convexity celerity");
        FIELD_MAKE_OPTIONAL(convexityBackboneScaling);
        FIELD(vhtHelper, "");
        FIELD_MAKE_TRANSIENT(vhtHelper);
        FIELD(done, "");
        FIELD_MAKE_TRANSIENT(done);
        FIELD(gridType, "type of vol interpolation: NONE, LINEAR or SPLINE");
        FIELD_MAKE_OPTIONAL(gridType);
        FIELD(numGrid, "number of grid points");
        FIELD_MAKE_OPTIONAL(numGrid);
        FIELD(TruncationStd, "trancation of standard normal variable");
        FIELD_MAKE_OPTIONAL(TruncationStd);
        FIELD(useGrid, "");
        FIELD_MAKE_TRANSIENT(useGrid);
        FIELD(basis, "");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(basis);
        FIELD(skewFactorExpiries, "");
        FIELD_MAKE_TRANSIENT(skewFactorExpiries);
        FIELD(convexityFactorExpiries, "");
        FIELD_MAKE_TRANSIENT(convexityFactorExpiries);
        FIELD(skewCelerityFactorExpiries, "");
        FIELD_MAKE_TRANSIENT(skewCelerityFactorExpiries);
        FIELD(convexityCelerityFactorExpiries, "");
        FIELD_MAKE_TRANSIENT(convexityCelerityFactorExpiries);

        // add our fields and their ranges to central list
        Calibrator::IAdjustable::registerField(
            clazz, "skew",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "convexity",
            new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "skewCelerity",
            new Range(OpenBoundary(0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "convexityCelerity",
            new Range(OpenBoundary(0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "skewFactors",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "convexityFactors",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "skewCelerityFactors",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "convexityCelerityFactors",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
    }

    VolHyperTrig():
    VolBaseParam(TYPE),
    approxPade("EXACT"),
    strikeRef(0.0),
    tenorRef("1Y"),
    tenorTime("EOD"),
    skew(-0.02),
    convexity(0.002),
    skewCelerity(3.0),
    convexityCelerity(0.3),
    useTimeDependentCelerities(true),
    version(0),
    convexityPower(0.5),
    skewBackboneScaling(0.0),
    convexityBackboneScaling(0.0),
    done(false),
    gridType(GRID_TYPE_NONE),
    numGrid(1000),
    TruncationStd(5.0),
    useGrid(false)
        {}

    static IObject* defaultCtor(){
        return new VolHyperTrig();
    }

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache() {
        static const string method("VolHyperTrig::buildCache");
        try{
            // in order to create modified input files for cases (parameter mapping)
            if (done) {
                return;
            }

            // step I: information in backbone has to be passed to vhtHelper
            VolSurfaceConstSP backbone(getBackboneSurface());
            // step II: create vhtHelper and pass info of vht and backbone
            if (version==0) {
                // Dummy try/catch to avoid vc6.opt crash
				try {
                    vhtHelper = HelperOriginalSP(new HelperOriginal(*this, backbone));
                } catch (...) { throw; }
            } else if (version==1) {
                vhtHelper = HelperSP(new Helper(*this, backbone));
            } else if (version==2) {
                vhtHelper = HelperFullSP(new HelperFull(*this, backbone));
            } else {
                throw ModelException("no correct version for class VolHyperTrig");
            }            

#ifdef CREATE_MODIFIED_XML_FILES
            // Override initial VHT fields with implied values
            // modify: skewFactorBMs, skewFactors, convexityFactorBMs, convexityFactors
            int nbDates = vhtHelper->backboneSkewCelerities.size();
            skewFactors.resize(nbDates);
            convexityFactors.resize(nbDates);
            skewFactorBMs.resize(nbDates);
            convexityFactorBMs.resize(nbDates);

            DoubleArray skewTermStructure(nbDates);
            DoubleArray convexityTermStructure(nbDates);

            double skewHelper = 0.0;
            double convexityHelper = 0.0;
            double skewParam = vhtHelper->backboneSkewParams[vhtHelper->dateRefIdx];
            double convexityParam = vhtHelper->backboneConvexityParams[vhtHelper->dateRefIdx];
            int i = 0;
            for (i = 0; i < nbDates; i++) {
                skewHelper += (i==0) ? vhtHelper->atmVars[0]*vhtHelper->atmVars[0] * vhtHelper->backboneSkewCelerities[0]
                    : (vhtHelper->atmVars[i]*vhtHelper->atmVars[i] - vhtHelper->atmVars[i-1]*vhtHelper->atmVars[i-1]) * vhtHelper->backboneSkewCelerities[i];
                skewTermStructure[i] = skewHelper /
                    (4.0 * vhtHelper->tradYears[i] * vhtHelper->atmVars[i] * sqrt(vhtHelper->atmVars[i]/vhtHelper->tradYears[i]) / skewParam);

                convexityHelper += (i==0) ? vhtHelper->atmVars[0]*vhtHelper->atmVars[0]*vhtHelper->atmVars[0]
                                                * vhtHelper->backboneConvexityCelerities[0]*vhtHelper->backboneConvexityCelerities[0]
                    : (vhtHelper->atmVars[i]*vhtHelper->atmVars[i]*vhtHelper->atmVars[i] - vhtHelper->atmVars[i-1]*vhtHelper->atmVars[i-1]*vhtHelper->atmVars[i-1])
                        * (vhtHelper->backboneConvexityCelerities[i]*vhtHelper->backboneConvexityCelerities[i]);
                double convexityTermStructureTemp = convexityHelper /
                    (6.0 * vhtHelper->tradYears[i] * vhtHelper->atmVars[i]*vhtHelper->atmVars[i] / convexityParam);
                convexityTermStructure[i] =
                    (convexityTermStructureTemp - skewTermStructure[i]*skewTermStructure[i])
                        / sqrt(vhtHelper->atmVars[i]/vhtHelper->tradYears[i]);
            }
            for (i = 0; i < nbDates; i++) {
                skewFactors[i] = skewTermStructure[i] / skewTermStructure[vhtHelper->dateRefIdx] *
                    sqrt(vhtHelper->tradYears[i] / vhtHelper->tradYears[vhtHelper->dateRefIdx]);
                convexityFactors[i] = convexityTermStructure[i] / convexityTermStructure[vhtHelper->dateRefIdx] *
                    sqrt(vhtHelper->tradYears[i] / vhtHelper->tradYears[vhtHelper->dateRefIdx]);
                MaturityPeriodSP period(MaturityPeriod::dateSubtract(vhtHelper->dates[i], vhtHelper->baseDate));
                string periodString = period->toString();
                skewFactorBMs[i] = periodString;
                convexityFactorBMs[i] = periodString;
            }
            done = true;
#endif
        }
        catch(exception& e){
          throw ModelException(e, method);
        }
    }

	/*** Build the parameterised vol and cache any values **/
    void buildCache(VolSurfaceSP volSurfaceForBackbone) {
        static const string method("VolHyperTrig::buildCache");
        try{
            // create vhtHelper and pass info of vht and backbone
            if (version==0) {
                // Dummy try/catch to avoid vc6.opt crash
				try {
                    vhtHelper = HelperOriginalSP(new HelperOriginal(*this, volSurfaceForBackbone));
                } catch (...) { throw; }
            } else if (version==1) {
                vhtHelper = HelperSP(new Helper(*this, volSurfaceForBackbone));
            } else if (version==2) {
                vhtHelper = HelperFullSP(new HelperFull(*this, volSurfaceForBackbone));
            } else {
                throw ModelException("no correct version for class VolHyperTrig");
            }            
        }
        catch(exception& e){
          throw ModelException(e, method);
        }
    }

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields){
        update();
    }

    /** Called after (calibrator) adjustments have been made to fields */
    void update(){
        if(version==2) {
            vhtHelper->update(skew,
                              convexity,
                              skewCelerity,
                              convexityCelerity,
                              skewFactors,
                              convexityFactors,
                              skewCelerityFactors,
                              convexityCelerityFactors);
        } else {
            vhtHelper->update(skew,
                              convexity,
                              skewFactors,
                              convexityFactors,
                              skewCelerity,
                              convexityCelerity);
        }
    }

public:
#ifdef CREATE_MODIFIED_XML_FILES
    void write(const string& tag, Writer* writer) const {
        CVolBaseParamSurface::write(tag, writer);
    }

    IObject* clone() const {
        VolHyperTrig* ptr = const_cast<VolHyperTrig*>(this);
        return ptr;
    }
#endif

};

CClassConstSP const VolHyperTrig::TYPE =
CClass::registerClassLoadMethod("VolHyperTrig", typeid(VolHyperTrig), load);

CClassConstSP const VolHyperTrig::VolHyperTrigVolParam::TYPE =
CClass::registerClassLoadMethod("VolHyperTrig::VolHyperTrigVolParam", typeid(VolHyperTrigVolParam), load);

CClassConstSP const VolHyperTrig::HelperFull::TYPE =
CClass::registerClassLoadMethod("HelperFull", typeid(VolHyperTrig::HelperFull), VolHyperTrig::HelperFull::load);

CClassConstSP const VolHyperTrig::Helper::TYPE =
CClass::registerClassLoadMethod("Helper", typeid(VolHyperTrig::Helper), VolHyperTrig::Helper::load);

CClassConstSP const VolHyperTrig::HelperOriginal::TYPE =
CClass::registerClassLoadMethod("HelperOriginal", typeid(VolHyperTrig::HelperOriginal), VolHyperTrig::HelperOriginal::load);

// external symbol to allow class to be forced to be linked in
bool VolHyperTrigLoad(){
    return (VolHyperTrig::TYPE != 0);
}

typedef smartPtr<VolHyperTrig> VolHyperTrigSP;

const string VolHyperTrig::GRID_TYPE_NONE      =   "none";
const string VolHyperTrig::GRID_TYPE_LINEAR    =   "linear";
const string VolHyperTrig::GRID_TYPE_SPLINE    =   "spline";
const string VolHyperTrig::SPEEDUP_METHOD_GRID =   "grid";

/** skewParam parameterised CVolBase that supports BS and DVF views of the world */
class VolHyperTrigDual: public VolBaseParam,
                        virtual public IVolProcessed,
                        virtual public IVolatilityBS,
                        virtual public IVolatilityDVF,
                        virtual public DeltaSurface::IShift,
                        public virtual Calibrator::IAdjustable{
private:
    class HelperDual: public VolHyperTrigPrimer{
    public:
        static CClassConstSP const TYPE;

        // constructor
        HelperDual(): VolHyperTrigPrimer(TYPE) {}
        HelperDual(const VolHyperTrigDual&  vhtDual,
                   VolSurfaceConstSP        backbone): VolHyperTrigPrimer(TYPE,
                                                                          // registered fields
                                                                          vhtDual.approxPade,
                                                                          vhtDual.tenorRef,
                                                                          vhtDual.strikeRefRatioBMs,
                                                                          vhtDual.strikeRefRatios,
                                                                          vhtDual.strikeRef,
                                                                          vhtDual.skewFactorBMs,
                                                                          vhtDual.skewFactors,
                                                                          vhtDual.skew,
                                                                          vhtDual.convexityFactorBMs,
                                                                          vhtDual.convexityFactors,
                                                                          vhtDual.convexity,
                                                                          vhtDual.tenorTime,
                                                                          vhtDual.useTimeDependentCelerities,
                                                                          vhtDual.convexityPower,
                                                                          vhtDual.skewBackboneScaling,
                                                                          vhtDual.convexityBackboneScaling,
                                                                          // additional fields
                                                                          backbone,
                                                                          vhtDual.gridType,
                                                                          vhtDual.numGrid,
                                                                          vhtDual.TruncationStd,
                                                                          vhtDual.useGrid,
                                                                          vhtDual.basis,
                                                                          vhtDual.markedDate) {
            static const string method("VolHyperTrigDual::HelperDual::HelperDual");
            // step I:   general update already in constructor of VolHyperTrigPrimer via call to validatePop2Object
            // step II: specific update function of HelperDual class (diff for VHT and VHTDual)
            update(skew,
                   convexity,
                   skewFactors,
                   convexityFactors,
                   vhtDual.levelPlusInf,
                   vhtDual.levelMinusInf);
        }

        void update(double             skewMod,
                    double             convexityMod,
                    DoubleArray        skewFactorsMod,
                    DoubleArray        convexityFactorsMod,
                    double             levelPlusInf,
                    double             levelMinusInf) {
            static const string method("VolHyperTrigDual::HelperDual::update");
            try{
                // step I:  general update -> delegated to VolHyperTrigDualPrimer (the same for VHT and VHTdual)
                skew = skewMod;
                convexity = convexityMod;
                skewFactors = skewFactorsMod;
                convexityFactors = convexityFactorsMod;
                updateParent();

                // step II: specific update for VHTDual
                double skewParam = (levelPlusInf - levelMinusInf) / 2.0;
                double convexityParam = (levelPlusInf + levelMinusInf - 2.0) / 2.0;
                // floor to zero so that calibrator does not fail
                if(Maths::isNegative(convexityParam)) {
                    convexityParam = 0.0;
                }
                if (useTimeDependentCelerities) {
                    backboneSkewCelerities[0] = Maths::isZero(skewParam) ? 1.0
                        : 2.0 * vectorLinEq1[0] / (skewParam * atmVars[0]*atmVars[0]);
                    backboneConvexityCelerities[0] = Maths::isZero(convexityParam) ? 1.0
                        : sqrt(3.0 * vectorLinEq2[0] / (convexityParam * atmVars[0]*atmVars[0]*atmVars[0]));
                    backboneSkewParams[0] = skewParam;
                    backboneConvexityParams[0] = convexityParam;
                    for (int i = 1; i < tradYears.size(); i++) {
                        backboneSkewParams[i] = skewParam;
                        backboneConvexityParams[i] = convexityParam;
                        double thisAtmVar = atmVars[i];
                        double prevAtmVar = atmVars[i-1];
                        backboneSkewCelerities[i] = Maths::isZero(skewParam) ? 1.0
                            : 2.0 * (vectorLinEq1[i] - vectorLinEq1[i-1])
                                / (skewParam * (thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar));
                        if (Maths::isNegative(vectorLinEq2[i] - vectorLinEq2[i-1])) {
                            backboneConvexityCelerities[i] = 0.0;
                        } else {
                            backboneConvexityCelerities[i] = Maths::isZero(convexityParam) ? 1.0
                                : sqrt(3.0 * (vectorLinEq2[i] - vectorLinEq2[i-1])
                                    / (convexityParam *
                                    (thisAtmVar*thisAtmVar*thisAtmVar-prevAtmVar*prevAtmVar*prevAtmVar)));
                        }
                    }
                } else {
                    throw ModelException("useTimeDependentCelerties = false not possible in VHTDual");
                }
            } catch(exception& e) {
                    ModelException(e, "VolHyperTrigDual::update");
            }
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(HelperDual, clazz);
            SUPERCLASS(VolHyperTrigPrimer);
            EMPTY_SHELL_METHOD(defaultHelperDual);
        }
        static IObject* defaultHelperDual(){
            return new HelperDual();
        }
    };
    typedef smartPtr<HelperDual> HelperDualSP;
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class VolHyperTrigDualVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        VolHyperTrigDualVolParam(): CVolParam(TYPE){}

        ~VolHyperTrigDualVolParam(){}

        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const{
            // turn the vol into what we must have
            const VolHyperTrigDual* myVol = static_cast<const VolHyperTrigDual *>(vol);
            // then just pass through the parameterised vol
            myVol->computeImpVol(strikes, maturities, impV);
        }

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface* spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const{
            // turn the vol into what we must have
            const VolHyperTrigDual* myVol = static_cast<const VolHyperTrigDual *>(vol);
            // then just pass through the parameterised vol
            return myVol->spotVolSurfaceFromStrikes(strikes);
        }

    private:
        static void load(CClassSP& clazz){
            REGISTER(VolHyperTrigDualVolParam, clazz);
            SUPERCLASS(CVolParam);
            EMPTY_SHELL_METHOD(defaultCtor);
        }
        static IObject* defaultCtor(){
            return new VolHyperTrigDualVolParam();
        }
    };
public:
    friend class VolHyperTrigDualVolParam;
    friend class HelperDual;
    friend class ParameterOutputAddin;
    static CClassConstSP const TYPE;

    ~VolHyperTrigDual(){}

    void validatePop2Object(){
        static const string method("VolHyperTrigDual::validatePop2Object");
        try{
            // the main part is in validatePop2Object of VolHyperTrigPrimer
            Calibrator::IAdjustable::checkRange(this);
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
        return new VolHyperTrigDualVolParam();
    }

    string sensName(DeltaSurface* shift) const{
        return getName();
    }

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const{
        return CVolBaseParamSurface::getName();
    }

    bool sensShift(DeltaSurface* shift){
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)){
            // only bother if non zero
            // what to do - spot moves to spot + shift
            // simply move strike ref k -> k (S+dS)/S
            double spot    = shift->getSpot();
            double newSpot = spot * (1.0 + shiftSize);
            strikeRef *= newSpot/spot;
            // get new backbone
            // Dummy try/catch to avoid vc6.opt crash
            try { buildCache(); } catch (...) { throw; }
        }
        // DON'T tweak the backbone !!!
        return false; // no more tweaking required here
    }

    IVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                   const CAsset*      asset) const {
        // determine useGrid here
        const LocVolRequest*  request = dynamic_cast<const LocVolRequest*>(volRequest);
        if (request) {
            bool gridAllowed = !CString::equalsIgnoreCase(gridType, VolHyperTrigPrimer::VolGrid::GRID_TYPE_NONE);
            bool gridRequested = CString::equalsIgnoreCase(request->getSpeedup(), VolHyperTrigPrimer::VolGrid::SPEEDUP_METHOD_GRID);
            bool useGrid = gridAllowed && gridRequested;
            VolHyperTrigDual*   vhtDual = const_cast<VolHyperTrigDual*>(this);
            vhtDual->useGrid = useGrid;
            (vhtDual->vhtDualHelper.get())->useGrid = useGrid;
            if (useGrid)
                vhtDual->vhtDualHelper.get()->volGrid.get()->setupGrid(vhtDual->vhtDualHelper.get());
        }

        return VolBaseParam::getProcessedVol(volRequest, asset);
    }

private:

    void computeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const {
        static const string method("VolHyperTrigDual::computeImpVol");
        if ((maturities.size() != strikes.size()) ||
            (maturities.size() != impV.size())) {
            throw ModelException(method, "Size mismatch between strikes ("+
                                 Format::toString(strikes.size()) +
                                 "), maturities ("+
                                 Format::toString(maturities.size())+
                                 ") and impV ("+
                                 Format::toString(impV.size())+ ")");
        }

        for (int iMat = 0; iMat < maturities.size(); iMat++) {
            int nbStrikes = strikes[iMat].size();
            if (nbStrikes != impV[iMat].size()){
                throw ModelException(method, "Size mismatch between strikes"
                                     " & maturities for Mat " +
                                     maturities[iMat].toString() +
                                     " (n "+ Format::toString(iMat) + ")");
            }

            DateTimeArray dates = vhtDualHelper->dates;
            double tradYear = vhtDualHelper->timeMetric->yearFrac(vhtDualHelper->baseDate, maturities[iMat]);

            for (int iStrike = 0; iStrike < nbStrikes; iStrike++) {
                impV[iMat][iStrike] = calcVol(strikes[iMat][iStrike],
                                              tradYear);
            }
        }
    }

    // calculate strike ref
    double getStrikeRef(double tradYear) const {
        return vhtDualHelper->getStrikeRef(tradYear);
    }

    // given a strike, a trading time, calculate implied volatility
    double calcVol(double strike,
                   double tradYear) const{
        return vhtDualHelper->calcVol(strike,tradYear);
    }

    /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
    VolSurface* spotVolSurfaceFromStrikes(
        const CDoubleArray&   strikes) const{
        static const string method("VolHyperTrigDual::spotVolSurfaceFromStrikes");
        try{
            const VolSurface* backbone = getBackboneSurface();
            CDoubleMatrix matrix = vhtDualHelper->spotVolMatrixFromStrikes(strikes);
            /** for performance use special constructors */
            return (new VolSurface(backbone, strikes, matrix));
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }
    // reference to VolHyperTrigDualPrimer class and thus to helper class
    HelperDualSP          vhtDualHelper;

    // REGISTERED FIELDS
    string            approxPade;               // optional: 'exact' if closed-form, nber of Pade coeffs if approx

    DateTime          markedDate;           // The last marked date
    double            strikeRef;
    StringArray       strikeRefRatioBMs;        // polymorphic types not supported by IMS
    DoubleArray       strikeRefRatios;

    // ref point specifics
    string            tenorRef;                 // polymorphic types not supported by IMS
    string            tenorTime;                // optional; applies to factorSkewBMs too if any

    double            skew;
    double            convexity;
    double            levelPlusInf;
    double            levelMinusInf;

    // term structure
    StringArray       skewFactorBMs;            // polymorphic types not supported by IMS
    DoubleArray       skewFactors;
    StringArray       convexityFactorBMs;       // polymorphic types not supported by IMS
    DoubleArray       convexityFactors;

    bool              useTimeDependentCelerities;

    double            convexityPower;

    double            skewBackboneScaling; // $unregistered
    double            convexityBackboneScaling; // $unregistered

    // volgrid
    string            gridType;
    int               numGrid;
    double            TruncationStd;
    bool              useGrid;                  //transient

    // vol basis
    VarSwapBasisSP    basis;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VolHyperTrigDual, clazz);
        SUPERCLASS(VolBaseParam);
        IMPLEMENTS(IVolatilityBS);
        IMPLEMENTS(IVolatilityDVF);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(Calibrator::IAdjustable);
        EMPTY_SHELL_METHOD(defaultCtor);

        FIELD(approxPade, "Closed Form or Nber Pade Coefficients");
        FIELD_MAKE_OPTIONAL(approxPade);
        FIELD(markedDate, "Last marked date");
        FIELD_MAKE_OPTIONAL(markedDate);
        FIELD(strikeRef, "Reference strike");
        FIELD(strikeRefRatioBMs, "Reference strike ratio bencharks");
        FIELD_MAKE_OPTIONAL(strikeRefRatioBMs);
        FIELD(strikeRefRatios, "Reference strike ratios");
        FIELD_MAKE_OPTIONAL(strikeRefRatios);
        FIELD(tenorRef, "Reference tenor");
        FIELD(skew, "Reference skew");
        FIELD(convexity, "Reference convexity");
        FIELD(levelPlusInf, "Level for +infinite strike");
        FIELD(levelMinusInf, "Level for -infinite strike");
        FIELD(skewFactorBMs, "Skew Factor Bencharks");
        FIELD(skewFactors, "Skew Factors");
        FIELD(convexityFactorBMs, "Convexity Factor Bencharks");
        FIELD(convexityFactors, "Convexity Factors");
        FIELD(tenorTime, "Time of Day for Tenors");
        FIELD_MAKE_OPTIONAL(tenorTime);
        FIELD(useTimeDependentCelerities, "use time dependent celerities");
        FIELD_MAKE_OPTIONAL(useTimeDependentCelerities);
        FIELD(convexityPower, "Propagation power for convexity factors");
        FIELD_MAKE_OPTIONAL(convexityPower);
        FIELD(vhtDualHelper, "");
        FIELD_MAKE_TRANSIENT(vhtDualHelper);
        FIELD(gridType, "type of vol interpolation: NONE, LINEAR or SPLINE");
        FIELD_MAKE_OPTIONAL(gridType);
        FIELD(numGrid, "number of grid points");
        FIELD_MAKE_OPTIONAL(numGrid);
        FIELD(TruncationStd, "trancation of standard normal variable");
        FIELD_MAKE_OPTIONAL(TruncationStd);
        FIELD(useGrid, "");
        FIELD_MAKE_TRANSIENT(useGrid);
        FIELD(basis, "");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(basis);


        // add our fields and their ranges to central list
        Calibrator::IAdjustable::registerField(
            clazz, "skew",
            new InfiniteRange());
        Calibrator::IAdjustable::registerField(
            clazz, "convexity",
            new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "levelPlusInf",
            new Range(OpenBoundary(0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "levelMinusInf",
            new Range(OpenBoundary(0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "skewFactors",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
        Calibrator::IAdjustable::registerField(
            clazz, "convexityFactors",
            new Range(OpenBoundary(0.0), Infinity(Infinity::Plus)));
    }

    VolHyperTrigDual():
    VolBaseParam(TYPE),
    approxPade("EXACT"),
    strikeRef(0.0),
    tenorRef("1Y"),
    tenorTime("EOD"),
    skew(-0.02),
    convexity(0.002),
    levelPlusInf(5.0),
    levelMinusInf(10.0),
    useTimeDependentCelerities(true),
    gridType(VolHyperTrigPrimer::VolGrid::GRID_TYPE_NONE),
    numGrid(1000),
    TruncationStd(5.0),
    useGrid(false){}

    static IObject* defaultCtor(){
        return new VolHyperTrigDual();
    }

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache() {
        static const string method("VolHyperTrigDual::buildCache");
        try{
            // step I: information in backbone has to be passed to vhtHelperDual
            VolSurfaceConstSP backbone(getBackboneSurface());
            // step II: create vhtHelperDual and pass info of vht and backbone
            vhtDualHelper = HelperDualSP(new HelperDual(*this, backbone));
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields){
        update();
    }

    /** Called after (calibrator) adjustments have been made to fields */
    void update(){
        vhtDualHelper->update(skew,
                              convexity,
                              skewFactors,
                              convexityFactors,
                              levelPlusInf,
                              levelMinusInf);
    }
};

CClassConstSP const VolHyperTrigDual::TYPE =
CClass::registerClassLoadMethod("VolHyperTrigDual", typeid(VolHyperTrigDual), load);

CClassConstSP const VolHyperTrigDual::VolHyperTrigDualVolParam::TYPE =
CClass::registerClassLoadMethod("VolHyperTrigDual::VolHyperTrigDualVolParam", typeid(VolHyperTrigDualVolParam), load);

CClassConstSP const VolHyperTrigDual::HelperDual::TYPE =
CClass::registerClassLoadMethod("HelperDual", typeid(VolHyperTrigDual::HelperDual), VolHyperTrigDual::HelperDual::load);

// external symbol to allow class to be forced to be linked in
bool VolHyperTrigDualLinkIn(){
    return (VolHyperTrigDual::TYPE != 0);
}

typedef smartPtr<VolHyperTrigDual> VolHyperTrigDualSP;

/* Addin StrikeRefRatioMakerAddin */
class StrikeRefRatioMakerAddin: public CObject{
public:
    static CClassConstSP const TYPE;

private:
    CMarketDataSP       market;
    CAssetWrapper       asset;
    StringArray         strikeRefRatioBMs;
    string              time;

    // transients
    int                 time2;
    DateTime            baseDate;
    DateTimeArray       strikeRefRatioDates;

    virtual void validatePop2Object(){
        static const string method = "StrikeRefRatioMakerAddin::validatePop2Object";
        try{
            if (strikeRefRatioBMs.size() == 0){
                throw ModelException(method,
                                     "the number of reference strike ratio benchmarks should be positive");
            }
            time2 = DateTime::timeConvert(time);
            DateTime baseDate = market->GetReferenceDate();
            // ref strike ratios
            // if strike ref ratios exist, make sure they are increasing
            int nbStrikeRefRatios = strikeRefRatioBMs.size();
            strikeRefRatioDates.resize(nbStrikeRefRatios);
            int iStrikeRefRatio = 0;
            for (; iStrikeRefRatio < nbStrikeRefRatios; ++iStrikeRefRatio){
                MaturityPeriod period(strikeRefRatioBMs[iStrikeRefRatio]);
                strikeRefRatioDates[iStrikeRefRatio] = DateTime(period.toDate(baseDate).getDate(), time2);
                if (iStrikeRefRatio > 0
                    && strikeRefRatioDates[iStrikeRefRatio - 1] >= strikeRefRatioDates[iStrikeRefRatio]){
                    throw ModelException(method,
                                         "The reference strike ratio benchmarks must be increasing.\nCompare the "
                                         + Format::toString(iStrikeRefRatio)
                                         + "-th benchmark date ("
                                         + strikeRefRatioDates[iStrikeRefRatio - 1].toString()
                                         + ") with the "
                                         + Format::toString(iStrikeRefRatio + 1)
                                         + "-th benchmark date ("
                                         + strikeRefRatioDates[iStrikeRefRatio].toString()
                                         + ")");
                }
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    DoubleArraySP calcRatios(){
        static const string routine = "StrikeRefRatioMakerAddin::calcRatios";
        try {
            MarketDataFetcherSP mdf(new MarketDataFetcherLN("VolSurface"));
            NonPricingModel model(mdf);
            asset.getData(&model, market);
            int nbDates = strikeRefRatioDates.size();
            DoubleArraySP fwds(new DoubleArray(nbDates));
            asset->fwdValue(strikeRefRatioDates, *fwds);
            double spot = asset->fwdValue(baseDate);
            for (int iDate = 0; iDate < nbDates; ++iDate){
                (*fwds)[iDate] /= spot;
            }
            return fwds;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static IObjectSP calcRatios2(StrikeRefRatioMakerAddin* params) {
        return params->calcRatios();
    }

    /** for reflection */
    StrikeRefRatioMakerAddin():
    CObject(TYPE),
    time("EOD"),
    time2(0){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(StrikeRefRatioMakerAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(market, "market");
        FIELD(asset, "Asset wrapper");
        FIELD(strikeRefRatioBMs, "reference strike ratio benchmarks");
        FIELD(time, "time of day");
        FIELD_MAKE_OPTIONAL(time);
        // transients
        FIELD(time2, "");
        FIELD_MAKE_TRANSIENT(time2);
        FIELD(baseDate, "");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(strikeRefRatioDates, "");
        FIELD_MAKE_TRANSIENT(strikeRefRatioDates);

        Addin::registerClassObjectMethod("VOLHYPERTRIG_STRIKEREF_RATIO_CALC",
                                          Addin::MARKET,
                                          "Calculates the reference strike ratios needed for VolHyperTrig",
                                          TYPE,
                                          false,
                                          Addin::expandMulti,
                                          (Addin::ObjMethod*)calcRatios2);
    }

    static IObject* defaultCtor(){
        return new StrikeRefRatioMakerAddin();
    }
};

/* Addin ParameterOutputAddin */
class ParameterOutputAddin: public CObject{
public:
    static CClassConstSP const TYPE;

private:
    CMarketDataSP       market;
    CAssetWrapper       asset;
    string              volType;
    string              choice;
    double              strike;
    double              spotAtStart;

    DoubleArraySP calcParameters(){
        static const string routine = "ParameterOutputAddin::calcParameters";
        try {
            MarketDataFetcherSP mdf(new MarketDataFetcherLN(volType));
            NonPricingModel model(mdf);
            asset.getData(&model, market);
            VolHyperTrigPrimerSP vhtHelper;
            if (CString::equalsIgnoreCase(volType,"VolHyperTrig")) {
                CVolBaseSP rawVol(VolRequestRaw::copyVolBase(*asset.get()));
                VolHyperTrigSP vol(VolHyperTrigSP::dynamicCast(rawVol));
                // Dummy try/catch to avoid vc6.opt crash
                try { vol->buildCache(); } catch (...) { throw; }
                vhtHelper = vol->vhtHelper;
            } else if (CString::equalsIgnoreCase(volType,"VolHyperTrigDual")) {
                CVolBaseSP rawVol(VolRequestRaw::copyVolBase(*asset.get()));
                VolHyperTrigDualSP vol(VolHyperTrigDualSP::dynamicCast(rawVol));
                // Dummy try/catch to avoid vc6.opt crash
                try { vol->buildCache(); } catch (...) { throw; }
                vhtHelper = vol->vhtDualHelper;
            }
            int size = vhtHelper->backboneSkewParams.size();
            DoubleArraySP result(new DoubleArray(size));
            // parameter A
            if(CString::equalsIgnoreCase(choice,"A")) {
                *result = vhtHelper->backboneSkewParams;
            // parameter B
            } else if(CString::equalsIgnoreCase(choice,"B")) {
                *result = vhtHelper->backboneConvexityParams;
            // paramter C
            }else if(CString::equalsIgnoreCase(choice,"C")) {
                *result = vhtHelper->backboneSkewCelerities;
            // parameter D
            } else if(CString::equalsIgnoreCase(choice,"D")) {
                *result = vhtHelper->backboneConvexityCelerities;
            // parameter X
            } else if(CString::equalsIgnoreCase(choice,"X")) {
                for (int i = 0; i < size; i++) {
                    (*result)[i] = vhtHelper->skewFactorInterp->value(vhtHelper->tradYears[i]);
                }
            // parameter Y
            } else if(CString::equalsIgnoreCase(choice,"Y")) {
                for (int i = 0; i < size; i++) {
                    (*result)[i] = vhtHelper->convexityFactorInterp->value(vhtHelper->tradYears[i]);
                }
            // input for skew (term structure)
            } else if(CString::equalsIgnoreCase(choice,"skewInput")) {
                static const double ln100To90 = -log(0.9);
                double dAtmVol_dx = vhtHelper->skew / ln100To90;
                for (int i = 0; i < size; i++) {
                    (*result)[i] = dAtmVol_dx * vhtHelper->skewFactorInterp->value(vhtHelper->tradYears[i])
                                    / sqrt(vhtHelper->tradYears[i] / vhtHelper->tradYearRef) * (-log(0.9));
                }
            // output for skew (term structure) for various strikes
            } else if(CString::equalsIgnoreCase(choice,"skewOutput")) {
                *result = computeSkewOutput(vhtHelper, strike);
            // output for skew factors (term struuctures) for various strikes
            } else if(CString::equalsIgnoreCase(choice,"skewFactorOutput")) {
                DoubleArray skewOutput = computeSkewOutput(vhtHelper, strike);
                for (int i = 0; i < size; i++) {
                    (*result)[i] = skewOutput[i] / skewOutput[vhtHelper->dateRefIdx]
                        * sqrt(vhtHelper->tradYears[i] / vhtHelper->tradYears[vhtHelper->dateRefIdx]);
                }
            // input for convexity (term structure)
            } else if(CString::equalsIgnoreCase(choice,"convexityInput")) {
                static const double ln100To90 = -log(0.9);
                double d2AtmVol_dx2 = vhtHelper->convexity / (ln100To90 * ln100To90);
                for (int i = 0; i < size; i++) {
                    (*result)[i] = d2AtmVol_dx2 * vhtHelper->convexityFactorInterp->value(vhtHelper->tradYears[i])
                                        / sqrt(vhtHelper->tradYears[i] / vhtHelper->tradYearRef) * (-log(0.9))*(-log(0.9));
                }
            // output for convexity (term structure) for various strikes
            } else if(CString::equalsIgnoreCase(choice,"convexityOutput")) {
               *result = computeConvexityOutput(vhtHelper, strike);
            // output for convexity factors (term structure) for various strikes
            } else if(CString::equalsIgnoreCase(choice,"convexityFactorOutput")) {
                DoubleArray convexityOutput = computeConvexityOutput(vhtHelper, strike);
                for (int i = 0; i < size; i++) {
                    (*result)[i] = convexityOutput[i] / convexityOutput[vhtHelper->dateRefIdx]
                        * pow(vhtHelper->tradYears[i] / vhtHelper->tradYears[vhtHelper->dateRefIdx], vhtHelper->convexityPower);
                }
            } else if(CString::equalsIgnoreCase(choice,"levels")) {
                // levelPlusInf
                (*result)[0] = 1 + vhtHelper->backboneSkewParams[vhtHelper->dateRefIdx]
                    + vhtHelper->backboneConvexityParams[vhtHelper->dateRefIdx];
                (*result)[1] = 1 - vhtHelper->backboneSkewParams[vhtHelper->dateRefIdx]
                    + vhtHelper->backboneConvexityParams[vhtHelper->dateRefIdx];
            } else if(CString::equalsIgnoreCase(choice, "atmVols")) {
                for (int i = 0; i < size; i++) {
                    (*result)[i] = vhtHelper->calcVol(spotAtStart, vhtHelper->tradYears[i]);
                }
            } else if(CString::equalsIgnoreCase(choice, "strikeRefRatios")) {
                for (int i = 0; i < size; i++) {
                    (*result)[i] = vhtHelper->getStrikeRef(vhtHelper->tradYears[i]);
                }
            }
            return result;
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    DoubleArray computeSkewOutput(VolHyperTrigPrimerSP vhtHelper, double strike) {
        int size = vhtHelper->backboneSkewParams.size();
        DoubleArray result(size);
        for (int i = 0; i < size; i++) {
            double tradYear = vhtHelper->tradYears[i];
            double spotNow = vhtHelper->getStrikeRef(tradYear);
            double volAtStrikeDown = vhtHelper->calcVol(spotNow * strike, tradYear);
            double volAtStrikeUp = vhtHelper->calcVol(spotNow / strike, tradYear);
            result[i] = (volAtStrikeUp - volAtStrikeDown) / 2.0 / (-log(strike))*(-log(0.9));
        }
        return result;
    }

    DoubleArray computeConvexityOutput(VolHyperTrigPrimerSP vhtHelper, double strike) {
        int size = vhtHelper->backboneSkewParams.size();
        DoubleArray result(size);
        for (int i = 0; i < size; i++) {
           double tradYear = vhtHelper->tradYears[i];
           double spotNow = vhtHelper->getStrikeRef(tradYear);
           double volAtStrikeDown = vhtHelper->calcVol(spotNow * strike, tradYear);
           double volAtATM = vhtHelper->calcVol(spotNow, tradYear);
           double volAtStrikeUp = vhtHelper->calcVol(spotNow / strike, tradYear);
           result[i] = (volAtStrikeUp - 2.0 * volAtATM + volAtStrikeDown)
               / (-log(strike)) / (-log(strike)) * (-log(0.9)) * (-log(0.9));
        }
        return result;
    }

    static IObjectSP callCalcParameters(ParameterOutputAddin* params) {
        return params->calcParameters();
    }

    /** for reflection */
    ParameterOutputAddin():
    CObject(TYPE),
    strike(0.9) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ParameterOutputAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultParameterOutputAddin);
        FIELD(market, "market");
        FIELD(asset, "asset");
        FIELD(volType, "volType");
        FIELD(choice, "choice");
        FIELD(strike, "strike");
        FIELD_MAKE_OPTIONAL(strike);
        FIELD(spotAtStart, "spotAtStart");
        FIELD_MAKE_OPTIONAL(spotAtStart);

        Addin::registerClassObjectMethod("VOLHYPERTRIG_PARAMETER_OUTPUT",
                                          Addin::MARKET,
                                          "returns parameters",
                                          TYPE,
                                          false,
                                          Addin::expandMulti,
                                          (Addin::ObjMethod*)callCalcParameters);
    }

    static IObject* defaultParameterOutputAddin(){
        return new ParameterOutputAddin();
    }
};

/* ConvertAtmVolsAddin */
class ConvertAtmVolsAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    MarketDataSP        market;
    CAssetWrapper       asset;
    double              spotAtStart;
    DateTimeArray         benchMarkDates;
    StringArray            strikeRefRatioBMs;
    DoubleArray            strikeRefRatios;
    string              volType; // optional (default = VolPreferred)

private:
    DoubleArraySP convertAtmVols(){
        static const string method = "ConvertAtmVolsAddin::convertAtmVols";
        try {
            DateTime valueDate;
            market->GetReferenceDate(valueDate);

            // need a model to get data out of the market
            MarketDataFetcherSP mdf(new MarketDataFetcherLN(volType));
            NonPricingModel model(mdf);
            asset.getData(&model, market);

            string name = asset->getName();

            // need a vol request in order to get a time metric
            LinearStrikeVolRequest volRequest(spotAtStart, valueDate, valueDate, false);
            CVolProcessedBSSP vol(asset->getProcessedVol(&volRequest));
            TimeMetricConstSP timeMetric = vol->GetTimeMetric();

            // check inputs for bench mark dates
            int sizeNew = benchMarkDates.size();
            if (sizeNew == 0) {
                throw ModelException("benchmarks of vol surface are missing");
            }

        DoubleArray benchMarkTradYears(sizeNew);

            int i = 0;
            for (i = 0; i < sizeNew; ++i){
                if (i > 0 && benchMarkDates[i-1] >= benchMarkDates[i]){
                    throw ModelException(method,
                                         "The benchmark dates must be increasing.\nCompare the "
                                         + Format::toString(i)
                                         + "-th benchmark date ("
                                         + benchMarkDates[i-1].toString()
                                         + ") with the "
                                         + Format::toString(i+1)
                                         + "-th benchmark date ("
                                         + benchMarkDates[i].toString()
                                         + ")");
                }
                benchMarkTradYears[i] = timeMetric->yearFrac(valueDate, benchMarkDates[i]);
            }

            int tenorTimeInt = DateTime::timeConvert("EOD"); // arbitrary
            // check input for strike ref ratios
            if (strikeRefRatios.size() != strikeRefRatioBMs.size()){
                throw ModelException(method,
                                     "The number of reference strike ratios ("
                                     + Format::toString(strikeRefRatios.size())
                                     + ") must be equal to the number of reference strike ratio benchmarks ("
                                     + Format::toString(strikeRefRatioBMs.size())
                                     + ")");
            }

            // if there are no strike ref ratios, then we are finished
            DoubleArraySP result(new DoubleArray(sizeNew));
            if ((strikeRefRatioBMs.size() == 0) | (strikeRefRatios.size() == 0)) {
                for (i = 0; i < sizeNew; ++i) {
                    LinearStrikeVolRequest volRequest(spotAtStart, valueDate, benchMarkDates[i], false);
                    CVolProcessedBSSP vol(asset->getProcessedVol(&volRequest));
                    (*result)[i] = vol->CalcVol(valueDate, benchMarkDates[i]);
                }
                return result;
            }

                        int size = strikeRefRatios.size();
            Maths::checkPositive(strikeRefRatios, "strikeRefRatios");

            DateTimeArray strikeRefRatioDates(size);
            DoubleArray strikeRefRatioTradYears(size);

            for (i = 0; i < strikeRefRatioBMs.size(); ++i){
                MaturityPeriod period(strikeRefRatioBMs[i]);
                strikeRefRatioDates[i] = DateTime(period.toDate(valueDate).getDate(), tenorTimeInt);
                if (i > 0 && strikeRefRatioDates[i-1] >= strikeRefRatioDates[i]){
                    throw ModelException(method,
                                         "The strike ref ratio benchmarks must be increasing.\nCompare the "
                                         + Format::toString(i)
                                         + "-th benchmark date ("
                                         + strikeRefRatioDates[i-1].toString()
                                         + ") with the "
                                         + Format::toString(i+1)
                                         + "-th benchmark date ("
                                         + strikeRefRatioDates[i].toString()
                                         + ")");
                }
                strikeRefRatioTradYears[i] = timeMetric->yearFrac(valueDate, strikeRefRatioDates[i]);
            }

            // create a linear interpolator for strike refs
            LinearInterpolantConstSP    strikeRefRatioInterp;
            LinearInterpolator inter;
            strikeRefRatioInterp = LinearInterpolantConstSP::dynamicCast(
                                        inter.computeInterp(strikeRefRatioTradYears,strikeRefRatios));

            // now calc the vol atm at the strike ref ratios
            for (i=0; i < sizeNew; ++i){
                // determine current spot at strike ref ratio (ATM at strike ref ratio)
                int position = i;
                double currentStrikeRefRatio =
                    strikeRefRatioInterp->valueWithGuess(benchMarkTradYears[i], position);
                double currentLevel = spotAtStart * currentStrikeRefRatio;

                LinearStrikeVolRequest volRequest(currentLevel, valueDate, benchMarkDates[i], false);
                CVolProcessedBSSP vol(asset->getProcessedVol(&volRequest));
                (*result)[i] = vol->CalcVol(valueDate, benchMarkDates[i]);
            }

            return result;

        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    static IObjectSP callConvertAtmVols(ConvertAtmVolsAddin* params) {
        return params->convertAtmVols();
    }

    /** for reflection */
    ConvertAtmVolsAddin():
    CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ConvertAtmVolsAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConvertAtmVolsAddin);
        FIELD(market, "");
        FIELD(asset, "");
        FIELD(spotAtStart, "");
        FIELD(benchMarkDates, "");
        FIELD(strikeRefRatioBMs, "");
        FIELD(strikeRefRatios, "");
        FIELD(volType, "");
        FIELD_MAKE_OPTIONAL(volType);

        Addin::registerClassObjectMethod("VOLHYPERTRIG_CONVERT_ATM_VOLS",
                                          Addin::MARKET,
                                          "returns atm vols at strike ref ratios",
                                          TYPE,
                                          false,
                                          Addin::expandMulti,
                                          (Addin::ObjMethod*)callConvertAtmVols);
    }

    static IObject* defaultConvertAtmVolsAddin(){
        return new ConvertAtmVolsAddin();
    }
};

/**  DeltaBasedSkewConvexityLeastSquareFit objFunc class */
class DeltaBasedSkewConvexityLeastSquareFit: public Calibrator::ObjFuncLeastSquare{
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object(){
        static const string method = "DeltaBasedSkewConvexityLeastSquareFit::validatePop2Object";
        try {
            // check size of inputs
            int sizeBM = BMs.size();

            if (sizeBM == 0) {
                throw ModelException("benchmarks of skew and convexity are missing");
            }

            if (!convexityOnly){
                int sizeSkew = deltaBasedSkewTgt.size();
                if(sizeSkew == 0){
                    throw ModelException("skew targets are missing");
                } else if (sizeSkew != sizeBM) {
                    throw ModelException("sizes of benchmarks and skew targets must be equal");
                }
            }

            int sizeConvexity = deltaBasedConvexityTgt.size();
            if (sizeConvexity == 0) {
                throw ModelException("convexity targets are missing");
            } else if (sizeConvexity != sizeBM) {
                throw ModelException("sizes of benchmarks and convexity targets must be equal");
            }

            if (Maths::isNegative(deltaForSkew)||Maths::isPositive(deltaForSkew-1.0)){
                throw ModelException("deltaForSkew must be between 0 and 100%.");
            }
            if (Maths::isNegative(deltaForConvexity)||Maths::isPositive(deltaForConvexity-1.0)){
                throw ModelException("deltaForConvexity must be between 0 and 100%.");
            }

            if (convexityOnly) {
                nbFuncs = sizeBM;
            } else {
                nbFuncs = 2 * sizeBM;
            }

            if (volType == "") {
                volType = "VolHyperTrig";
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    // extra validation
    void validate() {
        static const string method = "DeltaBasedSkewConvexityLeastSquareFit::validate";
        try {
            // check dates
            dates.resize(BMs.size());
            int     tenorTime = DateTime::timeConvert("EOD"); //arbitrary
            for (int iBM = 0; iBM < BMs.size(); ++iBM){
                MaturityPeriod period(BMs[iBM]);
                dates[iBM] = DateTime(period.toDate(valueDate).getDate(), tenorTime);
                if (iBM > 0 && dates[iBM-1] >= dates[iBM]){
                    throw ModelException(method,
                                         "The benchmark dates must be increasing.\nCompare the "
                                         + Format::toString(iBM)
                                         + "-th benchmark date ("
                                         + dates[iBM-1].toString()
                                         + ") with the "
                                         + Format::toString(iBM+1)
                                         + "-th benchmark date ("
                                         + dates[iBM].toString()
                                         + ")");
                }
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    void getMarket(MarketData* market){
        static const string method = "DeltaBasedSkewConvexityLeastSquareFit::getMarket";
        try{
            // need a model to get data out of the market
            MarketDataFetcherSP mdf(new MarketDataFetcherLN(volType));
            NonPricingModel model(mdf);
            CAsset::getAssetMarketData(&model,
                           market,
                           CAsset::CCY_TREATMENT_NONE,
                           discount,
                           asset);
            discount.getData(&model, market);
            market->GetReferenceDate(valueDate);
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    IObjectSP getAdjustableGroup(){
        return IObjectSP::attachToRef(this);
    }

    int getNbFuncs() const{
        return nbFuncs;
    }

    void calcValue(CDoubleArray& funcvals) const {
        static const string method = "DeltaBasedSkewConvexityLeastSquareFit::calcValue";
        try{
            calcDeltaBasedSkewConvexity();
            int iFunc = 0;
            for (int iBM = 0; iBM < BMs.size(); ++iBM){
                if(!convexityOnly){
                    funcvals[iFunc++] = deltaBasedSkew[iBM] - deltaBasedSkewTgt[iBM];
                }
                funcvals[iFunc++] = deltaBasedConvexity[iBM] - deltaBasedConvexityTgt[iBM];
            }
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }

    // for reflection
    DeltaBasedSkewConvexityLeastSquareFit():
    Calibrator::ObjFuncLeastSquare(TYPE),
    nbFuncs(2){}

    DeltaBasedSkewConvexityLeastSquareFit(CAssetWrapper           asset,
                                        YieldCurveWrapper       discount,
                                        string                  volType,
                                        double                  deltaForSkew,
                                        double                  deltaForConvexity,
                                        bool                    twoSideSkew,
                                        bool                    convexityOnly,
                                        StringArray             BMs,
                                        DoubleArray             deltaBasedSkewTgt,
                                        DoubleArray             deltaBasedConvexityTgt):
    Calibrator::ObjFuncLeastSquare(TYPE),
    asset(asset),
    discount(discount),
    volType(volType),
    deltaForSkew(deltaForSkew),
    deltaForConvexity(deltaForConvexity),
    twoSideSkew(twoSideSkew),
    convexityOnly(convexityOnly),
    BMs(BMs),
    deltaBasedSkewTgt(deltaBasedSkewTgt),
    deltaBasedConvexityTgt(deltaBasedConvexityTgt){
        validatePop2Object();
    }

    // calc delta based Skew and Convexity
    void calcDeltaBasedSkewConvexity() const {
        static const string method = "DeltaBasedSkewConvexityLeastSquareFit::calcDeltaBasedSkewConvexity";
        try {
            int size = BMs.size();
            strikes.resize(4);
            for (int i=0; i<4; i++) {
                strikes[i].resize(size);
            }
            deltaBasedSkew.resize(size);
            deltaBasedConvexity.resize(size);

             // Delta implied strike
            string className("CVanilla::DeltaImpliedStrikeMakerLN");
            IDeltaToStrikeMakerSP maker(dynamic_cast<IDeltaToStrikeMaker*>(CClass::forName(className)->newInstance()));
            if(!maker) {
                throw ModelException(method, "Failed to create CVanilla::DeltaImpliedStrikeMakerLN");
            }

            DoubleArray fwds(size);
            asset->fwdValue(dates, fwds);

            double spot = asset->getSpot();
            double strikeAbsAcc = spot * 0.01 / 100.0;   // 0.01% of current spot
            double lowerStrike = spot * 0.3;
            double upperStrike = spot * 3.0;
            InstrumentSettlementSP  instSetllement(new CashSettlePeriod(0));

            int iBM = 0;
            for (; iBM < size; ++iBM){
                // lowStrikes for put
                IDeltaToStrikeMaker::IDeltaToStrikeConstSP calc1(maker->make(
                    valueDate,
                    dates[iBM],
                    false,
                    asset.get(),
                    discount.get(),
                    instSetllement.get(),
                    Delta::DEFAULT_SHIFT,
                    -deltaForConvexity,
                    volType,
                    true));

                strikes[ConvexityDeltaPut][iBM] = calc1->calcStrike(lowerStrike, upperStrike, strikeAbsAcc);

                // highStrikes for call
                IDeltaToStrikeMaker::IDeltaToStrikeConstSP calc2(maker->make(
                    valueDate,
                    dates[iBM],
                    true,
                    asset.get(),
                    discount.get(),
                    instSetllement.get(),
                    Delta::DEFAULT_SHIFT,
                    deltaForConvexity,
                    volType,
                    true));

                strikes[ConvexityDeltaCall][iBM] = calc2->calcStrike(lowerStrike, upperStrike, strikeAbsAcc);

                if(!convexityOnly){
                    if (Maths::isZero(deltaForSkew-deltaForConvexity)){
                        strikes[SkewDeltaPut][iBM] = strikes[ConvexityDeltaPut][iBM];
                        strikes[SkewDeltaCall][iBM] = strikes[ConvexityDeltaCall][iBM];
                    } else {
                        // lowStrikes for put
                        IDeltaToStrikeMaker::IDeltaToStrikeConstSP calc1(maker->make(
                            valueDate,
                            dates[iBM],
                            false,
                            asset.get(),
                            discount.get(),
                            instSetllement.get(),
                            Delta::DEFAULT_SHIFT,
                            -deltaForSkew,
                            volType,
                            true));

                        strikes[SkewDeltaPut][iBM] = calc1->calcStrike(lowerStrike, upperStrike, strikeAbsAcc);

                        // highStrikes for call
                        IDeltaToStrikeMaker::IDeltaToStrikeConstSP calc2(maker->make(
                            valueDate,
                            dates[iBM],
                            true,
                            asset.get(),
                            discount.get(),
                            instSetllement.get(),
                            Delta::DEFAULT_SHIFT,
                            deltaForSkew,
                            volType,
                            true));

                        strikes[SkewDeltaCall][iBM] = calc2->calcStrike(lowerStrike, upperStrike, strikeAbsAcc);
                    }
                }
            }

            for (iBM = 0 ; iBM < size; ++iBM){
                double atmVol = calcVol(fwds[iBM],dates[iBM]);
                double lowVolSkew = calcVol(strikes[SkewDeltaPut][iBM],dates[iBM]);
                double highVolSkew = calcVol(strikes[SkewDeltaCall][iBM],dates[iBM]);
                double lowVolConv = calcVol(strikes[ConvexityDeltaPut][iBM],dates[iBM]);
                double highVolConv = calcVol(strikes[ConvexityDeltaCall][iBM],dates[iBM]);

                if (twoSideSkew) {
                    deltaBasedSkew[iBM] = highVolSkew - lowVolSkew;
                } else {
                    deltaBasedSkew[iBM] = atmVol - lowVolSkew;
                }
                deltaBasedConvexity[iBM] = lowVolConv + highVolConv - 2*atmVol;
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

private:
    // calc vol given strike and date
    double calcVol(double     strike,
                   DateTime   date) const{
        static const string method = "DeltaBasedSkewConvexityLeastSquareFit::calcVol";
        LinearStrikeVolRequest volRequest(strike, valueDate, date, false);
        volRequest.allowNegativeFwdVar(true);
        CVolProcessedBSSP vol(asset->getProcessedVol(&volRequest));
        return vol->CalcVol(valueDate, date);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DeltaBasedSkewConvexityLeastSquareFit, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDeltaBasedSkewConvexityLeastSquareFit);

        FIELD(asset, "");
        FIELD(discount, "");
        FIELD(volType, "");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(deltaForSkew, "");
        FIELD(deltaForConvexity, "");
        FIELD(twoSideSkew, "");
        FIELD(convexityOnly, "");
        FIELD(BMs, "");
        FIELD(deltaBasedSkewTgt, "");
        FIELD(deltaBasedConvexityTgt, "");

        // transients
        FIELD(valueDate, "");
        FIELD_MAKE_TRANSIENT(valueDate);
        FIELD(dates, "");
        FIELD_MAKE_TRANSIENT(dates);
        FIELD(strikes, "");
        FIELD_MAKE_TRANSIENT(strikes);
        FIELD(deltaBasedSkew, "");
        FIELD_MAKE_TRANSIENT(deltaBasedSkew);
        FIELD(deltaBasedConvexity, "");
        FIELD_MAKE_TRANSIENT(deltaBasedConvexity);
        FIELD(nbFuncs, "");
        FIELD_MAKE_TRANSIENT(nbFuncs);
    }

    static IObject* defaultDeltaBasedSkewConvexityLeastSquareFit(){
        return new DeltaBasedSkewConvexityLeastSquareFit();
    }

    static const int SkewDeltaPut;
    static const int SkewDeltaCall;
    static const int ConvexityDeltaPut;
    static const int ConvexityDeltaCall;

    //registered fields
    CAssetWrapper           asset;
    YieldCurveWrapper       discount;
    string                  volType;

    double                  deltaForSkew;
    double                  deltaForConvexity;
    bool                    twoSideSkew;
    bool                    convexityOnly;
    StringArray             BMs;
    DoubleArray             deltaBasedSkewTgt;
    DoubleArray             deltaBasedConvexityTgt;

    // transcient fields
    int                         nbFuncs;
    DateTime                    valueDate;
    DateTimeArray               dates;
    mutable DoubleArrayArray    strikes;
    mutable DoubleArray         deltaBasedSkew;
    mutable DoubleArray         deltaBasedConvexity;
};
typedef smartPtr<DeltaBasedSkewConvexityLeastSquareFit> DeltaBasedSkewConvexityLeastSquareFitSP;

const int DeltaBasedSkewConvexityLeastSquareFit::SkewDeltaPut              = 0;
const int DeltaBasedSkewConvexityLeastSquareFit::SkewDeltaCall             = 1;
const int DeltaBasedSkewConvexityLeastSquareFit::ConvexityDeltaPut         = 2;
const int DeltaBasedSkewConvexityLeastSquareFit::ConvexityDeltaCall        = 3;

/* DeltaBasedSkewConvexityAddin */
class DeltaBasedSkewConvexityAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    MarketDataSP            market;
    CAssetWrapper           asset;
    YieldCurveWrapper       discount;
    double                  deltaForSkew;           // delta to compute delta based skew
    double                  deltaForConvexity;      // delta to compute delta based convexity
    bool                    twoSideSkew;
    bool                    convexityOnly;
    StringArray             BMs;
    DoubleArray             deltaBasedSkewTgt;
    DoubleArray             deltaBasedConvexityTgt;
    OptimizerNDSP           optimizer;

    // transcient fields
    string                  volType;
    string                  refBM;
    StringArray             fieldNames;             // fields to be calibrated
    CBoolArray              arrayFields;            // whether the field is a double or a double array
    IntArray                arrayIndexBM;           // index of BMs elements in skewBMs

    virtual void validatePop2Object(){
        static const string method = "DeltaBasedSkewConvexityAddin::validatePop2Object";
        try {
            // check whether "1Y" point needs to be matched
            volType = "VolHyperTrig";
            refBM = "1Y";
            bool    hasRef = false;
            int iBM = 0;
            for(; iBM < BMs.size(); iBM++){
                if(CString::equalsIgnoreCase(refBM, BMs[iBM])){
                    hasRef = true;
                }
            }

            // set fieldNames and arrayFields
            if (hasRef){
                if (convexityOnly){
                    fieldNames.resize(2);
                    arrayFields.resize(2);
                    fieldNames[0] = "convexityFactors";
                    fieldNames[1] = "convexity";
                    arrayFields[0] = true;
                    arrayFields[1] = false;
                } else {
                    fieldNames.resize(4);
                    arrayFields.resize(4);
                    fieldNames[0] = "skewFactors";
                    fieldNames[1] = "convexityFactors";
                    fieldNames[2] = "skew";
                    fieldNames[3] = "convexity";
                    arrayFields[0] = true;
                    arrayFields[1] = true;
                    arrayFields[2] = false;
                    arrayFields[3] = false;
                }
            } else {
                if (convexityOnly){
                    fieldNames.resize(1);
                    arrayFields.resize(1);
                    fieldNames[0] = "convexityFactors";
                    arrayFields[0] = true;
                } else {
                    fieldNames.resize(2);
                    arrayFields.resize(2);
                    fieldNames[0] = "skewFactors";
                    fieldNames[1] = "convexityFactors";
                    arrayFields[0] = true;
                    arrayFields[1] = true;
                }
            }

            // need a model to get data out of the market
            MarketDataFetcherSP mdf(new MarketDataFetcherLN(volType));
            NonPricingModel model(mdf);
            asset.getData(&model, market);
            string          name = asset->getName();
            MarketObjectSP  obj = market->GetData(name,VolHyperTrig::TYPE);
            VolHyperTrigSP  vht = VolHyperTrigSP::dynamicCast(obj);
            StringArray     skewFactorBMs = vht->getSkewFactorBMs();
            StringArray     convexityFactorBMs = vht->getConvexityFactorBMs();

            // set arrayIndexBM
            // assume same benchmarks for skew and convexity factors
            if(skewFactorBMs.size() != convexityFactorBMs.size()) {
                throw ModelException("skewFactor and convexityFactor must have the same benchmarks");
            }
            int iSkewFactorBM = 0;
            for (; iSkewFactorBM < skewFactorBMs.size(); iSkewFactorBM++){
                if(!CString::equalsIgnoreCase(skewFactorBMs[iSkewFactorBM],convexityFactorBMs[iSkewFactorBM])){
                    throw ModelException("skewFactors and convexityFactors must have the same benchmarks");
                }
            }

            arrayIndexBM.resize(BMs.size(),skewFactorBMs.size());
            for (iBM = 0;iBM < BMs.size(); iBM++){
                if(!CString::equalsIgnoreCase(BMs[iBM],refBM)){
                    for(iSkewFactorBM = 0; iSkewFactorBM < skewFactorBMs.size(); iSkewFactorBM++){
                        if(CString::equalsIgnoreCase(BMs[iBM],skewFactorBMs[iSkewFactorBM])){
                            arrayIndexBM[iBM] = iSkewFactorBM;
                        }
                    }
                    if (arrayIndexBM[iBM] == skewFactorBMs.size()){
                        throw ModelException("BMs must be a subset of skewFactorBMs");
                    }
                }
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    IObjectSP convertSkewConvexity(){
        static const string method("DeltaBasedSkewConvexityAddin::convertSkewConvexity");

        // calibrator
        Calibrator          calibrator(optimizer);

        // objFunc
        DeltaBasedSkewConvexityLeastSquareFit      objFunc(asset,
                                                           discount,
                                                           volType,
                                                           deltaForSkew,
                                                           deltaForConvexity,
                                                           twoSideSkew,
                                                           convexityOnly,
                                                           BMs,
                                                           deltaBasedSkewTgt,
                                                           deltaBasedConvexityTgt);

        // if there's a market, give objFunc a chance to get it
        // otherwise, trust (!) that objFunc already has all the market
        // data it needs
        if (market.get()){
            objFunc.getMarket(market.get());
        }

        // further validation irrespective of getMarket having been
        // called or not
        objFunc.validate();

        // ids
        string name = asset->getName();
        Calibrator::InstanceIDArray ids;
        int iField = 0;
        ids.reserve(fieldNames.size() * BMs.size());
        for (; iField < fieldNames.size(); iField++){
            if(arrayFields[iField]) {
                int iBM = 0;
                for (; iBM < BMs.size(); iBM++){
                    // currently all doubles at the moment
                    if(!CString::equalsIgnoreCase(BMs[iBM],refBM)){
                        ids.push_back(Calibrator::InstanceIDSP(
                            new Calibrator::InstanceIDDb(volType,
                                                         name,
                                                         fieldNames[iField],
                                                         false,
                                                         1.0,
                                                         arrayIndexBM[iBM])));
                    }
                }
            } else {
                ids.push_back(Calibrator::InstanceIDSP(
                new Calibrator::InstanceIDDb(volType,
                                             name,
                                             fieldNames[iField],
                                             false,
                                             0.01,
                                             0)));
            }
        }

        // run the calibrator
        return calibrator.run(objFunc, ids);
    }

    static IObjectSP convert(DeltaBasedSkewConvexityAddin* params) {
        return params->convertSkewConvexity();
    }

    /** for reflection */
    DeltaBasedSkewConvexityAddin():
    CObject(TYPE) {}

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DeltaBasedSkewConvexityAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDeltaBasedSkewConvexityAddin);

        FIELD(market, "");
        FIELD(asset, "");
        FIELD(discount, "");
        FIELD(deltaForSkew, "");
        FIELD(deltaForConvexity, "");
        FIELD(twoSideSkew, "");
        FIELD(convexityOnly, "");
        FIELD(BMs, "");
        FIELD(deltaBasedSkewTgt, "");
        FIELD(deltaBasedConvexityTgt, "");
        FIELD(optimizer, "");

        // transients
        FIELD(volType, "");
        FIELD_MAKE_TRANSIENT(volType);
        FIELD(refBM, "");
        FIELD_MAKE_TRANSIENT(refBM);
        FIELD(fieldNames, "");
        FIELD_MAKE_TRANSIENT(fieldNames);
        FIELD(arrayFields, "");
        FIELD_MAKE_TRANSIENT(arrayFields);
        FIELD(arrayIndexBM, "");
        FIELD_MAKE_TRANSIENT(arrayIndexBM);

        Addin::registerClassObjectMethod("VOLHYPERTRIG_CONVERT_SKEW_CONVEXITY",
                                          Addin::MARKET,
                                          "returns VHT object matching delta based skew and convexity targets",
                                          TYPE,
                                          true,
                                          Addin::returnHandle,
                                          (Addin::ObjMethod*)convert);
    }

    static IObject* defaultDeltaBasedSkewConvexityAddin(){
        return new DeltaBasedSkewConvexityAddin();
    }
};


CClassConstSP const StrikeRefRatioMakerAddin::TYPE =
CClass::registerClassLoadMethod("StrikeRefRatioMakerAddin", typeid(StrikeRefRatioMakerAddin), load);

CClassConstSP const ParameterOutputAddin::TYPE =
CClass::registerClassLoadMethod("ParameterOutputAddin", typeid(ParameterOutputAddin), load);

CClassConstSP const ConvertAtmVolsAddin::TYPE =
CClass::registerClassLoadMethod("ConvertAtmVolsAddin", typeid(ConvertAtmVolsAddin), load);

CClassConstSP const DeltaBasedSkewConvexityLeastSquareFit::TYPE =
CClass::registerClassLoadMethod("DeltaBasedSkewConvexityLeastSquareFit", typeid(DeltaBasedSkewConvexityLeastSquareFit), load);

CClassConstSP const DeltaBasedSkewConvexityAddin::TYPE =
CClass::registerClassLoadMethod("DeltaBasedSkewConvexityAddin", typeid(DeltaBasedSkewConvexityAddin), load);

DRLIB_END_NAMESPACE
