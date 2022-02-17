//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : HitSample.cpp
//
//   Description : Class for computing probabilities and sampling of hitting times etc.
//
//   Date        : February 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/HitSample.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/RootFinder.hpp"
#include "edginc/Addin.hpp"


DRLIB_BEGIN_NAMESPACE


CClassConstSP const HitSample::TYPE = CClass::registerInterfaceLoadMethod(
    "HitSample", typeid(HitSample), HitSample::load);


void HitSample::load(CClassSP &clazz){
    REGISTER_INTERFACE(HitSample, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}


///////////////////////////////////////////////////////////

// Initialize static member variables
const double HitDistributionBM::INVERSION_TOLERANCE = 0.0001;
const double HitDistributionBM::HUGE_HITTING_TIME = 9999999999.0;


HitDistributionBM::HitDistributionBM(): CObject(TYPE), hitProb(0.0), 
commonNum(0.0), num1(0.0), num2(0.0), hitStart(false),
hitEnd(false), zeroVol(false), zeroLength(false) {}


HitDistributionBM::HitDistributionBM(double valueStart, double valueEnd, 
                                     double barrierStart, double barrierEnd, 
                                     double vol, double length, bool isUp):
CObject(TYPE), valueStart(valueStart), valueEnd(valueEnd), barrierStart(barrierStart), 
barrierEnd(barrierEnd), vol(vol), length(length), isUp(isUp),
hitProb(0.0), commonNum(0.0), num1(0.0), num2(0.0), hitStart(false),
hitEnd(false), zeroVol(false), zeroLength(false) {
    validatePop2Object();
}


HitDistributionBM::HitDistributionBM(double barrierStart, double barrierEnd, 
                                     double vol, double length, bool isUp):
CObject(TYPE), valueStart(0.0), valueEnd(0.0), barrierStart(barrierStart), 
barrierEnd(barrierEnd), vol(vol), length(length), isUp(isUp),
hitProb(0.0), commonNum(0.0), num1(0.0), num2(0.0), hitStart(false),
hitEnd(false), zeroVol(false), zeroLength(false) {
    validatePop2Object();
}

    
    double HitDistributionBM::hitTimeProb(double t) const {
    // First deal with degenerate cases

    // Just to define the function for all t < 0.0
    // if( Maths::isNegative(t) ) {
    if( t < 0.0 ) {
        return 0.0;
    }
    
    // hitProb is correct for hit at start or t > length
    if(hitStart || zeroLength || t >= length
        /*!Maths::isNegative(t - length) */) {
        return hitProb;
    }
    
    // No hit start, no hit end, zero vol
    if(zeroVol) {
        // Hit at end with zero vol
        if(hitEnd) {
            double linearWeight = (valueEnd - barrierEnd) / num1;
            double linearTime = (1.0 - linearWeight) * length;
            // double dist = !Maths::isNegative(t - linearTime) ? 1.0 : 0.0;
            double dist = (t <= linearTime) ? 0.0 : 1.0;
            return dist;
        } else {
            return 0.0;
        }
    }

    // Non degenerate cases now
    
    // Zero probability of exiting instantaneously in diffusion (cover t < 0.0 too)
    // if(!Maths::isPositive(t)) {
    if(t == 0.0) {
        return 0.0;
    }

    // Now deal with non-degenerate case
    double denom = vol * sqrt(t * length * (length - t));
        
    double w    = (t * num1 - commonNum) / denom;
    double wbar = (t * num2 - commonNum) / denom;

    if(!isUp) {
        // Down barrier
        w    = -w;
        wbar = -wbar;
    }

    double prob = N1(w) + hitProb * N1(wbar);
    return prob;
}
    

double HitDistributionBM::hitTimeGivenHitProb(double t) const {
    double condProb = hitTimeProb(t);
    if (!hitEnd && hitProb != 0.0) {
        // Scale distribution so that it integrates to 1.0
        condProb /= hitProb;
    }
    
    return condProb;
}


double HitDistributionBM::getHittingProb() const {
    return hitProb;
}


int HitDistributionBM::hitNoHitSample(double u) const {
    // Return u <= hitProb
    // return !Maths::isPositive(u - hitProb);
    if(hitProb == 0.0) {
        return 0;
    } else {
        return (u <= hitProb) ? 1 : 0;
    }
}


double HitDistributionBM::hitTimeSample(double u) const {
    static const string routine = "HitDistributionBM::hitTimeSample";
    
    try {
        // if(Maths::isPositive(u - hitProb)) {
        if(u > hitProb) {
            return HUGE_HITTING_TIME;
        }
        
        // Brackets
        double tLow = 0.0;
        double tHigh = length;
        
        // Create wrappers
        if(!probTime) {
            HitDistWrapper hitTimeProbTemp(this, &HitDistributionBM::hitTimeProb);
            probTime = HitDistSolveWrapperSP(new 
                HitDistSolveWrapper(hitTimeProbTemp, -0.5));
        }

        // Call helper function
        return hitTimeSampleHelper(probTime, tLow, tHigh, u);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


double HitDistributionBM::hitTimeGivenHitSample(double u) const {
    static const string routine = "HitDistributionBM::hitTimeGivenHitSample";
    
    try {
        // Brackets
        double tLow = 0.0;
        double tHigh = length;
        
        if(!probCondTime) {    
            HitDistWrapper hitTimeGivenHitProbTemp(this, &HitDistributionBM::hitTimeGivenHitProb);
            probCondTime = HitDistSolveWrapperSP(new
                HitDistSolveWrapper(hitTimeGivenHitProbTemp, -0.5));
        }
        
        // Call helper function
        return hitTimeSampleHelper(probCondTime, tLow, tHigh, u);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void HitDistributionBM::updateBridge(double valueStart, double valueEnd) {
    this->valueStart = valueStart;
    this->valueEnd   = valueEnd;

    // Recompute transient fields
    validatePop2Object();
}


#if 0
IObject* HitDistributionBM::clone() const {
    static const string routine = "HitDistributionBM::clone";
    try {
        IObject* Copy = CObject::clone(); // call parent
        HitDistributionBM* copy = dynamic_cast<HitDistributionBM*>(Copy);
        if(!copy) {
                throw ModelException("Clone method failed");
        }
        copy->probTime = probTime;
        copy->probCondTime = probCondTime;

        return copy;
    } catch(exception& e) {
            throw ModelException(e, routine);
    }
}
#endif


bool HitDistributionBM::hitEndpoint(bool start) const {
    double value = start ? valueStart : valueEnd;
    double barrier = start ? barrierStart : barrierEnd;

    // Up barrier and barrier <= value OR
    // Dn barrier and barrier >= value
    /*bool hit =  isUp && !Maths::isPositive(barrier - value) ||
               !isUp && !Maths::isPositive(value - barrier);*/
    bool hit =  isUp && value >= barrier ||
               !isUp && value <= barrier;

    return hit;
}


void HitDistributionBM::validatePop2Object() {
    static const string routine = "HitDistributionBM::validatePop2Object";

    // Check for negative vol and length of bridge
    if(Maths::isNegative(vol)) {
        throw ModelException(routine, "Volatility must be non-negative");
    }

    if(!Maths::isPositive(length)) {
        throw ModelException(routine, "The length of the BB must be positive");
    }

    // Floor vol and length to zero
    if(Maths::isZero(vol)) {
        vol = 0.0;
    }

    if(Maths::isZero(length)) {
        length = 0.0;
    }
    
    // Deduce if we are dealing with zero vol or time
    zeroVol = Maths::isZero(vol);
    zeroLength = Maths::isZero(length);

    // Update all transient fields that can be affected by moving the endpoints
    update();
}


void HitDistributionBM::update() {
    // Update transient fields when endpoints of BB change
    
    // Deduce instant hit
    hitStart = hitEndpoint(true);
    hitEnd   = hitEndpoint(false);
    
    // Compute hitting probability
    if( hitStart || hitEnd) {
        hitProb = 1.0;   // Hit immediately or by the end
    } else if (zeroLength || zeroVol) {
        hitProb = 0.0;   // Cannot hit with zero vol or length if not hit at edges
    } else {
        // The hitting probability for a Brownian Bridge
        // The two () are of the same sign and the exponential negative
        double prob = exp( -2.0 * (barrierStart - valueStart) * 
                                  (barrierEnd - valueEnd) / (Maths::square(vol) * length));
        // hitProb = Maths::max(PROBABILITY_THRESHOLD, prob);
        hitProb = prob;
    }

    // Common numerator
    commonNum = length * (barrierStart - valueStart);

    // Numerator elements
    num1 = valueEnd + barrierStart - valueStart - barrierEnd;
    num2 = barrierStart + barrierEnd - valueStart - valueEnd;
}


void HitDistributionBM::load(CClassSP& clazz) {
    REGISTER(HitDistributionBM, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(HitSample);
    EMPTY_SHELL_METHOD(defaultHitDistributionBM);
    // Fields
    FIELD(valueStart,    "Value of BM at start");
    FIELD(valueEnd,      "Value of BM at end");
    FIELD(barrierStart,  "Barrier at start");
    FIELD(barrierEnd,    "Barrier at end");
    FIELD(vol,           "Annualized implied vol");
    FIELD(length,        "Length of BB");
    FIELD(isUp,          "Flag for up barrier");
    // Transient fields
    FIELD(hitProb,       "Hitting probability within interval");
    FIELD_MAKE_TRANSIENT(hitProb);
    FIELD(commonNum,     "Common numerator");
    FIELD_MAKE_TRANSIENT(commonNum);
    FIELD(num1,          "Numerator 1");
    FIELD_MAKE_TRANSIENT(num1);
    FIELD(num2,          "Numerator 2");
    FIELD_MAKE_TRANSIENT(num2);
    FIELD(hitStart,      "Hit at start flag");
    FIELD_MAKE_TRANSIENT(hitStart);
    FIELD(hitEnd,        "Hit at end flag");
    FIELD_MAKE_TRANSIENT(hitEnd);
    FIELD(zeroVol,       "Zero volatility flag");
    FIELD_MAKE_TRANSIENT(zeroVol);
    FIELD(zeroLength,    "Zero length");
    FIELD_MAKE_TRANSIENT(zeroLength);
    clazz->setPublic();
}


IObject* HitDistributionBM::defaultHitDistributionBM() {
    return new HitDistributionBM();
}


double HitDistributionBM::hitTimeSampleHelper(HitDistSolveWrapperSP& wrapper,
                                        double tLow, double tHigh, double u) const {
    static const string routine = "HitDistributionBM::hitTimeSampleHelper";
    
    // Hit start: throw error
    if(hitProb == 0.0) {
        return HUGE_HITTING_TIME;
    }
    
    if(hitStart) {
        throw ModelException(routine, "Barrier has already been breached");
    }
    
    if(zeroVol) {
        if(hitEnd) {
            // Hitend return linear time
            double linearWeight = (valueEnd - barrierEnd) / num1;
            double linearTime = (1.0 - linearWeight) * length;
            return linearTime;        
        } else {
            // No hit
            return HUGE_HITTING_TIME;
        }
    }
    
    // Zero trading time means never hit
    if(zeroLength) {
        return HUGE_HITTING_TIME;
    }

    // Now deal with non-degenerate case
    
    // Set target
    wrapper->resetShift(-u);
    
    // Bracket and solve with Brent
    double solution = ZBrent_solve(*wrapper, tLow, tHigh, INVERSION_TOLERANCE);
    return solution;
}


CClassConstSP const HitDistributionBM::TYPE = CClass::registerClassLoadMethod(
    "HitDistributionBM", typeid(HitDistributionBM), HitDistributionBM::load);


///////////////////////////////////////////////////////////

/** ADDIN method for computing hitting distributions */
class HitDistributionBMAddin: public CObject{
    static CClassConstSP const TYPE;

    HitDistributionBMSP hitDist;        //!< Hitting distribution object
    DoubleArray         maturities;     //!< Requested maturities

    /** Computes the hitting time probabilities for the BB for various maturities */
    static IObjectSP hitTimeProb(HitDistributionBMAddin* params) {
        static const string routine = "HitDistributionBMAddin::hitTimeProb";
        
        try {
            int size = params->maturities.size();
            DoubleArraySP result(new DoubleArray(size));
            for(int i = 0; i < size; i++) {
                (*result)[i] = params->hitDist->hitTimeProb(params->maturities[i]);
            }
            
            return result;

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** Computes the hitting time probabilities for the BB for various maturities 
        conditional on a hit value */
    static IObjectSP hitTimeGivenHitProb(HitDistributionBMAddin* params) {
        static const string routine = "HitDistributionBMAddin::hitTimeGivenHitProb";
        try {
            int size = params->maturities.size();
            DoubleArraySP result(new DoubleArray(size));
            for(int i = 0; i < size; i++) {
                (*result)[i] = params->hitDist->hitTimeGivenHitProb(params->maturities[i]);
            }
            
            return result;

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    HitDistributionBMAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(HitDistributionBMAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultHitDistributionBMAddin);
        FIELD(hitDist,           "Hitting distribution");
        FIELD(maturities, "Maturities");

        Addin::registerClassObjectMethod("HITTING_TIME_PROB_BB",
                                         Addin::UTILITIES,
                                         "Hitting time distribution for BB",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)hitTimeProb);

        Addin::registerClassObjectMethod("HITTING_TIME_GIVENHIT_PROB_BB",
                                         Addin::UTILITIES,
                                         "Hitting time distribution for BB given a hit",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)hitTimeGivenHitProb);
    }

    static IObject* defaultHitDistributionBMAddin(){
        return new HitDistributionBMAddin();
    }
};

CClassConstSP const HitDistributionBMAddin::TYPE = CClass::registerClassLoadMethod(
"HitDistributionBMAddin", typeid(HitDistributionBMAddin), HitDistributionBMAddin::load);

///////////////////////////////////////////////////////////

/** ADDIN method for sampling hitting times */
class HitSampleBMAddin: public CObject{
    static CClassConstSP const TYPE;

    HitDistributionBMConstSP hitDist;        //!< HitSample object
    DoubleArray              uniforms;       //!< Uniform random numbers

    /** Samples the hitting time for the BB */
    static IObjectSP hitTimeSample(HitSampleBMAddin* params) {
        static const string routine = "HitSampleBMAddin::hitTimeSample";
        try {
            int size = params->uniforms.size();
            DoubleArraySP result(new DoubleArray(size));
            for(int i = 0; i < size; i++) {
                (*result)[i] = params->hitDist->hitTimeSample(params->uniforms[i]);
            }
            
            return result;

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** Samples the hitting time for the BB conditional on a hit value */
    static IObjectSP hitTimeGivenHitSample(HitSampleBMAddin* params) {
        static const string routine = "HitSampleBMAddin::hitTimeGivenHitSample";
        try {
            int size = params->uniforms.size();
            DoubleArraySP result(new DoubleArray(size));
            for(int i = 0; i < size; i++) {
                (*result)[i] = params->hitDist->hitTimeGivenHitSample(params->uniforms[i]);
            }
            
            return result;

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** Samples the hitting time for the BB conditional on a hit value */
    static IObjectSP hitNoHitSample(HitSampleBMAddin* params) {
        static const string routine = "HitSampleBMAddin::hitNoHitSample";
        try {
            int size = params->uniforms.size();
            DoubleArraySP result(new DoubleArray(size));
            for(int i = 0; i < size; i++) {
                (*result)[i] = params->hitDist->hitNoHitSample(params->uniforms[i]);
            }
            
            return result;

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** Basic validation */
    virtual void validatePop2Object() {
        static const string routine = "HitSampleBMAddin::validatePop2Object";
        
        int size = uniforms.size();
        for(int i = 0; i < size; i++) {
            if(Maths::isNegative(uniforms[i]) || 
               Maths::isPositive(uniforms[i] - 1.000001)) {
                throw ModelException(routine, "Uniform random numbers must be"
                                              " between 0.0 and 1.0");
            }
        }
    }


    /** for reflection */
    HitSampleBMAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(HitSampleBMAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultHitSampleBMAddin);
        FIELD(hitDist,           "Hitting distribution");
        FIELD(uniforms, "Uniform random numbers");

        Addin::registerClassObjectMethod("HITTING_TIME_SAMPLE_BB",
                                         Addin::UTILITIES,
                                         "Hitting time sample for BB",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)hitTimeSample);

        Addin::registerClassObjectMethod("HITTING_TIME_GIVENHIT_SAMPLE_BB",
                                         Addin::UTILITIES,
                                         "Hitting time sample for BB given a hit",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)hitTimeGivenHitSample);

        Addin::registerClassObjectMethod("HIT_NOHIT_SAMPLE_BB",
                                         Addin::UTILITIES,
                                         "Hit/NoHit sample for BB",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)hitNoHitSample);
    }

    static IObject* defaultHitSampleBMAddin(){
        return new HitSampleBMAddin();
    }
};

CClassConstSP const HitSampleBMAddin::TYPE = CClass::registerClassLoadMethod(
"HitSampleBMAddin", typeid(HitSampleBMAddin), HitSampleBMAddin::load);


////////////////////////////////////////////////////////////////////////


HitNoHitBB::HitNoHitBB(double barrierStart, double barrierEnd, 
                       double var, double length, bool isUp,
                       const DividendArray& divArray, const DoubleArray& divTimes):
valueStart(0.0), valueEnd(0.0), barrierStart(barrierStart), 
barrierEnd(barrierEnd), var(var), length(length), isUp(isUp),
hitProb(0.0) {
    validate();
    computeAdjustments(divArray, divTimes);
}


HitNoHitBB::HitNoHitBB(double var, double length, bool isUp,
                       const DividendArray& divArray, const DoubleArray& divTimes):
valueStart(0.0), valueEnd(0.0), barrierStart(0.0), 
barrierEnd(0.0), var(var), length(length), isUp(isUp), hitProb(0.0) {
    validate();
    computeAdjustments(divArray, divTimes);
}



void HitNoHitBB::updateBridge(double valueStart, double valueEnd) {
    this->valueStart = valueStart;
    this->valueEnd   = valueEnd;

    // Figure out if endpoints have hit
    if( ( isUp && (valueStart >= barrierStart || valueEnd >= barrierEnd ) ) ||
        (!isUp && (valueStart <= barrierStart || valueEnd <= barrierEnd ) ) ) {
        isHit = true;
    } else {
        isHit = false;
    }
    
    if(isHit) {
        // Hit immediately or by end
        hitProb = 1.0;
        return;
    } else if(zeroVar) {
        // No sure hit and zero var i.e. cannot hit
        hitProb = 0.0;
        return;
    } else {
        // The non degenerate hitting probability for a BB
        hitProb = exp( -2.0 * (barrierStart - valueStart) * 
                              (barrierEnd - valueEnd) / var);
        return;
    }
}


void HitNoHitBB::updateBridge(double barrierStart, double barrierEnd, 
                              double valueStart, double valueEnd) {
    // Update barrier and call update spot method
    this->barrierStart = barrierStart;
    this->barrierEnd = barrierEnd;
    updateBridge(valueStart, valueEnd);
}

int HitNoHitBB::updateAndSample(double valueStart, double valueEnd, double uniform) {
    updateBridge(barrierStart, barrierEnd, valueStart, valueEnd);
    if(areEndpointsHit()) {
        // For bridge that hits for sure, don't bother
        return 1;
    } else {
        updateBridgeOneDiv();
        return hitNoHitSample(uniform);
    }
}

int HitNoHitBB::updateAndSample(double barrierStart, double barrierEnd, double valueStart, double valueEnd, double uniform) {
    
    this->barrierStart = barrierStart;
    this->barrierEnd = barrierEnd;
    return updateAndSample(valueStart, valueEnd, uniform);
}

int HitNoHitBB::hitNoHitSample(double u) const {
    return (u <= hitProb) ? 1 : 0;
}


bool HitNoHitBB::areEndpointsHit() const {
    return isHit;
}


double HitNoHitBB::getHitProb() const {
    return hitProb;
}

#if 0

double HitNoHitBB::fwdBBDivAdjustment(const DividendArray& divArray, const DoubleArray& divTimes) const {
    // Compute dividend adjustment to be applied to initial condition of BB
    int nbDivs = divArray.size();
    if(!nbDivs) {
        return 0.0;
    }

//#if 0    
    double adjFwd   =  0.0;
    double tmpValueEnd = 2.0 * barrierEnd - valueEnd;
    for(int iDiv = 0; iDiv < nbDivs; iDiv++) {
        // Compute hitting time probability for BB starting at x and ending at (2.0 * barrierEnd - y > barrierEnd)
        double fwdDivWeight = 1.0 - hitTimeProb(valueStart, tmpValueEnd, 
            barrierStart, barrierEnd, sqrtVar, length, isUp, divTimes[iDiv]);
        adjFwd -= fwdDivWeight * log(1.0 - divArray[iDiv].getDivAmount());
    }
//#else
    double adjFwd   =  0.0;
    for(int iDiv = 0; iDiv < nbDivs; iDiv++) {
        adjFwd -=  (1.0 - divTimes[iDiv] / length) * log(1.0 - divArray[iDiv].getDivAmount());
    }
//#endif

    return adjFwd;
}


double HitNoHitBB::revBBDivAdjustment(const DividendArray& divArray, const DoubleArray& divTimes) const {
    // Compute dividend adjustment to be applied to initial condition of BB
    int nbDivs = divArray.size();
    if(!nbDivs) {
        return 0.0;
    }
//#if 0      
        double newValueStart   = 0.0;
        double newValueEnd     = valueStart - valueEnd;
        double newBarrierStart = barrierEnd - valueEnd;
        double newBarrierEnd   = barrierStart - valueEnd;

        double adjRev   = 0.0;
        double tmpValueEnd = 2.0 * newBarrierEnd - newValueEnd;
        for(int iDiv = 0; iDiv < nbDivs; iDiv++) {
	        // Compute hitting time probability for BB starting at 0.0 and ending at (2.0 * barrierEnd - y > barrierEnd)
	        double revDivWeight = 1.0 - hitTimeProb(newValueStart, tmpValueEnd,
		        newBarrierStart, newBarrierEnd, sqrtVar, length, isUp, length - divTimes[iDiv]);
	        adjRev += revDivWeight * log(1.0 - divArray[iDiv].getDivAmount());
        }
//#else
        double adjRev = 0.0;
        for(int iDiv = 0; iDiv < nbDivs; iDiv++) {
            adjRev +=  (divTimes[iDiv] / length) * log(1.0 - divArray[iDiv].getDivAmount());
        }
//#endif

    return adjRev;
}

double HitNoHitBB::hitTimeProb(double x, double y, double bStart, double bEnd, 
                               double sqrtVar, double length, bool isUp, double u) {
    // Trivial BB
    if(Maths::isZero(sqrtVar) || Maths::isZero(length)) {
        return 1.0;
    }
    
    // Boundary conditions
    if(u <= DBL_EPSILON) {
        return 0.0;
    } else if(u >= length - DBL_EPSILON) {
        return 1.0;
    }
    
    // The interior solution
    double denom = sqrt(u * (length - u)) * sqrtVar;
    double commonNum = length * (bStart - x);
    double w    = ( u * (y + bStart - x - bEnd) - commonNum ) / denom;
    double wbar = ( u * (bStart + bEnd - x - y) - commonNum ) / denom;
    if(!isUp) {
        w    = -w;
        wbar = -wbar;
    }
    double prob = N1(w) + N1(wbar);
    return prob;
}

#endif

void HitNoHitBB::validate() {
    static const string routine = "HitNoHitBB::validate";

    // Check for negative vol and length of bridge
    if(Maths::isNegative(var)) {
        throw ModelException(routine, "Variance must be non-negative");
    }

    // Floor vol and length to zero
    if(Maths::isZero(var)) {
        var = 0.0;
        sqrtVar = 0.0;
        zeroVar = true;
    } else {
        sqrtVar = sqrt(var);
        zeroVar = false;
    }

    if(Maths::isZero(length)) {
        length = 0.0;
    }
}

void HitNoHitBB::updateBridgeOneDiv() {

    double adjvalueStart = valueStart - forAdj;
    double adjvalueEnd = valueEnd + revAdj;
    double OneDivTime = CTime;
    double OneDivAmount = intAdj;

    if( ( isUp && (adjvalueStart >= barrierStart || adjvalueEnd >= barrierEnd ) ) ||
        (!isUp && (adjvalueStart <= barrierStart || adjvalueEnd <= barrierEnd ) ) ) {
        isHit = true;
    } else {
        isHit = false;
    }
    
    if(isHit) {
        // Hit immediately or by end
        hitProb = 1.0;
        return;
    } else if(zeroVar) {
        // No sure hit and zero var i.e. cannot hit
        hitProb = 0.0;
        return;
    } else {
        if (Maths::isPositive(OneDivAmount)) {

            // The non degenerate hitting probability for a BB with one dividend
            double bbvar = sqrt(var*OneDivTime*(length - OneDivTime));
            if (Maths::isZero(bbvar)) hitProb = 0.0;
            else {
                double argp0 = (length-OneDivTime)*(adjvalueStart-barrierStart)+OneDivTime*(adjvalueEnd+OneDivAmount-barrierEnd);
                double argp1 = (length-OneDivTime)*(adjvalueStart-barrierStart)+OneDivTime*(-adjvalueEnd-OneDivAmount+barrierEnd);
                double argp2 = (length-OneDivTime)*(-adjvalueStart+barrierStart)+OneDivTime*(adjvalueEnd-OneDivAmount-barrierEnd);
                double argp12 = (length-OneDivTime)*(adjvalueStart-barrierStart)+OneDivTime*(adjvalueEnd-OneDivAmount-barrierEnd);
                double p0 = N1(argp0/bbvar);
                double p1 = exp(-2.0*(adjvalueStart-barrierStart)*(OneDivAmount+adjvalueEnd-barrierEnd)/var)*N1(argp1/bbvar);
                double p2 = exp(-2.0*(adjvalueStart-barrierStart-OneDivAmount)*(adjvalueEnd-barrierEnd)/var)*N1(argp2/bbvar);
                double p12 = exp(-2.0*OneDivAmount*(adjvalueStart-adjvalueEnd-barrierStart+barrierEnd)/var)*N1(argp12/bbvar);
                hitProb = p0+p1+p2-p12;
            }
            return;
        }
        else {
            // The non degenerate hitting probability for a BB with no dividend
            hitProb = exp( -2.0 * (barrierStart - adjvalueStart)*(barrierEnd - adjvalueEnd) / var);
            return;
        }
    }
}

void HitNoHitBB::computeAdjustments(const DividendArray& divArray, const DoubleArray& divTimes) {

    forAdj = 0.0;
    revAdj = 0.0;
    intAdj = 0.0;
    CTime = 0.5*length;
    
    int nbDivs = divArray.size();
    if(!nbDivs) {
        return;
    }    
    
    double Sum = 0.0;
    double CT = 0.0;
    int iDiv;
    for(iDiv = 0; iDiv < nbDivs; iDiv++) {
        CT +=  -log(1.0 - divArray[iDiv].getDivAmount())*computeWeight(divTimes[iDiv],CTime,2)*divTimes[iDiv];
        Sum += -log(1.0 - divArray[iDiv].getDivAmount())*computeWeight(divTimes[iDiv],CTime,2);
    }
    if (Maths::isPositive(Sum)) CTime = CT/Sum;
    
    double fAdj = 0.0;
    double rAdj = 0.0;
    double iAdj = 0.0;
    for(iDiv = 0; iDiv < nbDivs; iDiv++) {
        fAdj +=  -log(1.0 - divArray[iDiv].getDivAmount())*computeWeight(divTimes[iDiv],CTime,0);
        rAdj +=  -log(1.0 - divArray[iDiv].getDivAmount())*computeWeight(divTimes[iDiv],CTime,1);
        iAdj +=  -log(1.0 - divArray[iDiv].getDivAmount())*computeWeight(divTimes[iDiv],CTime,2);
    }
    
    forAdj = fAdj;
    revAdj = rAdj;
    intAdj = iAdj;    
    return;
}

double HitNoHitBB::computeWeight(double t, double CT, int Wtype){
    
    double fWeight, rWeight, iWeight;

    if (Maths::isZero(length)) return 0;

    if ((Maths::isZero(CT))||(Maths::isZero(length-CT))){
        if (Maths::isZero(CT)){
            fWeight = 1.0;
            rWeight = 0.0;
            iWeight = 0.0;  
        }
        if (Maths::isZero(length-CT)){
            fWeight = 0.0;
            rWeight = 1.0;
            iWeight = 0.0;  
        }
    } else {
        if (t<=CT) {
            fWeight = Maths::square(1-t/CT)*(3-2*(1-t/CT));
            rWeight = 0.0;
            iWeight = 1.0-fWeight;
        } else {
            fWeight = 0.0;
            rWeight = Maths::square((t-CT)/(length-CT))*(3-2*(t-CT)/(length-CT));
            iWeight = 1.0-rWeight;
        }
    }

    if(Wtype==0) {
        return fWeight;
    } else if(Wtype==1) {
        return rWeight;
    } else {
        return iWeight;
    }
}

DRLIB_END_NAMESPACE
