//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolParamTweak.cpp
//
//   Description : Interface for Tweaked Parametrized Vol Surface
//
//   Author      : Jean-Noël Juston
//
//   Date        : 25 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/VolParamTweak.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

class TweakLattice: public CObject,
                    virtual public ITweakableWithRespectTo<VolParallel>,
                    virtual public ITweakableWithRespectTo<VolPointwise>,
                    virtual public VegaMatrix::IShift,
                    virtual public RootTimeVega::IShift,
                    virtual public VegaSkewParallel::IShift,
                    virtual public VegaSkewPointwise::IShift,
                    virtual public VolLevel::Shift,
                    virtual public VolParallelShift::Shift,
                    virtual public VolBenchmarkShift::Shift,
                    virtual public PowerVega::Shift,
                    virtual public VolRelativeShift::IShift,
                    virtual public VolAbsoluteShift::IShift {
//private:
public:
    friend class Pointwise;

    const CLatticeDouble*    strikes; // $unregistered
    const DateTimeArray*     maturities; // $unregistered
    const CVolParamTweak*    volTweak; // $unregistered
    const CVolBase*          vol; // $unregistered
    CLatticeDouble*          lattice; // $unregistered
    const TimeMetric*        timeMetric; // $unregistered
    static const string name; // = "TweakedVol"
public:
    static CClassConstSP const TYPE;

    TweakLattice(const CLatticeDouble&    strikes,
                 const DateTimeArray&     maturities,
                 const CVolParamTweak*    volTweak,
                 const CVolBase*          vol,                
                 CLatticeDouble&          impV,
                 const TimeMetric*        timeMetric): 
        CObject(TYPE), strikes(&strikes), 
        maturities(&maturities), volTweak(volTweak), vol(vol), lattice(&impV),
        timeMetric(timeMetric) {}

    void negVolMessage(const string& method, double shift) {
        throw ModelException(method,
                             "shift of " + Format::toString(shift) +
                             " makes vol -ve for " + vol->getName());
    }

    /** Implementation of ITweakableWithRespectTo<VolParallel>. Note not used in the normal
        manner but rather tweaks are applied at pricing time ie grid of
        implied vols is calculated and the shift is then applied to that */
    virtual string sensName(const VolParallel*) const{
        return name; // this method is never called
    }
    // actually does the work
    virtual TweakOutcome sensShift(const PropertyTweak<VolParallel>& tweak){
        // easy - just add shift amount to all implied vols
        bool checkForNegVol = tweak.coefficient < 0.0;
        for (int i = 0; i < lattice->size(); i++){
            for (int j = 0; j < (*lattice)[i].size(); j++){
                (*lattice)[i][j] += tweak.coefficient;
                if (checkForNegVol && Maths::isNegative((*lattice)[i][j])){
                    negVolMessage("TweakLattice::sensShift", tweak.coefficient);
                }
            }
        }
        return TweakOutcome(tweak.coefficient, false);
    }
    //// same thing for vega skew
    virtual string sensName(VegaSkewParallel* shift) const{
        return name; // this method is never called
    }
    //// actually does the work
    virtual bool sensShift(VegaSkewParallel* shift){
        for (int i = 0; i < lattice->size(); i++){
            for (int j = 0; j < (*lattice)[i].size(); j++){
                // rely on skewShift() to check for 0 strike
                double volShift = shift->skewShift((*strikes)[i][j]);
                (*lattice)[i][j] += volShift;
                if (Maths::isNegative((*lattice)[i][j])){
                    (*lattice)[i][j] = 0.0;
                    // too tough negVolMessage("TweakLattice::sensShift", volShift);
                }
            }
        }
        return false;
    }


    // The vol benchmark greeks/scenarios are all based around the same
    // complex loop of code determining the variance to be allocated from
    // shifting the benchmark to the time we want the vol for.
    // The portion that differs is getting the required shift size and
    // how the param vol gets adjusted (e.g. scenarios tend to have a
    // vol floor).
    // Rather than have the s ame code 3 times, each benchmark shift 
    // delegates the body of the code to benchmarkVolShifter and gives
    // it an appropriate BMarkVolShift which implements the differing
    // shift size and adjustment techniques

    // Base class for the shift size/adjust functions
    class BMarkVolShift {
    public:
        // get shift size for [i,j]th point on lattice
        virtual double shiftSize(TweakLattice* tl,
                                 int           i, 
                                 int           j) = 0;

        // add volToAdd to [i,j]th point on lattice
        virtual void addVol(TweakLattice* tl,
                            double        volToAdd,
                            int           i, 
                            int           j) = 0;
    };

    // For VegaPointwise
    class Pointwise: public BMarkVolShift {
    public:
        Pointwise(double coeff) : coeff(coeff) {}

        virtual double shiftSize(TweakLattice* tl, int i, int j) {
            return coeff;
        }

        virtual void addVol(TweakLattice* tl,
                            double        volToAdd,
                            int           i, 
                            int           j) {        
            (*tl->lattice)[i][j] += volToAdd;  // just add it                                     
        }
            
    private:
        double coeff;
    };

    // For VegaSkewPointwise
    class SkewPointwise: public BMarkVolShift {
    public:
        SkewPointwise(const VegaSkewPointwise* pwise) : pwise(pwise) {}

        virtual double shiftSize(TweakLattice* tl, int i, int j) {
            return pwise->skewShift((*tl->strikes)[i][j]);
        }

        virtual void addVol(TweakLattice* tl,
                            double        volToAdd,
                            int           i, 
                            int           j) {
            (*tl->lattice)[i][j] += volToAdd;  // just add it
            (*tl->lattice)[i][j] = Maths::max((*tl->lattice)[i][j], 0.0);
        }
            
    private:
        const VegaSkewPointwise* pwise;
    };

    // For VolBenchmarkShift
    class BenchmarkShift: public BMarkVolShift {
    public:
        BenchmarkShift(const VolBenchmarkShift* bm) : bm(bm) {}

        virtual double shiftSize(TweakLattice* tl, int i, int j) {
            return bm->getShiftSize();
        }

        virtual void addVol(TweakLattice* tl,
                            double        volToAdd,
                            int           i, 
                            int           j) {
            double shiftVol = (*tl->lattice)[i][j] + volToAdd;
            if( Maths::isNegative(volToAdd) ) {
                // Floor a downshift, but don't increase a low vol
                if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                    (*tl->lattice)[i][j]  = 
                        Maths::min((*tl->lattice)[i][j], 
                                   VolParallelShift::MIN_SPOT_VOL);
                }
                else {
                    (*tl->lattice)[i][j] = shiftVol;
                }                
            } else {
                // Don't floor an upshift
                (*tl->lattice)[i][j] = shiftVol;
            }
        }
            
    private:
        const VolBenchmarkShift* bm;
    };

    // For VegaMatrix
    class Matrix: public BMarkVolShift {
    public:
        Matrix(const VegaMatrix* vm) : vm(vm) {}

        virtual double shiftSize(TweakLattice* tl, int i, int j) {
            return vm->getShiftSize();
        }

        virtual void addVol(TweakLattice* tl,
                            double        volToAdd,
                            int           i, 
                            int           j) {   
            // having worked out our shift, decide if we actually want it

            // what strike are we trying to shift ?
            const DoubleArray& strikes = *(vm->getXAxisValues());
            int    idx = vm->getUpperIdx();

            double k = (*tl->strikes)[i][j];
 
            const double alpha   = 0.05;  // attenuation
            const double loBound = 0.9;
            const double hiBound = 1.1;

            double loK   = idx == 0 ? loBound*strikes[0] : strikes[idx-1];
            double thisK = strikes[idx];
            double hiK   = idx == strikes.size()-1 ? 
                hiBound*strikes[strikes.size()-1] : strikes[idx+1];
           
            double g;
            if (k <= thisK) {
                g = exp(log(alpha)*Maths::square((thisK-k)/(thisK-loK)));
            }
            else {
                g = exp(log(alpha)*Maths::square((k-thisK)/(hiK-thisK)));
            }            

            (*tl->lattice)[i][j] += volToAdd * g;
        }
            
    private:
        const VegaMatrix* vm;
    };

    bool benchmarkVolShifter(VectorShift* shift, BMarkVolShift* calculate){
        return benchmarkVolShifter(shift->getExpiry(), shift->getExpiries(),
                                   calculate);
    }

    bool benchmarkVolShifter(ExpiryConstSP expiry, ExpiryArrayConstSP expiries,
                             BMarkVolShift *calculate) {
        int e = expiry->search(expiries.get());

        return benchmarkVolShifter(
            ExpiryWindow(e > 0 ? (*expiries)[e-1] : ExpirySP(),
                         (*expiries)[e],
                         e < expiries->size()-1 ? (*expiries)[e+1] : ExpirySP()),
            calculate);
    }

    bool benchmarkVolShifter(const ExpiryWindow& expiries, BMarkVolShift* calculate){
        try {
            // turn expiries into dates
            const DateTime& valueDate = volTweak->valueDate;
            DateTime mid = expiries.expiry->toDate(valueDate);
            DateTime lo = !expiries.previous ?
                              valueDate : expiries.previous->toDate(valueDate);
            // flawed for vols > last benchmark date. Ignore for now
            DateTime hi = !expiries.next ? mid :
                                            expiries.next->toDate(valueDate);

            // calc trading time to three benchmarks
            double ttLo = timeMetric->yearFrac(valueDate, lo);
            double ttMid = timeMetric->yearFrac(valueDate, mid);
            double ttHi = timeMetric->yearFrac(valueDate, hi);

            // calculate vol at current benchmark
            CSliceDouble   slice(1);
            DateTimeArray  mat(1, mid);
            for (int i = 0; i < lattice->size(); i++){
                if ((*maturities)[i].isGreater(lo) && 
                    (!!expiries.previous || (*maturities)[i].isGreater(mid)) &&
                    hi.isGreater((*maturities)[i])){
                    // lies between lo and hi (excludes case of 1st bm)

                    // what's going on here then ?
                    // we want the additional vol to add to the implied vol
                    // (baseVol here) at the time we're interested in (wanted)
                    // such that the total variance is equal to the unshifted
                    // variance plus a linearly interpolated portion of extra
                    // variance that comes from bumping the benchmark vol

                    // the untweaked benchmark vol
                    CSliceDouble   volsAtBM((*strikes)[i].size());
                    volTweak->paramVol->ComputeImpVol(vol,
                                                      (*strikes)[i], 
                                                      mat, volsAtBM);

                    // the untweaked vol at the time we're interested in
                    CSliceDouble   baseVol((*strikes)[i].size());
                    DateTimeArray  wanted(1, (*maturities)[i]);

                    volTweak->paramVol->ComputeImpVol(vol,
                                                      (*strikes)[i], 
                                                      wanted, baseVol);


                    double tt = timeMetric->yearFrac(valueDate, 
                                                     (*maturities)[i]);
                    for (int j = 0; j < (*lattice)[i].size(); j++){
                        double shiftSize = calculate->shiftSize(this, i, j);

                        /* the untweaked variance at the time we're
                           interested in */
                        double baseVar = tt * baseVol[j] * baseVol[j];

                        // the extra variance at the benchmark due to tweaking
                        double extraVar = (2.0 * volsAtBM[j] + 
                                           shiftSize)*shiftSize*ttMid;
                        bool afterMid = (*maturities)[i].isGreater(mid);
                        double varToAdd;
                    
                        // linearly apportion extra variance to 
                        // time we're interested in
                        if (afterMid) {
                            varToAdd = extraVar*((ttHi - tt)/(ttHi - ttMid));
                        }
                        else {
                            varToAdd = extraVar*((tt - ttLo)/(ttMid - ttLo));
                        }

                        // back out equivalent vol shift
                        double volToAdd = sqrt((baseVar+varToAdd)/tt) - 
                            baseVol[j];

                        //and add it on
                        calculate->addVol(this, volToAdd, i, j);
                    }
                }
                else if ((!expiries.next && 
                         (*maturities)[i].isGreaterOrEqual(hi)) ||
                         (!expiries.previous && hi.isGreater((*maturities)[i]))) {
                    //                    cerr << "case 2 " << i << "\n";

                    // we're tweaking the last benchmark and our vols are
                    // >= the last benchmark OR we're in the interval between
                    // today and the first benchmark.
                    // So just add the vol shift
                    for (int j = 0; j < (*lattice)[i].size(); j++){
                        double shiftSize = calculate->shiftSize(this, i, j);
                        calculate->addVol(this, shiftSize, i, j);
                    }                    
                }
            }

//             {for (int i = 0; i < lattice->size(); ++i)
//                     for (int j = 0; j < (*lattice)[i].size(); ++j)
//                         cerr << (*lattice)[i][j] << " ";
//                 cerr << "\n";}

            return false;
        }
        catch (exception& e) {
            throw ModelException(e, "TweakLattice::benchmarkVolShifter",
                                 "Failed to shift " + expiries.expiry->toString() +
                                 " point on " + vol->getName());
        }
    }

    //// same thing for vega pointwise
    virtual string sensName(const VolPointwise*) const{
        return name;
    }
    virtual ExpiryWindowArrayConstSP sensQualifiers(const VolPointwise*) const{
        return ExpiryWindowArrayConstSP(); // this method is never called
    }
    TweakOutcome sensShift(const PropertyTweak<VolPointwise>& shift){
        Pointwise calculate(shift.coefficient);
        return TweakOutcome(shift.coefficient,
                            benchmarkVolShifter(*shift.qualifier, &calculate));
    }

    /** VegaMatrix Interface */
    string sensName(VegaMatrix* shift) const{
        return name;
    }
    ExpiryArrayConstSP sensExpiries(VegaMatrix* shift) const{
        return ExpiryArrayConstSP(); // this method is never called
    }
    bool sensShift(VegaMatrix* shift){
        Matrix calculate(shift);
        // spoof off a pointwise vega to reuse benchmarkVolShifter
        // for the time dependent shift. Matrix supplies the "are we
        // affected by this strike" logic
        return benchmarkVolShifter((*shift->getExpiries())[shift->getExpiryIdx()],
                                   shift->getExpiries(),
                                   &calculate);
    }
    /** RootTimeVega Interface */
    string sensName(RootTimeVega* shift) const{
        return name;
    }

    bool sensShift(RootTimeVega* shift){
        double shiftSize = shift->getShiftSize();
        bool checkForNegVol = shiftSize < 0.0;
        const DateTime& valueDate = volTweak->valueDate;
        for (int i = 0; i < lattice->size(); i++){
            double tt = timeMetric->yearFrac(valueDate, (*maturities)[i]);
            double volShift = shift->rtVegaShift(tt);
            for (int j = 0; j < (*lattice)[i].size(); j++){
                (*lattice)[i][j] += volShift;
                if (checkForNegVol && Maths::isNegative((*lattice)[i][j])){
                    negVolMessage("TweakLattice::sensShift", shiftSize);
                }
            }
        }
        return false; // none of our components has a vega type sensitivity
    }

    /** PowerVega Interface */
    string sensName(PowerVega* shift) const{
        return name;
    }

    bool sensShift(PowerVega* shift){
        const DateTime& valueDate = volTweak->valueDate;
        for (int i = 0; i < lattice->size(); i++){
            double tt = timeMetric->yearFrac(valueDate, (*maturities)[i]);
            for (int j = 0; j < (*lattice)[i].size(); j++){
                double shiftSize = shift->powerShift(tt);
                double shiftVol = (*lattice)[i][j] + shiftSize;
                if( Maths::isNegative(shiftSize) ) {
                    // Floor a downshift, but don't increase a low vol
                    if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                        (*lattice)[i][j]  = 
                            Maths::min((*lattice)[i][j], 
                                       VolParallelShift::MIN_SPOT_VOL);
                    }
                    else {
                        (*lattice)[i][j] = shiftVol;
                    }                 
                } else {
                    // Don't floor an upshift
                    (*lattice)[i][j] = shiftVol;
                }
            }
        }
        return false; // none of our components has a vega type sensitivity
    }

    /** VegaSkewPointwise Interface */
    string sensName(VegaSkewPointwise* shift) const{
        return name;
    }
    ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const {
        return ExpiryArrayConstSP(); // this method is never called
    }
    bool sensShift(VegaSkewPointwise* shift){
        SkewPointwise calculate(shift);
        return benchmarkVolShifter(shift, &calculate);
    }

    /** VolLevel */
    string sensName(VolLevel* shift) const {
        return name;
    }
    bool sensShift(VolLevel* shift) {
        // easy - return the supplied vol level
        double level = shift->getShiftSize();
        if (Maths::isNegative(level)) {
            throw ModelException("TweakLattice::sensShift",
                                 "vol level (" + 
                                 Format::toString(level) + 
                                 ") is -ve for " + name);
        }
        
        for (int i = 0; i < lattice->size(); i++){
            for (int j = 0; j < (*lattice)[i].size(); j++){
                (*lattice)[i][j] = level;
            }
        }
        return false;
    }

    // VolParallelShift
    string sensName(VolParallelShift* shift) const {
        return name;
    }
    bool sensShift(VolParallelShift* shift) {
        // easy - just add to the vol
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {      
            for (int i = 0; i < lattice->size(); i++){
                for (int j = 0; j < (*lattice)[i].size(); j++){
                    double shiftVol = (*lattice)[i][j] + shiftSize;
                    if( Maths::isNegative(shiftSize) ) {
                        // Floor a downshift, but don't increase a low vol
                        if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                            (*lattice)[i][j]  = 
                                Maths::min((*lattice)[i][j], 
                                           VolParallelShift::MIN_SPOT_VOL);
                        }
                        else {
                            (*lattice)[i][j] = shiftVol;
                        }      
                    } else {
                        // Don't floor an upshift
                        (*lattice)[i][j] = shiftVol;
                    }
                }
            }
        }
        return false;
    }

    /** VolBenchmarkShift Interface */
    string sensName(VolBenchmarkShift* shift) const{
        return name;  // not actually called
    }

    bool sensShift(VolBenchmarkShift* shift){
        BenchmarkShift calculate(shift);
        return benchmarkVolShifter(shift->getExpiry(), shift->getAllExpiries(),
                                   &calculate);
    }

    /** VolRelativeShift Interface */
    string sensName(VolRelativeShift* shift) const{
        return name;
    }

    bool sensShift(VolRelativeShift* shift){
        const DateTime& valueDate = volTweak->valueDate;
        for (int i = 0; i < lattice->size(); i++){
            double shiftSize = shift->shiftSize(valueDate,(*maturities)[i]);

            for (int j = 0; j < (*lattice)[i].size(); j++){
                double shiftVol = (*lattice)[i][j] + shiftSize;
                if( Maths::isNegative(shiftSize) ) {
                    // Floor a downshift, but don't increase a low vol
                    if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                        (*lattice)[i][j]  = 
                            Maths::min((*lattice)[i][j], 
                                       VolParallelShift::MIN_SPOT_VOL);
                    }
                    else {
                        (*lattice)[i][j] = shiftVol;
                    }                 
                } else {
                    // Don't floor an upshift
                    (*lattice)[i][j] = shiftVol;
                }
            }
        }
        return false; // none of our components has a vega type sensitivity
    }

    /** VolAbsoluteShift Interface */
    string sensName(VolAbsoluteShift* shift) const{
        return name;
    }

    bool sensShift(VolAbsoluteShift* shift){
        const DateTime& valueDate = volTweak->valueDate;
        for (int i = 0; i < lattice->size(); i++){
            double shiftSize = shift->shiftSize(valueDate,(*maturities)[i]);

            for (int j = 0; j < (*lattice)[i].size(); j++){
                double shiftVol = (*lattice)[i][j] + shiftSize;
                if( Maths::isNegative(shiftSize) ) {
                    // Floor a downshift, but don't increase a low vol
                    if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                        (*lattice)[i][j]  = 
                            Maths::min((*lattice)[i][j], 
                                       VolParallelShift::MIN_SPOT_VOL);
                    }
                    else {
                        (*lattice)[i][j] = shiftVol;
                    }                 
                } else {
                    // Don't floor an upshift
                    (*lattice)[i][j] = shiftVol;
                }                 
            }
        }
        return false; // none of our components has a vega type sensitivity
    }



private:
    static void myLoad(CClassSP& clazz){
        REGISTER(TweakLattice, clazz);
        SUPERCLASS(CObject);
    }

};
//// do minimum needed to make this class derived from IObject
CClassConstSP const TweakLattice::TYPE = CClass::registerClassLoadMethod(
    "TweakLattice", typeid(TweakLattice), myLoad);
        
const string TweakLattice::name = "TweakedVol";

/** Interface for Tweaked Parametrized Vol Surface */
CClassConstSP const CVolParamTweak::TYPE = CClass::registerClassLoadMethod(
    "VolParamTweak", typeid(CVolParamTweak), load);

void CVolParamTweak::load(CClassSP& clazz){
    REGISTER(CVolParamTweak, clazz);
    SUPERCLASS(CVolParam);
}

CVolParamTweak::CVolParamTweak(const CClassConstSP& clazz) : CVolParam(clazz){}

//// constructor: takes clone of shift, but not of paramVol
CVolParamTweak::CVolParamTweak(const CVolParamSP&         paramVol,
                               const DateTime&            valueDate,
                               const TimeMetricConstSP&   timeMetric,
                               const ITweakID*            shift):
    CVolParam(TYPE), paramVol(paramVol), valueDate(valueDate),
    timeMetric(timeMetric){
    const IObject& obj = dynamic_cast<const IObject&>(*shift);
    shiftSP.reset(obj.clone());
    ASSERT(theShift = dynamic_cast<ITweakID*>(shiftSP.get()));
}

VolSurface* CVolParamTweak::spotVolSurfaceFromStrikes(
    const CVolBase*       vol,
    const CDoubleArray&   strikes) const{
    static const string routine("CVolParamTweak::spotVolSurfaceFromStrikes");
    try{
        VolSurfaceSP volToShift(paramVol->
                                spotVolSurfaceFromStrikes(vol, strikes));
        theShift->shift(volToShift);
        return volToShift.release();
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** Implementation of method defined in CVolParam. */
void CVolParamTweak::ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const{
    // compute unperturbed vols
    paramVol->ComputeImpVol(vol, strikes, maturities, impV);
    // set up our tweaker - it has a pointer to the strikes & lattice
    TweakLattice tweaker(strikes, maturities, this, vol, impV, 
                         timeMetric.get()); 
    IObjectSP objToShift(IObjectSP::attachToRef(&tweaker));
    // do the business
    try{
        // calls relevant shift method on TweakLattice for given shift
        theShift->shift(objToShift);
    } catch (exception& e){
        throw ModelException(e, "CVolParamTweak::ComputeImpVol", "Shift for "+
                             shiftSP->getClass()->getName()+" failed for "
                             "parameterised vol");
    }
}


DRLIB_END_NAMESPACE
