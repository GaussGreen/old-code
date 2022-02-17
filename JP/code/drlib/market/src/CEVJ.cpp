//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CEVJ.hpp
//
//   Description : CEV+jump volatility
//
//   Date        : 30 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CEVJ.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/CEVJProcessed.hpp"
#include "edginc/Maths.hpp"

static const double MAX_CEVPOWER = 0.999, MIN_CEVPOWER = 0.001;
static const double MAX_DIFF_VOL = 3.0, MIN_DIFF_VOL = 0.01;
static const double MAX_JUMP_WIDTH = 0.5, MIN_JUMP_WIDTH = 0.0;
static const double MAX_JUMP_MEAN = 0.3, MIN_JUMP_MEAN = -0.3;
static const double MAX_JUMP_RATE = 20.0, MIN_JUMP_RATE = 0.0;

DRLIB_BEGIN_NAMESPACE

class CEVJHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CEVJ, clazz);
        SUPERCLASS(CVolBase);
        IMPLEMENTS(VolLevel::Shift);
        IMPLEMENTS(VolParallelShift::Shift);
        IMPLEMENTS(ITweakableWithRespectTo<VolParallel>);
        IMPLEMENTS(ITweakableWithRespectTo<VolPointwise>);
        IMPLEMENTS(JumpRatePointwise::IShift);
        IMPLEMENTS(JumpRateParallel::IShift);
        IMPLEMENTS(CEVPowerPointwise::IShift);
        IMPLEMENTS(CEVPowerParallel::IShift);
        IMPLEMENTS(Theta::IShift);
        IMPLEMENTS(VolAbsoluteShift::IShift);
        IMPLEMENTS(VolBenchmarkShift::Shift);
        IMPLEMENTS(PowerVega::Shift);
        EMPTY_SHELL_METHOD(defaultCEVJ);
        FIELD(name, "Vol identifier");
        FIELD(baseDate, "Base Date");
        FIELD_MAKE_TRANSIENT(baseDate);
        FIELD(timeMetric, "Time metric");
        FIELD_MAKE_TRANSIENT(timeMetric);
        FIELD(ATMVolArr, "ATM vol values");
        FIELD(BenchMarkStrg, "String Array (1W, 1M, 1Y...) for BenchMarks")
        FIELD(ATMVolBM, "ATM vol bench mark");
        FIELD_MAKE_TRANSIENT(ATMVolBM);
        FIELD(CEVPowerArr, "Lambda values");
        FIELD(JumpMeanArr, "jump mean value");
        FIELD(JumpRateArr, "jumps per year");
        FIELD(JumpWidthArr, "jump width");
        FIELD(JumpMode, "0=proportional(default), else=absolute");
        FIELD(SpotRef, "vol spot reference");
        FIELD(VolMode, "0=log-normal(N/A->1), 1=CEV+jump, 2=Lambda+jump(N/A->1),");
        FIELD(DEBUG_fix2ndMoment, "0=ON, elxe=Off");
        FIELD_MAKE_OPTIONAL(DEBUG_fix2ndMoment);

        ClassSetAcceptMethod(CEVJ::acceptValueDateCollector);

        Addin::registerConstructor("CEVJ_VOL",
                                   Addin::MARKET,
                                   "Creates a handle to CEV+jump vol object",
                                   CEVJ::TYPE);
    }

    static IObject* defaultCEVJ(){
        return new CEVJ();
    }
};

CClassConstSP const CEVJ::TYPE = CClass::registerClassLoadMethod(
               "CEVJ", typeid(CEVJ), CEVJHelper::load);

/** Returns name of vol */
string CEVJ::getName() const{
    return name;
}

CEVJ::CEVJ(): CVolBase(TYPE){
    DEBUG_fix2ndMoment = 1;
    VolMode = 1;
    JumpMode = 0;
    SpotRef = 0.0;

}

/** populate from market cache */
void CEVJ::getMarket(const IModel* model, const MarketData* market) {
    /** returns a constant reference to surface to be used for the backbone */
    VolSurfaceSP backbone = VolSurfaceSP(VolSurfaceSP::dynamicCast(
                market->GetData(getName(),VolSurface::TYPE)));
    backbone->getMarket(model,market);  // put holiday and so on into volsurface.
    baseDate = backbone->getBaseDate();
    timeMetric = backbone->getTimeMetric();
}

void CEVJ::validatePop2Object(){
    const static string method = "CEVJ::validatePop2Object";
    try {
        int size = ATMVolArr.size();

        if (size == 0)
            throw ModelException(method, "no diffusion vol supplied.");

        // Make ATMVolBM
        //DateTime valDatePM(baseDate.getDate(), DateTime::END_OF_DAY_TIME);
        ATMVolBM = ExpiryArraySP(new ExpiryArray(size));
        for (int i=0;i<size;i++)
            (*ATMVolBM)[i] = ExpirySP(new MaturityTimePeriod(BenchMarkStrg[i], DateTime::END_OF_DAY_TIME));
        
        if (size != BenchMarkStrg.size())
            throw ModelException(method, "the size of BenchMarkStrg must be "
                                 "the same as ATMVolArr value array.");
        if (size != CEVPowerArr.size())
            throw ModelException(method, "the size of bench marks must be "
                                 "the same as CEVPower value array.");
        if (size != JumpMeanArr.size())
            throw ModelException(method, "the size of bench marks must be "
                                 "the same as JumpMean value array.");
        if (size != JumpRateArr.size())
            throw ModelException(method, "the size of bench marks must be "
                                 "the same as JumpRate value array.");
        if (size != JumpWidthArr.size())
            throw ModelException(method, "the size of bench marks must be "
                                 "the same as JumpWidth value array.");
        if (Maths::isZero(SpotRef))
            throw ModelException(method, "SpotRef = 0 is wrong");
    }
    catch (exception& e) {
        throw ModelException(e,method,"Failed for vol ("+getName()+")");
    }   
}

/** Combines market and instrument data together to give a Processed Vol */
CVolProcessed* CEVJ::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      asset) const{

    const static string method = "CEVJ::getProcessedVol";
    int size, i;

    size = ATMVolArr.size();
    // validate diffusion vol
    for (i=0; i < size; i++)
    {//ATMVols
        if (ATMVolArr[i] > MAX_DIFF_VOL || ATMVolArr[i] < MIN_DIFF_VOL)
        {
            throw ModelException(method, "DiffVol out of range : " + 
                                 Format::toString(MIN_DIFF_VOL) +
                                 " to " + Format::toString(MAX_DIFF_VOL));
        }
        if (CEVPowerArr[i] > MAX_CEVPOWER || CEVPowerArr[i] < MIN_CEVPOWER)
        {
            throw ModelException(method, "CEVPower out of range : " + 
                                 Format::toString(MIN_CEVPOWER) +
                                 " to " + Format::toString(MAX_CEVPOWER));
        }
    }

    // validate jump parameters
    for (i=0; i < size; i++)
    {
        //JumpRate Arr
        if (JumpRateArr[i] > MAX_JUMP_RATE || JumpRateArr[i] < MIN_JUMP_RATE)
        {
            throw ModelException(method, "jump rate out of range : " +
                                 Format::toString(MIN_JUMP_RATE) +
                                 " to " + Format::toString(MAX_JUMP_RATE));
        }
        //JumpMean Arr
        if (JumpMeanArr[i] > MAX_JUMP_MEAN || JumpMeanArr[i] < MIN_JUMP_MEAN)
        {
            throw ModelException(method, "jump mean out of range : " + 
                                 Format::toString(MIN_JUMP_MEAN) +
                                 " to " + Format::toString(MAX_JUMP_MEAN));
        }
        //JumpWidth Arr
        if (JumpWidthArr[i] > MAX_JUMP_WIDTH || 
            JumpWidthArr[i] < MIN_JUMP_WIDTH)
        {
            throw ModelException(method, "jump width out of range : " +
                                 Format::toString(MIN_JUMP_WIDTH) +
                                 " to " + Format::toString(MAX_JUMP_WIDTH));
        }
    }
    
    CEVJProcessed* ptr = new CEVJProcessed(dynamic_cast<CEVJ*>(clone()));
    ptr->createExpiryDates();
    return ptr;
}

/** Shifts the object using given shift (see VolLevel::Shift)*/
bool CEVJ::sensShift(VolLevel* shift) {
    static const string method = "CEVJ::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (Maths::isNegative(shiftSize)) {
            throw ModelException(method,
                                 "vol level (" + 
                                 Format::toString(shiftSize) + 
                                 ") is -ve");
        }

        for (int j = 0; j < ATMVolArr.size(); j++) {
            ATMVolArr[j] = shiftSize;
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             "VolLevel scenario failed for " + getName());
    }
    return false;
}

/** Shifts the object using given shift */
TweakOutcome CEVJ::sensShift(const PropertyTweak<VolParallel>& tweak)
{
    VolParallelShift tmpShift(tweak.coefficient);
    return TweakOutcome(tweak.coefficient, sensShift(&tmpShift));
}

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryWindowArrayConstSP CEVJ::sensQualifiers(const VolPointwise*) const{
    return ExpiryWindow::series(ATMVolBM);
}

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryArrayConstSP CEVJ::sensExpiries(JumpRatePointwise* shift) const{
    return ATMVolBM;
}

/** Returns the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this vol */
ExpiryArrayConstSP CEVJ::sensExpiries(CEVPowerPointwise* shift) const{
    return ATMVolBM;
}

/** Shifts the object using given shift */
TweakOutcome CEVJ::sensShift(const PropertyTweak<VolPointwise>& shift)
{
    static const string method = "CEVJ::sensShift";
    try {
        if (!Maths::isZero(shift.coefficient)) {
            int expiryIdx = shift.qualifier->expiry->search(ATMVolBM.get());
            double shiftVol = ATMVolArr[expiryIdx] + shift.coefficient;
            if (shiftVol > MAX_DIFF_VOL || shiftVol < MIN_DIFF_VOL) {
                throw ModelException(method, "diffusion vol out of range "
                                     "when tweaking: "
                                     + Format::toString(shiftVol));
            } else {
                ATMVolArr[expiryIdx] = shiftVol;
            }
        }
        return TweakOutcome(shift.coefficient, false);
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             "Volatility tweaking failed for " + getName());
    }
}

/** Shifts the object using given shift */
bool CEVJ::sensShift(JumpRatePointwise* shift)
{
    static const string method = "CEVJ::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            int expiryIdx = shift->getExpiry()->search(ATMVolBM.get());
            double shiftJumpRate = JumpRateArr[expiryIdx] + shiftSize;
            if (shiftJumpRate > MAX_JUMP_RATE ||
                shiftJumpRate < MIN_JUMP_RATE){
                throw ModelException(method, 
                                     "jump rate out of range when tweaking: "
                                     + Format::toString(shiftJumpRate));
            } else {
                JumpRateArr[expiryIdx] = shiftJumpRate;
            }
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             JumpRatePointwise::NAME + 
                             " tweaking failed for " + getName());
    }
    return false;
}

/** Shifts the object using given shift */
bool CEVJ::sensShift(JumpRateParallel* shift)
{
    static const string method = "CEVJ::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            for (int j = 0; j < JumpRateArr.size(); j++)
            {
                double shiftJumpRate = JumpRateArr[j] + shiftSize;
                if (shiftJumpRate > MAX_JUMP_RATE ||
                    shiftJumpRate < MIN_JUMP_RATE) {
                    throw ModelException(method, "jump rate out of range "
                                         "when tweaking: "
                                         + Format::toString(shiftJumpRate));
                } else {
                    JumpRateArr[j] = shiftJumpRate;
                }
            }
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             JumpRateParallel::NAME + 
                             " tweaking failed for " + getName());
    }
    return false; 
}

/** Shifts the object using given shift */
bool CEVJ::sensShift(CEVPowerPointwise* shift)
{
    static const string method = "CEVJ::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            int expiryIdx = shift->getExpiry()->search(ATMVolBM.get());
            double shiftCEVPower = CEVPowerArr[expiryIdx] + shiftSize;
            if (shiftCEVPower > MAX_CEVPOWER || shiftCEVPower < MIN_CEVPOWER){
                throw ModelException(method, 
                                     "CEVPower out of range when tweaking: "
                                     + Format::toString(shiftCEVPower));
            }else {
                CEVPowerArr[expiryIdx] = shiftCEVPower;
            }
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             CEVPowerPointwise::NAME + 
                             " tweaking failed for " + getName());
    }
    return false;
}

/** Shifts the object using given shift */
bool CEVJ::sensShift(CEVPowerParallel* shift)
{
    static const string method = "CEVJ::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            for (int j = 0; j < CEVPowerArr.size(); j++)
            {
                double shiftCEVPower = CEVPowerArr[j] + shiftSize;
                if (shiftCEVPower > MAX_CEVPOWER|| 
                    shiftCEVPower < MIN_CEVPOWER){
                    throw ModelException(method, "CEVPower out of range "
                                         "when tweaking: "
                                         + Format::toString(shiftCEVPower));
                }else{
                    CEVPowerArr[j] = shiftCEVPower;
                }
            }
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             CEVPowerParallel::NAME + 
                             " tweaking failed for " + getName());
    }
    return false; 
}

/** Shifts the object using given shift (see VolParallelShift::Shift)*/
bool CEVJ::sensShift(VolParallelShift* shift)
{
    static const string method = "CEVJ::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            for (int j = 0; j < ATMVolArr.size(); j++)
            {
                double shiftVol = ATMVolArr[j] + shiftSize;
                if (shiftVol > MAX_DIFF_VOL || 
                    shiftVol < MIN_DIFF_VOL){
                    throw ModelException(method, "diffusion vol out of "
                                         "range when tweaking: "
                                         + Format::toString(shiftVol));
                }else{
                    ATMVolArr[j] = shiftVol;
                }
            }
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             "VolParallelShift tweaking failed for "+getName());
    }
    return false;
}

/** Shifts the object using given shift. */
bool CEVJ::sensShift(Theta* shift)
{
    baseDate = shift->rollDate(baseDate);
    return false;  // none of the components have Theta sensitivity
}

void CEVJ::acceptValueDateCollector(const CEVJ*          cevj, 
                                    CValueDateCollector* collector)
{
    collector->valueDateValidate(cevj->baseDate, 
                                 cevj->getName());
}

// Copied from FlatVol::getProcessedVol.
CVolProcessed* CEVJ::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      eqAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const 
{
//    throw ModelException("CEVJ::getProcessedVol", 
//                         "not implemented for currency struck vol");

    static const string routine("CEVJ::getProcessedVol");
    if (!fxAsset){
        throw ModelException(routine, "NULL fx asset");
    }
    if (!eqFXCorr){
        throw ModelException(routine, "NULL correlation");
    }
    /* vol(struck)= sqrt(vol(equity) * vol(equity) + vol(FX) * vol(FX)
       + 2 * correlation * vol(equity) * vol(FX))
    */
    try{
        return CEVJ::getProcessedVol(volRequest,eqAsset);
    } catch (exception& e){
        throw ModelException(e, routine);
    }

}

/** Returns name identifying this object for VolLevel */
string CEVJ::sensName(VolLevel* shift) const {
    return getName();
}
/** Returns name identifying this object for VolParallelShift */
string CEVJ::sensName(VolParallelShift* shift) const {
    return getName();
}

/** Implements VegaParallel  */
string CEVJ::sensName(const VolParallel*) const {
    return getName();
}

/** Implements JumpRateParallel  */
string CEVJ::sensName(JumpRateParallel* shift) const {
    return getName();
}

/** Implements CEVPowerParallel  */
string CEVJ::sensName(CEVPowerParallel* shift) const {
    return getName();
}

/** Implements VegaPointwise */
string CEVJ::sensName(const VolPointwise*) const{
    return getName();
}

/** Implements JumpRatePointwise */
string CEVJ::sensName(JumpRatePointwise* shift) const {
    return getName();
}
/** Implements CEVPowerPointwise */
string CEVJ::sensName(CEVPowerPointwise* shift) const{
    return getName();
}

/** Returns name identifying this object for VolAbsoluteShift */
string CEVJ::sensName(VolAbsoluteShift* shift) const {
    return getName();
}   

/** Shifts the object using given shift (see VolAbsoluteShift::IShift)*/
bool CEVJ::sensShift(VolAbsoluteShift* shift) {
    static const string method("CEVJ::sensShift");
    try {
        for (int i = 0; i < ATMVolBM->size(); i++) {
            double shiftSize = shift->shiftSize(baseDate,(*ATMVolBM)[i]->toDate(baseDate));
            double shiftVol = ATMVolArr[i] + shiftSize;
            if (shiftVol > MAX_DIFF_VOL || 
                shiftVol < MIN_DIFF_VOL){
                throw ModelException(method, "diffusion vol out of "
                                     "range when tweaking: "
                                     + Format::toString(shiftVol));
            }else{
                ATMVolArr[i] = shiftVol;
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method, 
                             "VolAbsoluteShift scenario failed for "+getName());
    }
    return false; // all done
}

/** Returns name identifying this object for VolBenchmarkShift */
string CEVJ::sensName(VolBenchmarkShift* shift) const {
    return getName();
}
/** Shifts the object using given shift (see VolBenchmarkShift::Shift)*/
bool CEVJ::sensShift(VolBenchmarkShift* shift) {
    static const string method = "CEVJ::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)) {
            int  i = 0;
            bool found = false;

            while (i < ATMVolBM->size() && !found) {
                found = shift->expiryEquals((*ATMVolBM)[i].get());
                i++;
            }

            if (!found) {
                throw ModelException(method, 
                                     "benchmark not found on vol surface " + 
                                     getName());
            }

            i--;  // we've stepped over it
                    
            // shift all strikes for this expiry
            double shiftVol = ATMVolArr[i] + shiftSize;
            if( Maths::isNegative(shiftSize) ) {
                // Floor a downshift, but don't increase a low vol
                if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                    ATMVolArr[i] = Maths::min(ATMVolArr[i], 
                                              VolParallelShift::MIN_SPOT_VOL);
                }
                else {
                    ATMVolArr[i] = shiftVol;
                } 
            } else {
                // Don't floor an upshift
                ATMVolArr[i] = shiftVol;
            }

            /* now remove -ve variance
            if (shift->lastShift()) {
                safeFwdVol();
            }*/
        }
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             "VolBenchmarkShift scenario failed for " + 
                             getName());
    }
    return false; // all done
}

/** Returns name identifying this object for PowerVega */
string CEVJ::sensName(PowerVega* shift) const {
    return getName();
}
/** Shifts the object using given shift (see PowerVega::Shift)*/
bool CEVJ::sensShift(PowerVega* shift) {
    static const string method = "CEVJ::sensShift";
    try {
        double shiftSize = shift->getShiftSize();
        if (!Maths::isZero(shiftSize)){
            /* only bother if non zero */
            
            for (int i = 0; i < ATMVolBM->size(); i++) 
            {/* Get trading time to the ith benchmark from cache*/
                double tt = timeMetric->yearFrac(baseDate, (*ATMVolBM)[i]->toDate(baseDate));
                double shiftSize = shift->powerShift(tt);
                double shiftVol = ATMVolArr[i] + shiftSize;
                if( Maths::isNegative(shiftSize) ) {
                    // Floor a downshift, but don't increase a low vol
                    if (shiftVol < VolParallelShift::MIN_SPOT_VOL) {
                        ATMVolArr[i]  = 
                            Maths::min(ATMVolArr[i], 
                                       VolParallelShift::MIN_SPOT_VOL);
                    }
                    else {
                        ATMVolArr[i] = shiftVol;
                    }                 
                } else {
                    // Don't floor an upshift
                    ATMVolArr[i] = shiftVol;
                }
            }
        }

            // now remove -ve variance
//            safeFwdVol();
    }
    catch (exception&e ) {
        throw ModelException(e, method, 
                             "PowerVega scenario failed for " + getName());
    }
    return false; // CEVJ doesn't have a vega type sensitivity ?
}


DRLIB_END_NAMESPACE
