//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : T1Delta.cpp
//
//   Description : Base class for family of T+1 delta shifts
//
//   Author      : Andrew J Swain
//
//   Date        : 21 September 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/T1Delta.hpp"
#include "edginc/ThetaFwdSpot.hpp"

DRLIB_BEGIN_NAMESPACE

/** constructor */
T1Delta:: T1Delta(CClassConstSP clazz, 
                  double        shift, 
                  int           offset, 
                  HolidaySP     hols, 
                  bool          useAssetFwd):
    Sensitivity(clazz), deltaShift(shift), offset(offset), hols(hols),
    useAssetFwd(useAssetFwd) {}

T1Delta::T1Delta(CClassConstSP clazz): Sensitivity(clazz),useAssetFwd(false) {}

/** for reflection */
T1Delta::T1Delta(): Sensitivity(TYPE),useAssetFwd(false) {}

/** populate from market cache */
void T1Delta::getMarket(const IModel* model, const MarketData* market) {
    hols.getData(model, market);
}

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one.  */
bool T1Delta::discreteShift() const{
    return true;
}
  
/** From INextDaySensitivity interface */
ThetaSP T1Delta::offsetRequired(const CInstrument* inst) const{
    if (useAssetFwd) {
        return ThetaSP(new ThetaFwdSpot(offset, HolidaySP(hols.getSP().clone())));
    }
    else {
        return ThetaSP(new Theta(offset, HolidaySP(hols.getSP().clone())));
    }
}
    
/** From INextDaySensitivity interface */
void T1Delta::writeResults(const SensitivitySP&      reqdSens,
                                Results*             dest,
                                const UntweakableSP& untweakable,
                                const Results*       src) const{
    // can just use simple implementation
    bool copyWholePacket = true; 
    Util::writeResults(this, reqdSens, copyWholePacket, dest, untweakable, src);
}

/** calculates given sensitivity - invoked by calculateSens */
void T1Delta::calculate(TweakGroup*  tweakGroup,
                        CResults*    results) {
    CClassConstSP thetaType = useAssetFwd ? ThetaFwdSpot::TYPE : Theta::TYPE;
    // do next day type of greeks together
    Theta::calculateNextDaySens(thetaType, control, tweakGroup, results);
}

class T1DeltaHelper {
public:
    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(T1Delta, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        IMPLEMENTS(Theta::INextDaySensitivity);
        FIELD(deltaShift, "delta shift size");
        FIELD(offset, "number of days to roll");
        FIELD(hols, "holidays");
        FIELD(useAssetFwd, "use asset fwd as next day spot");
        FIELD_MAKE_OPTIONAL(useAssetFwd);
    }
};

CClassConstSP const T1Delta::TYPE = CClass::registerClassLoadMethod(
    "T1Delta", typeid(T1Delta), T1DeltaHelper::load);

DRLIB_END_NAMESPACE
