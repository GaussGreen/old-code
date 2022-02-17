//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : T1Delta.hpp
//
//   Description : Base class for family of T+1 delta shifts
//
//   Author      : Andrew J Swain
//
//   Date        : 21 September 2004
//
//
//----------------------------------------------------------------------------

#ifndef T1DELTA_HPP
#define T1DELTA_HPP

#include "edginc/Holiday.hpp"
#include "edginc/SensControl.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

/** T+1 delta shift */
class RISKMGR_DLL T1Delta: public Sensitivity,
               public virtual Additive,
               public virtual Theta::INextDaySensitivity{
public:
    static CClassConstSP const TYPE;

    /** constructor */
    T1Delta(CClassConstSP clazz);
    T1Delta(CClassConstSP clazz, 
            double        shift, 
            int           offset, 
            HolidaySP     hols, 
            bool          useAssetFwd);

    /** identifies the name used storing associated results in the output 
        (echo pure virtual function in Sensitivity class) */
    virtual const string& getSensOutputName() const = 0;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns true */
    virtual bool discreteShift() const;

    /** Echo from INextDaySensitivity interface */
    virtual SensitivitySP requiredSensitivity(TweakGroup* tweakGroup) const = 0;
    
    /** From INextDaySensitivity interface */
    virtual ThetaSP offsetRequired(const CInstrument* inst) const;
    
    /** From INextDaySensitivity interface */
    virtual void writeResults(const SensitivitySP& reqdSens,
                              Results*             dest,
                              const UntweakableSP& untweakable,
                              const Results*       src) const;
protected:
    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*  tweakGroup,
                           CResults*    results);

    double         deltaShift;
    mutable int    offset; // mutable fpr DeltaTplusN
    HolidayWrapper hols;
    bool           useAssetFwd;

private:
    /** for reflection */
    T1Delta();
    friend class T1DeltaHelper;
};

DRLIB_END_NAMESPACE

#endif
