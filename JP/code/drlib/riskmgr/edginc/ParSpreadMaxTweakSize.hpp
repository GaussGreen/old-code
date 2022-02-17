//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : ParSpreadMaxTweakSize.hpp
//
//   Description : Pseudo sensitivity that uses the qualifier mechanism
//                 to return maximum tweak sizes that may be applied to
//                 the par spread curve.
//
//   Author      : Gordon Stephens
//
//   Date        : 22 June 2005
//

//
//----------------------------------------------------------------------------

#ifndef PAR_SPREAD_MAX_TWEAK_SIZE_HPP
#define PAR_SPREAD_MAX_TWEAK_SIZE_HPP

#include "edginc/TweakQualifierID.hpp"
#include "edginc/TweakNameResolver.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/IRiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL ParSpreadMaxTweakSize : public virtual ITweakQualifierID,
                              public virtual ITweakNameResolver{
public:
    //What an implementing class must provide
    class RISKMGR_DLL IShift {
    public:
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(ParSpreadMaxTweakSize* shift) const = 0;

        /** Return the array of expiries (ie maturities/benchmark dates) that
            need to be tweaked for this  yield curve */
        virtual DoubleArrayConstSP sensMaxTweakSize(
            ParSpreadMaxTweakSize* shift) const = 0;
    };

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    virtual CClassConstSP shiftInterface() const;

    /** Returns this */
    virtual ITweakNameResolver* nameResolver();

    /** returns the name identifying the market data to be shifted. */
    virtual OutputNameConstSP getMarketDataName() const;

    /** does the given name match the name identifying the market data
        to be shifted. Could consider being const, but has large downstream
        consequences. */
    virtual bool nameMatches(const OutputName& name,
                             IObjectConstSP    obj);

    /** Casts supplied obj to IShift and calls sensMaxTweakSize() */
    virtual IObjectConstSP qualifier(IObjectConstSP obj);

    /** Calculates appropriate shift sizes for each expiry given original
     * shift size */
    DoubleArrayConstSP calculateTweakSizes(
        const IObject*     tweakGroup,
        DoubleArrayConstSP origShiftSizes,
        OutputNameConstSP  name,
        const ExpiryArray& expiries);

    //-----------------------------
    //ParSpreadMaxTweakSize methods
    //-----------------------------
    IExpiryRiskPropertyConstSP wrt();

    ParSpreadMaxTweakSize(OutputNameConstSP name,
                          IExpiryRiskPropertyConstSP withRespectTo);

    virtual ~ParSpreadMaxTweakSize();

    static IExpiryRiskPropertyConstSP adapted(
               IExpiryRiskPropertyConstSP property);

private:
    //fields
    OutputNameConstSP name;
    IExpiryRiskPropertyConstSP withRespectTo; //holds information about which concrete tweak
    //wants this information, which allows some specialisation
};

DRLIB_END_NAMESPACE

#endif
