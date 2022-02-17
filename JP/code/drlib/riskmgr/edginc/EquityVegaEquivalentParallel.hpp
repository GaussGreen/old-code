//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : EquityVegaEquivalentParallel.hpp
//
//   Description : ATM Equity Vega derived from Asset vega 
//
//   Author      : André Segger
//
//   Date        : 04 March 2003
//
//
//----------------------------------------------------------------------------

#ifndef EQUITY_VEGA_EQUIVALENT_HPP
#define EQUITY_VEGA_EQUIVALENT_HPP
#include "edginc/Sensitivity.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for ParSpreadRhoEquivalentParalle - a scalar shift where the derivative is
    calculated via a one sided tweak operation */
class RISKMGR_DLL EquityVegaEquivalentParallel: public Sensitivity,
                                    public virtual Additive {
public:
    friend class EquityVegaEquivalentParallelHelper;
    static CClassConstSP const TYPE;
    const static string NAME;

    /** constructor with explicit shift size */
    EquityVegaEquivalentParallel(double shiftSize);

    /** identifies the name used for storing associated results in the output*/
    virtual const string& getSensOutputName() const;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one (return value: false) */
    virtual bool discreteShift() const;

    virtual void calculate(TweakGroup*  tweakGroup,
                           Results*     results);

private:
    double shiftSize;
    /** for reflection */
    EquityVegaEquivalentParallel();
    EquityVegaEquivalentParallel(const EquityVegaEquivalentParallel &rhs);
    EquityVegaEquivalentParallel& operator=(const EquityVegaEquivalentParallel& rhs);
};

typedef smartConstPtr<EquityVegaEquivalentParallel> EquityVegaEquivalentParallelConstSP;
typedef smartPtr<EquityVegaEquivalentParallel> EquityVegaEquivalentParallelSP;

DRLIB_END_NAMESPACE

#endif
