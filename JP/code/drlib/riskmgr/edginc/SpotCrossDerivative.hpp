//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpotCrossDerivative.hpp
//
//   Description : Interface for sensitivities that involve tweaking two
//   (or more?) pieces of market data (eg FXCrossGamma, CrossGamma) in order 
//   to tweak a cross derivative. Used by 'Quick X Gamma' in Monte Carlo
//
//   Author      : Mark A Robson
//
//   Date        : 8 November 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_CROSS_DERIV_HPP
#define EDR_CROSS_DERIV_HPP

#include "edginc/ScalarShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for sensitivities that involve tweaking two (or more?)
    pieces of market data such that an asset spot price is shifted (eg
    FXCrossGamma, CrossGamma) in order to tweak a cross
    derivative. Used by 'Quick X Gamma' in Monte Carlo */
class RISKMGR_DLL ISpotCrossDerivative: public virtual IObject  {
public:
    static CClassConstSP const TYPE;
    /** Returns the shifts which define the 'area of operation of the
     shifts'. This means the names of what is going to be shifted, the
     largest shift size for each shift and the nature of the shift
     (which is captured by the type of the SensControl). The names
     should be stored so that they can be rerieved using the names()
     method. The shifts returned must be distinct (ie not SP to the
     same object). Additionally the names within the ScalarShifts must
     not have blanks or duplicates*/
    virtual ScalarShiftArray definingShifts(CInstrument*   instrument,
                                            IModel*        model,
                                            Control*       control) const = 0;

    virtual ~ISpotCrossDerivative(); // empty

};
typedef smartPtr<ISpotCrossDerivative> ISpotCrossDerivativeSP;

DRLIB_END_NAMESPACE
#endif
