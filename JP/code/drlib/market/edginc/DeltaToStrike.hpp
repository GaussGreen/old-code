//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaToStrike.hpp
//
//   Description : Delta to Strike Interface
//
//   Author      : Manos Venardos
//
//   Date        : 5 October 2001
//
//
//----------------------------------------------------------------------------

#ifndef DELTA_TO_STRIKE_HPP
#define DELTA_TO_STRIKE_HPP

#include "edginc/Object.hpp"
#include "edginc/Class.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for class that can create an IDeltaToStrike */
class MARKET_DLL IDeltaToStrikeMaker: virtual public IObject {
public:
    /** Interface for converting a delta denominated strike to an absolute */
    class MARKET_DLL IDeltaToStrike {
    public:
        /** Computation method */
        virtual double calcStrike(double lowerStrike,
                                  double upperStrike,
                                  double strikeAbsAcc) const = 0;

        /** Virtual destructoe */
        virtual ~IDeltaToStrike();
    };

    DECLARE_REF_COUNT(IDeltaToStrike);

    static CClassConstSP const TYPE;

    /** Factory method */
    virtual const IDeltaToStrike* make(const DateTime&             valueDate,
                                       const DateTime&             maturityDate,
                                       bool                        isCall,
                                       const CAsset*               asset,
                                       const YieldCurve*           discount,
                                       const InstrumentSettlement* settle,
                                       double                      deltaShiftSize,
                                       double                      tgtDelta,
                                       const string&               volType,
                                       bool                        allowNegativeFwdVar) = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartConstPtr<IDeltaToStrikeMaker> IDeltaToStrikeMakerConstSP;
typedef smartPtr<IDeltaToStrikeMaker> IDeltaToStrikeMakerSP;

DRLIB_END_NAMESPACE

#endif
