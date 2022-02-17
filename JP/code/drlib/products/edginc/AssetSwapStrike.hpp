//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : AssetSwapStrike.hpp
//
//   Description : base class for asset swap strike objects 
//
//   Date        :02 June 2003 
//
//
//----------------------------------------------------------------------------

#ifndef EDR_ASSET_SWAP_STRIKE_HPP
#define EDR_ASSET_SWAP_STRIKE_HPP

#include "edginc/Object.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ConvBond.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL AssetSwapStrike: public CObject {
public:
    static CClassConstSP const TYPE;

    /** calculate the asset swap strike on a given date */
    virtual double getStrike(const DateTime&    baseDate,
                             ConvBondConstSP    cvb,
                             YieldCurveConstSP  discountCurve) const = 0;

    /** returns true if the asset swap is exercisable on the baseDate */
    virtual bool   isExercisable(const DateTime& baseDate) const = 0;

    virtual ~AssetSwapStrike() {};

protected:
    AssetSwapStrike(CClassConstSP clazz): CObject(clazz) {}
private:
    AssetSwapStrike();

};

typedef smartPtr<AssetSwapStrike>      AssetSwapStrikeSP;
typedef smartConstPtr<AssetSwapStrike> AssetSwapStrikeConstSP;


DRLIB_END_NAMESPACE

#endif // EDR_ASSET_SWAP_STRIKE_HPP

