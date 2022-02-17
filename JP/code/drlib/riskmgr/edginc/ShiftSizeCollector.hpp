//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ShiftSizeCollector.hpp
//
//   Description : Delta Shift size collector class
//
//   Author      : André Segger
//
//   Date        : 20 Jun 2001
//
//
//----------------------------------------------------------------------------

#ifndef SHIFTSIZE_COLLECT_HPP
#define SHIFTSIZE_COLLECT_HPP
#include <string>
#include <set>
#include "edginc/smartPtr.hpp"
#include "edginc/Collector.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Delta.hpp"

DRLIB_BEGIN_NAMESPACE

/** A class to help determine a good delta shift size */
class RISKMGR_DLL ShiftSizeCollector: public CObject, public virtual ICollector {
public:
    friend class ShiftSizeCollHelper;
    static CClassConstSP const TYPE;

    enum TAdjustmentType
    {
        FWD_START_ADJUSTMENT = 0, /* means instrument will be interpolating
                                     at the spot price */
        SPOT_START_ADJUSTMENT     /* means instrument will be interpolating
                                     at the strike */
    };

    /** Vols should call this with their asset's current spot price 
        through the IVolDeltaShiftSize interface. 
        If the information is not available
        then this routine should not be called. This method then
        alters the shift size in view of the information. */
    void adjustShiftSize(const DoubleArraySP& strikeList,
                         const string&        assetName,
                         const double         spot);

    /** Creates a shift size collector that can be passed to an asset's
        accept method */
    ShiftSizeCollector(const Delta*          sensControl,
                       const double          strike,
                       const TAdjustmentType adjType);

    /** returns the altered shift size */
    double          getShiftSize() const;

    /** sets the shift size value */
    void            setShiftSize(double newShiftSize);

    /** returns the adjustment type */
    TAdjustmentType getAdjustmentType() const;

    /** returns the strike*/
    double getStrike() const;

    /** returns the original sens control */
    DeltaConstSP  getSensControl() const;

private:
    DeltaConstSP        sensCtrl; // $unregistered
    double              shiftSize; // $unregistered
    double              strike; // $unregistered
    TAdjustmentType     adjType; // $unregistered

    static void load(CClassSP& clazz);
    ShiftSizeCollector(const ShiftSizeCollector& rhs);
    ShiftSizeCollector& operator=(const ShiftSizeCollector& rhs);
};

typedef smartPtr<ShiftSizeCollector>       ShiftSizeCollectorSP;
typedef smartPtr<const ShiftSizeCollector> ShiftSizeCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
