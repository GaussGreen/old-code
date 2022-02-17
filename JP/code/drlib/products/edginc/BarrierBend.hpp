//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BarrierBend.hpp
//
//   Description : 
//
//   Date        : Q Hou Nov 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_BarrierBend_HPP
#define EDR_BarrierBend_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Array.hpp"

DRLIB_BEGIN_NAMESPACE


/*************************************************************/
// Allow some flexibility in styles of barrier bending information

class PRODUCTS_DLL IBarrierBend {
public:
    virtual ~IBarrierBend() {};

    // get bending for simDate from the schedule bent for bendDate
    virtual double getBending(int endDateIdx,
                              int simDateIdx) const = 0;

    // get start date index for the schedule bent for bendDate, and
    // get pointer to array of bending for schedule bent for bendDate
    // array size is (bendDateIdx - startDateIdx + 1)
    virtual bool getBending(int endDateIdx, 
                           DoubleArray &bendAmt) const = 0;

    // get array of end date index for bent schedule originating at the startDate
    virtual IntArray getBendingEndByStart(int startDateIdx) const = 0;

    // get array of end date index for bent schedule covering the startDate
    virtual IntArray getBendingEndByCover(int simDateIdx) const = 0;

    // get bend period start given period end
    virtual int getBendingStartByEnd(int endDateIdx) const = 0;

    // query for start/end index of the bent periods. sorted and no duplicate
    virtual IntArray getBendingStart() const = 0;

    virtual IntArray getBendingEnd() const = 0;
};

typedef refCountPtr<IBarrierBend> IBarrierBendSP;


// Build IBarrierBend via a Maker class.
class PRODUCTS_DLL IBarrierBendMaker : virtual public IObject {
public:
    static CClassConstSP const TYPE;
    
    // value date is used to allow removing historical bending
    virtual IBarrierBendSP getBarrierBend(const DateTime& valueDate, const DateTimeArray &simDates) const = 0;

    virtual bool trivial() const=0;

    virtual ~IBarrierBendMaker() {};

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IBarrierBendMaker> IBarrierBendMakerSP;

DRLIB_END_NAMESPACE

#endif

