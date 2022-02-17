//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BucketShift.hpp
//
//   Description : Bucket shift interface 
//
//   Author      : Stephen Hope
//
//   Date        : 14 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_BUCKET_SHIFT_H
#define EDG_BUCKET_SHIFT_H
#include "edginc/ScalarShift.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

/**  */
class RISKMGR_DLL BucketShift: public ScalarShift{
public:    
    friend class BucketShiftHelper;
    static CClassConstSP const TYPE;
    
    /** Returns the expiry which is currently being tweaked */
    ExpiryConstSP getExpiry()const;

    /** Returns the expiries which are to be tweaked. */
    ExpiryArrayConstSP getExpiries() const;

    /** returns the start- and end-dates for the current bucket */
    void getBucketDates(DateTime& startDate, DateTime& endDate) const;

    /** implements a one sided bucket derivative for each instance of the
        market data which is sensitive to this SensControl */
    void calculate(TweakGroup*      tweakGroup,
                   CResults*        results);

protected:
    /** Returns the expiries which are to be tweaked */
    ExpiryArrayConstSP calculateExpiries();

   /** returns a sorted expiryArray which excludes any past date */
    void getRevisedExpiries(const DateTime& valueDate,
                            ExpiryArraySP& revisedList,
                            ExpiryArraySP& pastDateList);

    /** Create a bucket shift of type clazz,
        which uses outputName (eg MU_S) to identify results and
        with given shiftSize and expiry dates */
    BucketShift(const CClassConstSP& clazz,
                const string&        sensName,
                const double&        shiftSize,
                const ExpiryArray*   expiries);


    /** for reflection */
    BucketShift(const CClassConstSP& clazz,
                const string&        sensName);

private:
    ExpiryConstSP  expiry;
    ExpiryArraySP  expiries;
    DateTime       bucketStartDate;
    DateTime       bucketEndDate;
    BucketShift(const BucketShift& rhs);
    BucketShift& operator=(const BucketShift& rhs);
};

DRLIB_END_NAMESPACE

#endif
