//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BorrowCurve.hpp
//
//   Description : borrowing curve
//
//   Author      : Stephen Hope
//
//   Date        : 21 Feb 2001
//
//
//---------------------------------------------------------------------------


#ifndef EDGBORROWCURVE_HPP
#define EDGBORROWCURVE_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RhoBorrowParallel.hpp"
#include "edginc/RhoBorrowPointwise.hpp"
#include "edginc/PrivateObject.hpp"
#include "edginc/PublicObject.hpp"
#include "edginc/Theta.hpp"
#include "edginc/BorrowParallelShift.hpp"
#include "edginc/BorrowLevel.hpp"
using namespace std;

DRLIB_BEGIN_NAMESPACE


/** The borrowing curve contains a time series of rates.
    Default methods are provided for interpolation
    on the curve but note that 'time of day' is not considered in the 
    calculation */
class MARKET_DLL BorrowCurve : public CObject, 
                    public virtual IPrivateObject,
                    public virtual Theta::IShift
{
public:
    static CClassConstSP const TYPE;
    friend class BorrowCurveHelper;

    /** Constructor */
    BorrowCurve(const DateTime&     baseDate,
                const ExpiryArray*  expiries,
                const DoubleArray&  rates,
                const string&       name,
                const int           basis = CompoundBasis::CONTINUOUS,
                DayCountConvention* dayCountConv = NULL);
    
    virtual ~BorrowCurve();

    /**
     * Get data from market
     *
     * Note: BorrowCurve is not a MarketObject (maybe it should be).
     */

    virtual void getMarket(const IModel* extModel, const MarketData* market);

    /** returns a linear interpolation at the required date */
    double interpAtDate(const DateTime& interpDate)const;

    double continuousFwdRate(const DateTime& startDate,
                             const DateTime& endDate) const;
    
    /** returns the forward borrowing rate between two dates. */
    double fwdRate(const DateTime&  startDate,
                   const DateTime&  endDate)const;

    /** returns the forward borrowing rate between two dates. Day
        count convention and compounding basis can be chosen */
    double fwdRate(const DateTime&           startDate,
                   const DateTime&           endDate,
                   const DayCountConvention* dcc,
                   const int                 fwdBasis)const;

    /** Returns the basis of the borrowing curve */
    int getBasis() const;

    const DoubleArray&   getRates() const { return rates; }
    const DateTimeArray& getDates() const { return dates; }

    /** Returns the name of the borrow curve - used to determine
        whether to tweak the object */
    string sensName(RhoBorrowParallel* shift)const;

    /** Shifts the object using given shift */    
    bool sensShift(RhoBorrowParallel* shift);

    /** Restores the object to its original form */
    void sensRestore(RhoBorrowParallel* shift);

    /** Returns the name of the borrow curve - used to determine whether 
        to tweak the object */
    string sensName(RhoBorrowPointwise* shift)const;
    
    /** Return the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this  borrow curve */
    ExpiryArrayConstSP sensExpiries(RhoBorrowPointwise* shift)const;
    
    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    bool sensShift(RhoBorrowPointwise* shift);

    /** Restores the object to its original form */
    void sensRestore(RhoBorrowPointwise* shift);  
    
    /** Shifts the object using given shift. */
    bool sensShift(Theta* shift);

    /** Returns the name of the borrow curve - used to determine
        whether to tweak the object */
    string sensName(BorrowParallelShift* shift)const;

    /** Shifts the object using given shift */    
    bool sensShift(BorrowParallelShift* shift);

    /** Returns the name of the borrow curve - used to determine
        whether to tweak the object */
    string sensName(BorrowLevel* shift)const;

    /** Shifts the object using given shift */    
    bool sensShift(BorrowLevel* shift);

    virtual IPublicObject* toPublicObject() const;

    /** are all the rates zero ? */
    bool isZero() const;        

    /** Returns value date */
    const DateTime& getBaseDate() const;

    /** define the public interface to this class */
    class MARKET_DLL Interface : public CObject, public IPublicObject {
    public:
        static CClassConstSP const TYPE;
        friend class BorrowCurveInterfaceHelper;
        virtual IPrivateObject* toPrivateObject() const;

        Interface(
            const DateTime&    baseDate,
            const ExpiryArray* expiries,
            const DoubleArray& rates,
            const string&      name,
            const int          basis,
            const string&      dayCount);

    private:
        Interface();

        DateTime       baseDate;
        ExpiryArraySP  expiries;
        DoubleArray    rates;
        string         name;
        int            basis;
        string         dayCount;  
    };

private:
    BorrowCurve();
    BorrowCurve(const BorrowCurve& rhs);
    BorrowCurve& operator=(const BorrowCurve& rhs);
    virtual void validatePop2Object();
    void buildDateCache();

    DateTime       baseDate;
    ExpiryArraySP  expiries;
    DoubleArray    rates;
    string         name;
    int            basis;
    DayCountConventionConstSP dayCountConv;

    DateTimeArray dates; // cached - internally derived from expiries

    // internal use for interpolation
    mutable int loBound; // $unregistered
    mutable int hiBound; // $unregistered
};

typedef smartConstPtr<BorrowCurve> BorrowCurveConstSP;
typedef smartPtr<BorrowCurve> BorrowCurveSP;

typedef smartConstPtr<BorrowCurve> BorrowCurveConstSP;
typedef smartPtr<BorrowCurve> BorrowCurveSP;

DRLIB_END_NAMESPACE
#endif




