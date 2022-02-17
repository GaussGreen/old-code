//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BenchmarkDate.hpp
//
//   Description : Defines fixed date expiries used to define yield curve & 
//                 vol surface points
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef BENCHMARKDATE_HPP
#define BENCHMARKDATE_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Expiry.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** Defines fixed date expiries used to define yield curve & 
    vol surface points. BenchmarkDates are immutable - see explanation
    in Expiry.hpp.
    DO NOT ADD ANY NON CONST METHODS 
*/
class TOOLKIT_DLL BenchmarkDate : public Expiry {
public:
    static CClassConstSP const TYPE;
    friend class BenchmarkDateHelper;

    BenchmarkDate(const DateTime& date);

    /** Override clone method for performance */
    virtual IObject* clone() const;

    virtual ~BenchmarkDate();

    /** return this expiry as a string */
    virtual string toString() const;

    /** return this expiry as an absolute date */
    virtual DateTime toDate(const DateTime& date) const;

    /** specific to BenchmarkDate - return this expiry as an absolute date */
    DateTime toDate() const;

    /** Returns true if given object matches this */
    virtual bool equalTo(const IObject* expiry) const;

    /** Returns true if given expiry matches this */
    virtual bool equals(const Expiry* expiry) const;

    /** returns a hashcode for the object based upon the date and time.
        Overridden for performance */
    virtual int hashCode() const;

    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;
    /** populate an empty object in from reader  */
    virtual void import(Reader::Node* elem, Reader* reader); 

private:
    BenchmarkDate();
    BenchmarkDate(const BenchmarkDate& bm);
    // fields
    const DateTime date;
};

typedef smartConstPtr<BenchmarkDate> BenchmarkDateConstSP;
typedef smartPtr<BenchmarkDate> BenchmarkDateSP;
#ifndef QLIB_BENCHMARKDATE_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<BenchmarkDate>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<BenchmarkDate>);
#endif

DRLIB_END_NAMESPACE
#endif
