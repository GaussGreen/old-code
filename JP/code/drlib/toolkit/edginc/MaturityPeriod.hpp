//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MaturityPeriod.hpp
//
//   Description : Defines floating expiries used to define yield curve & 
//                 vol surface points e.g. 1M, 5Y
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef MATURITYPERIOD_HPP
#define MATURITYPERIOD_HPP

#include "edginc/Expiry.hpp"

using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** Defines floating expiries used to define yield curve & 
    vol surface points e.g. 1M, 5Y
    MaturityPeriods are immutable - see explanation in Expiry.hpp.
    DO NOT ADD ANY NON CONST METHODS 
*/

class TOOLKIT_DLL  MaturityPeriod: public Expiry {
public:
    static CClassConstSP const TYPE;
    friend class MaturityPeriodHelper;

    MaturityPeriod(const string& period);
    MaturityPeriod(int count, const string& interval);
	MaturityPeriod(int frequency);

    virtual ~MaturityPeriod();

    /** Returns true if given object matches this */
    bool equalTo(const IObject* expiry) const;

    /** return this expiry as a string */
    virtual string toString() const;

    /** return this expiry as an absolute date */
    virtual DateTime toDate(const DateTime& date) const;

    /** return an expiry as an absolute date */
    static DateTime toDate(int count, const string& interval, const DateTime& date);

    /** as above but takes a const char* - avoid construction of string */
    static DateTime toDate(int count, const char* interval, 
                           const DateTime& date);

    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;
    /** populate an empty object in from reader  */
    virtual void import(Reader::Node* elem, Reader* reader); 

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Override clone method to copy our extra data over */
    virtual IObject* clone() const;

    /** Populates interval and count fields */
    void validatePop2Object();

    /** Returns true if given expiry matches this */
    virtual bool equals(const Expiry* expiry) const;
    
    /** returns a hashcode for the object based upon the period.
        Overridden for performance */
    virtual int hashCode() const;

    /** for use with SwapTool methods */
    void decompose(int& count, string& interval) const;

    /** subtract two dates (d1 - d2) and return difference as a period */
    static MaturityPeriod* dateSubtract(const DateTime& d1,const DateTime& d2);

    /** convert to equivalent number of years */
    double toYears() const;

	/** convert to a number of months. Has to be an exact number of months */
	int toMonths() const;

	/** convert to an annual frequency. Has to be an exact number of months */
	int annualFrequency() const;

	/** convert to an annual frequency */
	int approxAnnualFrequency() const;

    /** test if count is 0 */
    bool isZeroLength() const;

    /** Answers whether this produces dates purely by rolling some number of months */
    static bool isMonthBasedInterval(const string& interval);

private:
    MaturityPeriod();
    MaturityPeriod(const MaturityPeriod& matPeriod);
    // fields
    string period;    // $unregistered
    string interval;  // $unregistered
    int    count;     // $unregistered
};

typedef smartConstPtr<MaturityPeriod> MaturityPeriodConstSP;
typedef smartPtr<MaturityPeriod> MaturityPeriodSP;
typedef array<MaturityPeriodSP, MaturityPeriod> MaturityPeriodArray;
typedef smartPtr<MaturityPeriodArray> MaturityPeriodArraySP;

#ifndef QLIB_MATURITYPERIOD_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<MaturityPeriod>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<MaturityPeriod>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<MaturityPeriodSP _COMMA_ MaturityPeriod>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<MaturityPeriodArray>);
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<MaturityPeriodSP>(
                    MaturityPeriodSP* t));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<MaturityPeriodSP>(MaturityPeriodSP* t,
                                                        IObjectSP o));
#endif

DRLIB_END_NAMESPACE
#endif
