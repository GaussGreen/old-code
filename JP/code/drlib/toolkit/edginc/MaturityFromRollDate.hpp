//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MaturityFromRollDate.hpp
//
//   Description : Defines floating expiries with toDate relative to given set of roll dates
//                 default roll dates are credit IMM dates. Can be used for single name and index swaps
//
//   Author      : Matthias Arnsdorf
//
//   Date        : 11/11/2005
//
//
//----------------------------------------------------------------------------

#ifndef MATURITY_FROM_ROLL_DATE_HPP
#define MATURITY_FROM_ROLL_DATE_HPP

#include "edginc/MaturityPeriod.hpp"


using namespace std;    // string

DRLIB_BEGIN_NAMESPACE

/** Defines floating expiries e.g. 1M, 5Y but relative to a set of IMM dates
    i.e.: toDate(aDate) returns frist Imm date on or after aDate+period.
*/

class TOOLKIT_DLL  MaturityFromRollDate: public Expiry{
public:
    static CClassConstSP const TYPE;
    friend class MaturityFromRollDateHelper;

    MaturityFromRollDate(const string& period);
    MaturityFromRollDate(int count, const string& interval);

    virtual ~MaturityFromRollDate();

    
    /** return this expiry as a string */
    virtual string toString() const;

    /** return this expiry as an absolute date 
	This is the first roll date on or AFTER date+period*/
    virtual DateTime toDate(const DateTime& date) const;

	   
    virtual IObject* clone() const;

    /** Populates interval and count fields */
    void validatePop2Object();

    /** Returns true if given expiry matches this */
    virtual bool equals(const Expiry* expiry) const;
    
   /** get first Roll date on or after date*/
   DateTime getNextRollDate(const DateTime & date) const;

   /** get first Roll date on or after date given rollDate and month */
   static DateTime getNextRollDate(const DateTime & date, int _rollDay, const IntArray & _rollMonths);

private:
    MaturityFromRollDate();
    MaturityFromRollDate(const MaturityFromRollDate& matPeriod);
    // fields
	MaturityPeriodSP maturityPeriod; // $unregistered
    const string period;

	/** Roll day i.e 20 (default)*/
	const int rollDay;
	/** months in which have roll date i.e. 3,6,9,12 for Mar,Jun,Sep,Dec (default) */
	const IntArray rollMonths; // $unregistered
    
};

typedef smartConstPtr<MaturityFromRollDate> MaturityFromRollDateConstSP;
typedef smartPtr<MaturityFromRollDate> MaturityFromRollDateSP;
typedef array<MaturityFromRollDateSP, MaturityFromRollDate> MaturityFromRollDateArray;
typedef smartPtr<MaturityFromRollDateArray> MaturityFromRollDateArraySP;



DRLIB_END_NAMESPACE
#endif
