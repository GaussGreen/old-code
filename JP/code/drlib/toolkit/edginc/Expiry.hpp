//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Expiry.hpp
//
//   Description : Defines interface to expiries used to define yield curve & 
//                 vol surface points
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef EXPIRY_HPP
#define EXPIRY_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** Defines interface to expiries used to define yield curve & 
    vol surface points etc. Note that expiries are immutable, ie, once created
    they cannot be changed - this fact is used to boost performance when
    copying expiries

    DO NOT ADD ANY NON CONST METHODS 

    (except for validatePop2Object which is
    essentially part of the construction of the object)
*/
class Expiry;
typedef smartConstPtr<Expiry> ExpiryConstSP;
typedef smartPtr<Expiry> ExpirySP;

typedef array<ExpirySP, Expiry> ExpiryArray;
typedef smartPtr<ExpiryArray> ExpiryArraySP;
typedef smartConstPtr<ExpiryArray> ExpiryArrayConstSP;

typedef array<ExpiryArraySP, ExpiryArray> ExpiryArrayArray;
typedef smartPtr<ExpiryArrayArray> ExpiryArrayArraySP;
typedef smartConstPtr<ExpiryArrayArray> ExpiryArrayArrayConstSP;

class TOOLKIT_DLL Expiry : public CObject {
public:
    static CClassConstSP const TYPE;
    friend class ExpiryHelper;

    virtual ~Expiry();

    /** return this expiry as a string */
    virtual string toString() const = 0;

    /** return this expiry as an absolute date */
    virtual DateTime toDate(const DateTime& date) const = 0;

    /** Returns true if given expiry matches this exactly (including type
        of expiry match) */
    virtual bool equals(const Expiry* expiry) const = 0;

    /** Returns true if toDate(base) would match for both expiries */
    bool equals(const DateTime& base, const Expiry*   expiry) const;

    /** Returns index of this expiry in supplied array which is of type
        const ExpiryArray*    */
    int search(const array<smartPtr<Expiry>, Expiry>* expiries) const;

    /** Returns index of this expiry in supplied array which is of type
        const ExpiryArray*. The supplied base date is used in the
        comparison so that it is possible that eg MaturityPeriods
        are contained in arrays of BenchmarkDates */
    int search(const DateTime& base,
               const array<smartPtr<Expiry>, Expiry>* expiries) const;

    /** Report if the part ExpiryArray is a subset of the full ExpiryArray,
        using the base DateTime to convert into dates if required, ie, so
        that it is possible to compare MaturityPeriods and BenchmarkDates */
    static bool isSubset(const DateTime& base,
                         const array<smartPtr<Expiry>, Expiry>& full, 
                         const array<smartPtr<Expiry>, Expiry>& part);

    /** Report if the part ExpiryArray is a subset of the full DateTimeArray,
        using the base DateTime to convert expiries to dates */
    static bool isSubset(const DateTime& base,
                         const DateTimeArray& full, 
                         const array<smartPtr<Expiry>, Expiry>& part);

    /** Compares two arrays of expiries. The supplied base date is used in the
        comparison so that it is possible that an array of eg MaturityPeriods
        could equal an array of BenchmarkDates */
    static bool equals(const DateTime& base,
                       const array<smartPtr<Expiry>, Expiry>* expiries1,
                       const array<smartPtr<Expiry>, Expiry>* expiries2);

    /** Compares two arrays of expiries (including type of expiry match)*/
    static bool equals(const array<smartPtr<Expiry>, Expiry>* expiries1,
                       const array<smartPtr<Expiry>, Expiry>* expiries2);

    /** Merges two arrays of expiries. Where the toDate() method gives the
        same value (in conjunction with the supplied base) only the expiry
        from expiries1 is retained. Both arrays must be in ascending order with
        no duplicates. The results are not defined if this is not the case.
        Also both arrays must not contain null expiries */
    static const array<smartPtr<Expiry>, Expiry> merge(
        const DateTime& base,
        const array<smartPtr<Expiry>, Expiry>* expiries1,
        const array<smartPtr<Expiry>, Expiry>* expiries2);

    /** is an array of expiries strictly increasing ? */
    static bool isIncreasing(const DateTime& base,
                             const array<smartPtr<Expiry>, Expiry>* expiries);

    //// Creates an Expiry object from a string.
    static ExpiryConstSP createExpiry(const string& input);

protected:
    Expiry(CClassConstSP clazz);
};


#ifndef QLIB_EXPIRY_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<Expiry>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<Expiry>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<ExpirySP _COMMA_ Expiry>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<ExpiryArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<ExpiryArray>);
// again, avoid code bloat by declaring common implementation of function
// templates as extern
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<Expiry>(Expiry* t, IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<ExpiryArray>(ExpiryArray* t, IObjectSP o));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<ExpirySP>(ExpirySP* t));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<ExpiryArraySP>(ExpiryArraySP* t));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<ExpirySP>(ExpirySP* t, IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<ExpiryArraySP>(ExpiryArraySP* t,
                                                     IObjectSP o));
#endif

DRLIB_END_NAMESPACE
#endif
