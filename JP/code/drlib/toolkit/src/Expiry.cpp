//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Expiry.cpp
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

#include "edginc/config.hpp"
#define QLIB_EXPIRY_CPP
#include "edginc/Expiry.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Null.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/BenchmarkDate.hpp"
#include <boost/tokenizer.hpp>

DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<Expiry>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<Expiry>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<ExpirySP _COMMA_ Expiry>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<ExpiryArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<ExpiryArray>);
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<Expiry>(Expiry* t, IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<ExpiryArray>(ExpiryArray* t,
                                                      IObjectSP o));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<ExpirySP>(ExpirySP* t));
INSTANTIATE_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<ExpiryArraySP>(
                         ExpiryArraySP* t));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<ExpirySP>(ExpirySP* t, IObjectSP o));
INSTANTIATE_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<ExpiryArraySP>(ExpiryArraySP* t,
                                                          IObjectSP o));

Expiry::~Expiry() {
    // empty
}

Expiry::Expiry(CClassConstSP clazz): CObject(clazz){}

/** Returns index of this expiry in supplied array which is of type
    const ExpiryArray*    */
int Expiry::search(const array<smartPtr<Expiry>, Expiry>* expiries) const{
    static const string routine = "Expiry::search";
    if (!expiries){
        throw ModelException(routine, "null array");
    }
    for (int idx = 0; idx < expiries->size(); idx++){
        if (equals((*expiries)[idx].get())){
            return idx;
        }
    }
    throw ModelException(routine,"Expiry ("+toString()+") not found");
}


/** Returns index of this expiry in supplied ExpiryArray.
    The supplied base date is used in the comparison so that it is 
    possible to search for eg MaturityPeriods in arrays of 
    BenchmarkDates */
int Expiry::search(const DateTime& base,
                   const array<smartPtr<Expiry>, Expiry>* expiries) const
{
    static const string routine = "Expiry::search(base...)";
    if (!expiries){
        throw ModelException(routine, "null array");
    }
    const CClassConstSP myClass = getClass();
    const DateTime& myDate = toDate(base);
    ExpirySP expiryIdx;
    int expiriesSize = expiries->size();

    for (int idx = 0; idx<expiriesSize; idx++){
        expiryIdx = (*expiries)[idx];
        if (myClass == expiryIdx->getClass()) {
            if (equals(expiryIdx.get())) {
                return idx;
            }
        }
        else {
            if (myDate.equals(expiryIdx->toDate(base))) {
                return idx;
            }
        }
    }
    throw ModelException(routine, "Expiry ("+toString()+") not found");
}


/** Report if the part ExpiryArray is a subset of the full ExpiryArray,
    using the base DateTime to convert into dates if required, ie, so
    that it is possible to compare MaturityPeriods and BenchmarkDates */
bool Expiry::isSubset(const DateTime& base,
                      const ExpiryArray& full, 
                      const ExpiryArray& part) 
{
    int pSize = part.size();
    int fSize = full.size();
    int f = 0;
    for (int p=0; p<pSize; p++) {
        // Search all elements of part in full sequentially (ie, they must be
        // in the same order - typically, in increasing order)
        for (/*f*/; f<fSize; f++) {
            if (full[f]->getClass() == part[p]->getClass()) {
                if (part[p]->equals(full[f].get())) {
                    break; // Get out of the for loop in f
                }
                // else skip expiries in full which are not present in part
            }
            else {
                const DateTime& fDate = full[f]->toDate(base);
                const DateTime& pDate = part[p]->toDate(base);
                if (fDate.equals(pDate)) {
                    break; // Get out of the for loop in f
                }
                // else skip expiries in full which are not present in part
            }
        }
        if (f == fSize) {
            // It means we did not find part[p] in full, ie, part is not a subset
            return false;
        }
    }
    return true;
}


/** Report if the part ExpiryArray is a subset of the full DateTimeArray,
    using the base DateTime to convert expiries to dates */
bool Expiry::isSubset(const DateTime& base,
                      const DateTimeArray& full, 
                      const ExpiryArray& part) 
{
    int pSize = part.size();
    int fSize = full.size();
    int f = 0;
    for (int p=0; p<pSize; p++) {
        DateTime pDate = part[p]->toDate(base);
        // Search all elements of part in full, sequentially (ie, they must be
        // in the same order - typically, in increasing order)
        for (/*f*/; f<fSize; f++) {
            if (pDate.equals(full[f])) {
                break; // Get out of the for loop in f
            }
            // else skip expiries in full which are not present in part
        }
        if (f == fSize) {
            // It means we did not find part[p] in full, ie, part is not a subset
            return false;
        }
    }
    return true;
}


/** Compares two arrays of expiries (including type
    of expiry match)*/
bool Expiry::equals(const array<smartPtr<Expiry>, Expiry>* expiries1,
                    const array<smartPtr<Expiry>, Expiry>* expiries2){
    static const string routine = "Expiry::equals";
    if (!expiries1 || !expiries2){
        throw ModelException(routine, "null array");
    }
    if (expiries1->size() != expiries2->size()){
        return false;
    }
    const ExpiryArray& expy1 = *expiries1;
    const ExpiryArray& expy2 = *expiries2;
    for (int i = 0; i < expy1.size(); i++){
        if (!expy1[i]->equals(expy2[i].get())){
            return false;
        }
    }
    return true;
}

/** Compares two arrays of expiries. The supplied base date is used in the
    comparison so that it is possible that an array of eg MaturityPeriods
    could equal an array of BenchmarkDates */
bool Expiry::equals(const DateTime& base,
                    const array<smartPtr<Expiry>, Expiry>* expiries1,
                    const array<smartPtr<Expiry>, Expiry>* expiries2){
    static const string routine = "Expiry::equals(base...)";
    if (!expiries1 || !expiries2){
        throw ModelException(routine, "null array");
    }
    if (expiries1->size() != expiries2->size()){
        return false;
    }
    const ExpiryArray& expy1 = *expiries1;
    const ExpiryArray& expy2 = *expiries2;
    for (int i = 0; i < expy1.size(); i++){
        if (expy1[i]->getClass() == expy2[i]->getClass()){
            if (!expy1[i]->equals(expy2[i].get())){
                return false;
            }
        } else {
            if (!expy1[i]->toDate(base).equals(expy2[i]->toDate(base))){
                return false;
            }
        }
    }
    return true;
}

/** Returns true if toDate(base) would match for both expiries */
bool Expiry::equals(const DateTime& base,
                    const Expiry*   expiry) const{
    if (!expiry){
        return false;
    }
    // optimise performance
    if (getClass() == expiry->getClass()){
        return equals(expiry);
    }
    return (toDate(base).equals(expiry->toDate(base)));
}


/** Merges two arrays of expiries. Where the toDate() method gives the
    same value (in conjunction with the supplied base) only the expiry
    from expiries1 is retained. Both arrays must be in ascending order with
    no duplicates. The results are not defined if this is not the case.
    Also both arrays must not contain null expiries */
const array<smartPtr<Expiry>, Expiry> Expiry::merge(
    const DateTime& base,
    const array<smartPtr<Expiry>, Expiry>* expiries1,
    const array<smartPtr<Expiry>, Expiry>* expiries2)
{
    static const string routine = "Expiry::merge";
    if (!expiries1 || !expiries2){
        throw ModelException(routine, "null array");
    }
    const ExpiryArray& expy1 = *expiries1;
    const ExpiryArray& expy2 = *expiries2;
    // create list at least big enough
    ExpiryArray list(0);
    list.reserve(expy1.size() + expy2.size());
    int pos1, pos2;
    for (pos1 = 0, pos2 = 0; pos1 < expy1.size() && pos2 < expy2.size();){
        /* performance (assuming that most of the time
           the expiries are similar) */
        if (expy1[pos1]->equals(base, expy2[pos2].get())){
            list.push_back(expy1[pos1++]);
            pos2++;
        } else if (expy1[pos1]->toDate(base).
                   isGreater(expy2[pos2]->toDate(base))){
            list.push_back(expy2[pos2++]);
        } else {
            list.push_back(expy1[pos1++]);
        }
    }
    return list;
}

/** is an array of expiries strictly increasing ? */
bool Expiry::isIncreasing(const DateTime& base,
                          const array<smartPtr<Expiry>, Expiry>* expiries) {
    static const string method("Expiry::isIncreasing");
    try {
        if (!expiries) {
            throw ModelException(method, "expiries is null");
        }
        if (expiries->empty()) {
            throw ModelException(method, "expiries is empty");
        }

        bool isInc = true;
        if (expiries->size() == 1) {
            return true;
        }

        DateTime lo = (*expiries)[0]->toDate(base);
        for (int i = 1; i < expiries->size() && isInc; i++) {
            DateTime hi = (*expiries)[i]->toDate(base);;
            isInc = hi > lo;
            lo = hi;
        }
        return isInc;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

class ExpiryHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Expiry, clazz);
        SUPERCLASS(CObject);
        clazz->enableCloneOptimisations(); // for derived types
        // no fields
    }

};

CClassConstSP const Expiry::TYPE = CClass::registerClassLoadMethod(
    "Expiry", typeid(Expiry), ExpiryHelper::load);

//template<> CClassConstSP const ExpiryArray::TYPE = CClass::registerClassLoadMethod(
//"ExpiryArray", typeid(ExpiryArray), load);
DEFINE_TEMPLATE_TYPE(ExpiryArray);

/** Addin for building handles to an expiry */
class ExpiryAddin: public CObject{
    static CClassConstSP const TYPE;

    /** the single parameter that our addin takes */
    ExpirySP  expiry;

    /** the 'addin function' - doesn't do anything clever but just allows us
        to use the infrastructure */
    static IObjectSP createExpiry(ExpiryAddin* params){
        if (!params->expiry){
            return CNull::create();
        }
        return params->expiry;
    }

    /** for reflection */
    ExpiryAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ExpiryAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultExpiryAddin);
        FIELD(expiry, "A benchmark date or an interval");
        FIELD_MAKE_OPTIONAL(expiry);
        Addin::registerClassObjectMethod(
            "EXPIRY",
            Addin::UTILITIES,
            "Constructs a handle to an expiry (either a benchmark date"
            " or an interval)",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)createExpiry);
    }

    static IObject* defaultExpiryAddin(){
        return new ExpiryAddin();
    }
    
};

CClassConstSP const ExpiryAddin::TYPE = CClass::registerClassLoadMethod(
    "ExpiryAddin", typeid(ExpiryAddin), load);

/** Addin for building handles to an expiry arrays */
class ExpiryArrayAddin: public CObject{
    static CClassConstSP const TYPE;

    /** the single parameter that our addin takes */
    ExpiryArraySP  expiryArray;

    /** the 'addin function' - doesn't do anything clever but just allows us
        to use the infrastructure */
    static IObjectSP createExpiryArray(ExpiryArrayAddin* params){
        if (!params->expiryArray){
            return CNull::create();
        }
        return params->expiryArray;
    }

    /** for reflection */
    ExpiryArrayAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ExpiryArrayAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultExpiryArrayAddin);
        FIELD(expiryArray, "An array of benchmark dates or intervals");
        FIELD_MAKE_OPTIONAL(expiryArray);
        Addin::registerClassObjectMethod(
            "EXPIRY_ARRAY",
            Addin::UTILITIES,
            "Constructs a handle to an expiry "
            "array (benchmark dates/intervals)",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)createExpiryArray);
    }

    static IObject* defaultExpiryArrayAddin(){
        return new ExpiryArrayAddin();
    }
    
};

CClassConstSP const ExpiryArrayAddin::TYPE = CClass::registerClassLoadMethod(
    "ExpiryArrayAddin", typeid(ExpiryArrayAddin), load);



//// Creates an Expiry object from a string.
ExpiryConstSP Expiry::createExpiry(const string& input) {
   string errorString = "Expiry::createExpiry(): Cannot create "
                        "Expiry from string '" + input + "'";
   ExpiryConstSP expiry;

   //Tokenize
   using namespace boost;
   char_delimiters_separator<char> sep(false, " "); 
   tokenizer<> tok(input, sep);
   vector<string> halves;
   std::copy(tok.begin(), tok.end(), back_inserter(halves));

   if (halves.size() == 2) { 
       //It's either a Benchmark date or a MaturityTimePeriod
       try {
           expiry = ExpiryConstSP(new MaturityTimePeriod(halves[0],  
                                                         DateTime::timeConvert(halves[1])));
       } catch (ModelException e) {
           // It wasn't a MaturityTimePeriod, it's a BenchmarkDate
           try {
              expiry = ExpiryConstSP(new BenchmarkDate(DateTime(halves[0], halves[1])));
           } catch (ModelException e) {
              throw ModelException(e, errorString);
           }
       }
   } else if (halves.size() == 1) {
       // It's a maturity period
       try { 
          expiry = ExpiryConstSP(new MaturityPeriod(halves[0]));
       } catch(ModelException e) {
          throw ModelException(e, errorString);
       }
   } else {
       //Don't know what this is.
       throw ModelException(errorString);
   }
   return expiry;
}



class TestCreateExpiry: public CObject{
 
public:
    static CClassConstSP const TYPE;
    string   input;
    DateTime date;
    
    string test(){
        return Expiry::createExpiry(input)->toDate(date).toString();
    }

    TestCreateExpiry(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(TestCreateExpiry, clazz);
        SUPERCLASS(CObject);
        FIELD(input, "input string");
        FIELD(date, "input date");
        EMPTY_SHELL_METHOD(defaultConstructor);
        Addin::registerStringMethod("EXPIRY_TEST",
                                    Addin::XL_TESTS,
                                    "Test createExpiry",
                                    &TestCreateExpiry::test);
    }

    static IObject* defaultConstructor(){
        return new TestCreateExpiry();
    }
};

CClassConstSP const TestCreateExpiry::TYPE = CClass::registerClassLoadMethod(
    "TestCreateExpiry", typeid(TestCreateExpiry), load);

DRLIB_END_NAMESPACE
