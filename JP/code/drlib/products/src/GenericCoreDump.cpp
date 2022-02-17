//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericCoreDump.cpp
//
//   Description : Generic that cores in interesting ways
//
//   Author      : Andrew J Swain
//
//   Date        : 21 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/ClosedFormLN.hpp"
#include <signal.h>

DRLIB_BEGIN_NAMESPACE

class GenericCoreDump: public Generic1Factor,
                       public CClosedFormLN::IIntoProduct {
public:
    static CClassConstSP const TYPE; 

    virtual void Validate() {}

/** Invoked when Class is 'loaded' */
static void load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(GenericCoreDump, clazz);
    SUPERCLASS(Generic1Factor);
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultGenericCoreDump);
    FIELD(mode, "mode");
}

static IObject* defaultGenericCoreDump(){
    return new GenericCoreDump();
}
    
private:
    friend class GenericCoreDumpClosedForm;

    GenericCoreDump():Generic1Factor(TYPE) {}; 
    GenericCoreDump(const GenericCoreDump& rhs);
    GenericCoreDump& operator=(const GenericCoreDump& rhs);

    void price(Control* control, CResults* results)const{
        static const string method = "GenericCoreDump::price";
        
        try {
            if (mode == "SEGV") {
                raise(SIGSEGV);
            }
            else if (mode == "FPE") {
                raise(SIGFPE);
            }
            else if (mode == "ILL") {
                raise(SIGILL);
            }
            else if (mode == "CORE") {
                double* x = 0;
                x[0]++;
            }
            else if (mode == "FAIL") {
                throw ModelException(method, mode); 
            }
            else {
                raise(SIGTERM);
            }            
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
    

/** Implementation of ClosedFormLN::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

private:
    string mode;
};


CClassConstSP const GenericCoreDump::TYPE = CClass::registerClassLoadMethod(
    "GenericCoreDump", typeid(GenericCoreDump), GenericCoreDump::load);


/** private class */
class GenericCoreDumpClosedForm: public CClosedFormLN::IProduct{
private:
    const GenericCoreDump* gcd; // a reference

public:
    GenericCoreDumpClosedForm(const GenericCoreDump* gcd): gcd(gcd){}

    void price(CClosedFormLN* model,
               Control*    control, 
               CResults*   results)const{
        gcd->price(control, results);
    }
};

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* GenericCoreDump::createProduct(CClosedFormLN* model) const
{
    return new GenericCoreDumpClosedForm(this);
}

// for class loading 
bool GenericCoreDumpLoad() {
    return (GenericCoreDump::TYPE != 0);
}

DRLIB_END_NAMESPACE
