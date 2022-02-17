//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Dividend.cpp
//
//   Description : Dividend representation
//
//   Author      : Stephen Hope
//
//   Date        : 25 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Dividend.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

// PUBLIC methods

const int Dividend::DIVIDEND_EXDIV_TIME = DateTime::EX_DIV_BEFORE_START_TIME;

/** Constructor */
Dividend::Dividend(const DateTime&     exDivDate,
                   const DateTime&     payDivDate,
                   const TDividendType divType,
                   const double        divAmount):
    CObject(TYPE), exDivDate(exDivDate), payDivDate(payDivDate), 
    divType(divType), divAmount(divAmount)
{
    validatePop2Object();
}

/** for reflection */
Dividend::Dividend(): CObject(TYPE), divType(AMOUNT), divAmount(0.0){
    // empty
}

/** Overridden for performance */
IObject* Dividend::clone() const{
    return new Dividend(*this);
}

/** Destructor Dividend::~Dividend() is in the header file for performance 
    reasons */

/* Access to data memebers */

/** return the ex-dividend date */
const DateTime& Dividend::getExDate()const
{
    return exDivDate;
}

/** return the payment date */
const DateTime& Dividend::getPayDate()const
{
    return payDivDate;
}

/** return the dividend amount */
double Dividend::getDivAmount()const
{
    return divAmount;
}

/** return the dividend type */
Dividend::TDividendType Dividend::getDivType()const
{
    return (TDividendType)divType;
}

/** set the dividend amount */
void Dividend::setDivAmount(double newAmount)
{
    divAmount = newAmount;
}

// PRIVATE methods

/** Validation checks 
    Conditions :
    ex-div date <= corresponding pay date
    dividend type = valid dividend type */
void Dividend::validatePop2Object() {
    static const string routine = "Dividend::validatePop2Object";
       
    if (exDivDate.getTime() != Dividend::DIVIDEND_EXDIV_TIME) {
        // force it
        // may be too slow - may have to modify original
        exDivDate = DateTime(exDivDate.getDate(),Dividend::DIVIDEND_EXDIV_TIME);
    }

    if (exDivDate.isGreater(payDivDate)) {
        throw ModelException(routine,
                             "exDivDate " + 
                             exDivDate.toString() +
                             " is after payment Date " + 
                             payDivDate.toString());
    }

    if (divType == Dividend::CONTINUOUS && 
        payDivDate.getTime() != Dividend::DIVIDEND_EXDIV_TIME) {
        // force this too
        payDivDate = DateTime(payDivDate.getDate(),Dividend::DIVIDEND_EXDIV_TIME);
    }

    if (divType != Dividend::AMOUNT &&
        divType != Dividend::PERCENT && 
        divType != Dividend::CONTINUOUS)
    {
        throw ModelException(routine,
                             "Unknown dividend type: (" +
                             Format::toString(divType) +
                             ") for dividend with exDivDate " +
                             exDivDate.toString() +
                             " and payDate " +
                             payDivDate.toString());
    }
}

/** tweak the dividend amount by shift size */
void Dividend::tweakDividend(double shiftSize){
    divAmount *= 1.0 + shiftSize;
}

/** Sets div type to 'dollar' (ie AMOUNT) and scales amount by scaling factor
    provided */
void Dividend::convertToDollar(double scalingFactor){
    divType = AMOUNT;
    divAmount *= scalingFactor;
}


void Dividend::convertToYield(double scalingFactor) {
    divType = PERCENT;
    divAmount *= scalingFactor;
}


/** Scaled dividend by supplied factor 
    ie divAmount -> divAmount* scalingFactor */
void Dividend::scale(double scalingFactor){
    divAmount *= scalingFactor;
}

/** Sets the pay date to supplied value */
void Dividend::setPayDate(const DateTime& payDate){
    payDivDate = payDate;
}

class DividendHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Dividend, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDividend);
        clazz->enableCloneOptimisations();
        // order the fields to reflect the usual order on the s/sheet
        FIELD(exDivDate, "ex-dividend date");
        FIELD(divAmount, "dividend amount");
        FIELD(divType, "dividend type");
        FIELD(payDivDate, "dividend payment date");
    }

    static IObject* defaultDividend(){
        return new Dividend();
    }
};

CClassConstSP const Dividend::TYPE = CClass::registerClassLoadMethod(
    "Dividend", typeid(Dividend), DividendHelper::load);

DEFINE_TEMPLATE_TYPE(DividendArray);

/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<Dividend>::toIObject(
    const Dividend& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<Dividend>::toIObject(Dividend& value){
    return DividendSP::attachToRef(&value);
}

/** Turns the IObjectSP into a Dividend */
const Dividend& arrayObjectCast<Dividend>::fromIObject(IObjectSP& value){
    Dividend *dtPtr = dynamic_cast<Dividend *>(value.get());
    if (!dtPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a Dividend");
    }
    return *dtPtr;
}

// explicit clone for arrays of Dividends - for performance
IObject* arrayClone<Dividend>::clone(const CArray* arrayToClone){
    const DividendArray& theArray = 
        static_cast<const DividendArray&>(*arrayToClone);
    return new DividendArray(theArray);
}

DRLIB_END_NAMESPACE


