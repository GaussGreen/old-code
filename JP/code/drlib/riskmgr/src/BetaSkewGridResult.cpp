//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : BetaSkewGridResult.cpp
//
//   Description : Result for pointwise type beta skew sensitivities
//
//   Author      : Antoine Gregoire
//
//   Date        : 17-Jun-2005
//
//
//   $Log: $
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/BetaSkewGridResult.hpp"


DRLIB_BEGIN_NAMESPACE

BetaSkewGridResult::~BetaSkewGridResult(){}

BetaSkewGridResult::BetaSkewGridResult():
    CObject(TYPE),
    point(DateTime(), 0.0),
    result(0){}

BetaSkewGridResult::BetaSkewGridResult(const BetaSkewGridPoint& point, double result): 
    CObject(TYPE), point(point), result(result) {}

/** returns the grid point associated with this result */
BetaSkewGridPoint BetaSkewGridResult::getGridPoint() const{
    return point;
}

/** returns the double detailing this result */
double BetaSkewGridResult::getResult() const{
    return result;
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void BetaSkewGridResult::outputWrite(
    const string& linePrefix,
    const string& prefix,
    ostream& stream) const
{
    ios::fmtflags oldFlags     = stream.flags();     // save settings
    long          oldPrecision = stream.precision(); // save settings
    if (result < 1 && result > -1){
        // 8 dp
        stream.flags(ios::fixed);
    } else {
        // 8 sf
        stream.unsetf(ios::fixed);
    }        
    stream.precision(8);

    stream << linePrefix << prefix << "_" << point.getMaturity().toString() 
           << "_" << point.getStrike() << ": " << result << endl;
    stream.flags(oldFlags); // restore settings
    stream.precision(oldPrecision);
}

/** scale by factor x */
void BetaSkewGridResult::scale(double x) {
    result *= x;
}

/** add BetaSkewGridResult object to this result (Implementation of
    CombinableResult) */
void BetaSkewGridResult::add(const CombinableResult& x, double scaleFactor) {
    // gcc bug: force to IObject before dynamic cast
    const BetaSkewGridResult& ex = dynamic_cast<const BetaSkewGridResult&>(static_cast<const IObject&>(x));
    if (point.getMaturity() != ex.point.getMaturity() ||
        point.getStrike() != ex.point.getStrike())
    {
        // shouldn't happen really
        throw ModelException("BetaSkewGridResult::add", 
                             "can't add [" + 
                             ex.point.getMaturity().toString() + "," + 
                             Format::toString(ex.point.getStrike()) + "] to [" +
                             point.getMaturity().toString() + "," + 
                             Format::toString(point.getStrike()) + "]");
    }
    result += scaleFactor * ex.result;           
}

/** adds the supplied value to this object */
void BetaSkewGridResult::addToResult(double value){
    result += value;
}

/** Invoked when class is 'loaded' */
void BetaSkewGridResult::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BetaSkewGridResult, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(CombinableResult);
    EMPTY_SHELL_METHOD(defaultBetaSkewGridResult);
    FIELD(point, "Indicates which grid point was tweaked ");
    FIELD(result, "The result of the expiry tweak");
}

IObject* BetaSkewGridResult::defaultBetaSkewGridResult(){
    return new BetaSkewGridResult();
}

CClassConstSP const BetaSkewGridResult::TYPE = 
    CClass::registerClassLoadMethod(
        "BetaSkewGridResult",
        typeid(BetaSkewGridResult),
        BetaSkewGridResult::load);

DEFINE_TEMPLATE_TYPE(BetaSkewGridResultArray);

DRLIB_END_NAMESPACE
