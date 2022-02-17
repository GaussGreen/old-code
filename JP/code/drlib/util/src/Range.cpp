//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Range.cpp
//
//   Description : 
//
//   Date        : 21 May 02
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_RANGE_CPP
#include "edginc/Range.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

void Boundary::load (CClassSP& clazz) {
	clazz->setPublic(); // make visible to EAS/spreadsheet
	REGISTER(Boundary, clazz);
	SUPERCLASS(CObject);
	EMPTY_SHELL_METHOD(defaultBoundary);
	
	FIELD(value,"Boundary value (if infinite set +/- 1 for +/- INFINITY)");
	FIELD(isinf, "Is the boundary infinite? [default = false]");
	FIELD_MAKE_OPTIONAL(isinf);
	FIELD(isclosed,"Is the boundary closed (If infinite will be set to false)? [default = true]")
	FIELD_MAKE_OPTIONAL(isclosed);
}

CClassConstSP const Boundary::TYPE = CClass::registerClassLoadMethod(
	"Boundary", 
	typeid(Boundary), 
	load);

DEFINE_TEMPLATE_TYPE(BoundaryArray);

/** private constructor */
Boundary::Boundary() : CObject(TYPE), isinf(false), isclosed(true)
{}

/** default 'constructor' */
IObject* Boundary::defaultBoundary() 
{
	return new Boundary();
}

// validate inputs immediately after the object is constructed
void Boundary::validatePop2Object()
{
	if(isinf)
	{
		isclosed = false;
		value = value > 0 ? 1.0 : -1.0;
	}
}

Boundary::Boundary(double value, bool isinf, bool isclosed):
CObject(TYPE),
value(value),
isinf(isinf),
isclosed(isclosed){}

Boundary::Boundary(const Boundary& rhs):
CObject(TYPE),
value(rhs.value),
isinf(rhs.isinf),
isclosed(rhs.isclosed){}

Boundary& Boundary::operator=(const Boundary& rhs){
    if (this == &rhs){
        return *this;
    }
    value = rhs.value;
    isinf = rhs.isinf;
    isclosed = rhs.isclosed;
    return *this;
}

/** Boundary classes */
bool Boundary::isGreater(const Boundary& rhs) const{
    if (isinf) {
        if (rhs.isinf) {
            return Maths::isPositive(value - rhs.value); /* value > rhs.value */
        }
        else if (value > 0.0) { /* sign = Plus */
            return true;
        }
        else {  /* sign = Minus */
            return false;
        }
    }
    else if (rhs.isinf) {
        if (rhs.value > 0.0) {  /* sign = Plus */
            return false;
        }
        else {  /* sign = Minus */
            return true;
        }
    }
    else {
        return Maths::isPositive(value - rhs.value); /* value > rhs.value */
    }
}

bool Boundary::isLess(const Boundary& rhs) const{
    if (isinf) {
        if (rhs.isinf) {
            return Maths::isNegative(value - rhs.value); /* value < rhs.value */
        }
        else if (value > 0.0) { /* sign = Plus */
            return false;
        }
        else {  /* sign = Minus */
            return true;
        }
    }
    else if (rhs.isinf) {
        if (rhs.value > 0.0) {  /* sign = Plus */
            return true;
        }
        else {  /* sign = Minus */
            return false;
        }
    }
    else {
        return Maths::isNegative(value - rhs.value); /* value < rhs.value */
    }
}

Boundary Boundary::subtract(double var) const{
    if (isinf){
        return Boundary(*this);
    }
    return Boundary(value - var, false, isclosed);
}

string Boundary::toString() const{
    string buffer;
    if (isinf){
        if (value > 0.0){
            buffer += "+";
        }
        else{
            buffer += "-";
        }
        buffer += "oo";
    }
    else{
        buffer += Format::toString("%.2f", value);
    }
    return buffer;
}

ClosedBoundary::ClosedBoundary(double value):
Boundary(value, false, true){}

OpenBoundary::OpenBoundary(const OpenBoundary& rhs): Boundary(rhs){}

OpenBoundary::OpenBoundary(double value):
Boundary(value, false, false){}

OpenBoundary::OpenBoundary(double sign, bool isinf):
Boundary(sign, isinf, false){}

Infinity::Infinity(Infinity::Sign sign):
OpenBoundary(static_cast<double>(sign), true){}

// RANGE ------------------------------------------------------------------------

void Range::load (CClassSP& clazz) {
	clazz->setPublic(); // make visible to EAS/spreadsheet
	REGISTER(Range, clazz);
	SUPERCLASS(CObject);
	EMPTY_SHELL_METHOD(defaultRange);

	FIELD(l, "Lower (left) boundary");
	FIELD(r, "Upper (right) boundary");
}

CClassConstSP const Range::TYPE = CClass::registerClassLoadMethod(
	"Range", 
	typeid(Range), 
	load);

DEFINE_TEMPLATE_TYPE(RangeArray);

/** private constructor */
Range::Range() : CObject(TYPE)
{}

/** default 'constructor' */
IObject* Range::defaultRange() 
{
	return new Range();
}


Range::~Range(){}

Range::Range(const Boundary& lower, const Boundary& upper):
CObject(TYPE),
l(lower), r(upper){}

Range::Range(double lower, bool lowerClosed,
             double upper, bool upperClosed):
CObject(TYPE),
    l(OpenBoundary(0)), r(OpenBoundary(0))
{
    if (lower == -HUGE_VAL) {
        l = Infinity(Infinity::Minus);
    }
    else if (lower == HUGE_VAL) {
        l = Infinity(Infinity::Plus);
    }
    else if (lowerClosed) {
        l = ClosedBoundary(lower);
    }
    else {
        l = OpenBoundary(lower);
    }

    if (upper == -HUGE_VAL) {
        r = Infinity(Infinity::Minus);
    }
    else if (upper == HUGE_VAL) {
        r = Infinity(Infinity::Plus);
    }
    else if (upperClosed) {
        r = ClosedBoundary(upper);
    }
    else {
        r = OpenBoundary(upper);
    }
}

Range::Range(const Range& rhs):
CObject(TYPE),
l(rhs.l), r(rhs.r){}

Range& Range::operator=(const Range& rhs){
    if (this == &rhs){
        return *this;
    }
    l = rhs.l;
    r = rhs.r;
    return *this;
}

bool Range::isNonEmpty() const{
    /* If l > r, the order is wrong */
    if (l.isGreater(r)){
        return false;
    }
    /* We now know l <= r. We must ensure that if l = r, then l and r are closed brackets */
    if (!r.isGreater(l)){    // l >= r and therefore l == r
        return (l.isClosedBracket() && r.isClosedBracket());
    }
    return true;
}

void Range::checkIsNonEmpty(const Range& interval){
    if (!interval.isNonEmpty()){
        throw ModelException("Range::checkIsNonEmpty",
                             interval.toString() + " does not define a non-empty interval");
    }
}

bool Range::isSingleton() const{
    /* If l > r, the order is wrong */
    if (l.isGreater(r)){
        return false;
    }
    /* If r > l, not a singleton */
    if (r.isGreater(l)){
        return false;
    }
    return (l.isClosedBracket() && r.isClosedBracket());
}

void Range::checkIsNotSingleton(const Range& interval){
    if (interval.isSingleton()){
        throw ModelException("Range::checkIsNotSingleton",
                             interval.toString() + " is a singleton");
    }
}

bool Range::isOpen() const{
    return (!l.isClosedBracket() && !r.isClosedBracket());
}

void Range::checkIsOpen(const Range& range){
    if (!range.isOpen()){
        throw ModelException("Range::checkIsOpen",
                             range.toString() + " is not open");
    }
}

bool Range::variableIsInRange(const Range& range, 
                              double var){
    const Boundary& l = range.l;
    const Boundary& r = range.r;
    // infinite interval
    if (l.isinf && r.isinf){  // (-infty, +infty)
        return true;
    }
    // semi infinite interval
    if (l.isinf){  // (-infty, upper)
        if (r.isclosed){
            return !Maths::isPositive(var - r.value);    // var <= r.value
        }
        return Maths::isPositive(r.value - var);    // var < r.value
    }
    if (r.isinf){  // (lower, +infty)
        if (l.isclosed){
            return !Maths::isPositive(l.value - var);    // l.value <= var
        }
        return Maths::isPositive(var - l.value);    // l.value < var
    }
    // finite interval
    if (l.isclosed){
        if (Maths::isPositive(l.value - var)){   // var < l.value
            return false;
        }
    }
    else if (!Maths::isPositive(var - l.value)){  // var <= l.value
        return false;
    }
    if (r.isclosed){
        if (Maths::isPositive(var - r.value)){   // var > r.value
            return false;
        }
    }
    else if (!Maths::isPositive(r.value - var)){  // var >= r.value
        return false;
    }
    return true;
}

void Range::checkVariableIsInRange(const Range& range, 
                                   double var,
                                   string varName){
    if (!variableIsInRange(range, var)){
        throw ModelException("Range::checkVariableIsInRange",
                             varName + " (" + Format::toString(var) + 
                             ") is not in range " + range.toString());
    }
}

Range Range::subtract(double var) const{
    return Range(l.subtract(var), r.subtract(var));
}

Range Range::subtract(const Range& range, double var){
    return range.subtract(var);
}

string Range::toString() const{
    string buffer;
    if (l.isClosedBracket()){
        buffer += "[";
    }
    else{
        buffer += "(";
    }
    buffer += l.toString();
    buffer += ", ";
    buffer += r.toString();
    if (r.isClosedBracket()){
        buffer += "]";
    }
    else{
        buffer += ")";
    }
    return buffer;
}

void Range::makeOpen(){
    l.isclosed = false;
    r.isclosed = false;
}

InfiniteRange::InfiniteRange():
Range(Infinity(Infinity::Minus), Infinity(Infinity::Plus)){}

RangeArraySP InfiniteRange::createInfiniteRangeArray(int size)
{ 
    RangeArraySP result( new RangeArray(size));
    for (int i=0; i<size; ++i)
    {
        ((*result)[i]).reset(new InfiniteRange());
    }
    return result;
}

DRLIB_END_NAMESPACE
