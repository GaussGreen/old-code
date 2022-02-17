//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : TimePoint2D.cpp
//
//   Description : Identifies a (DateTime, double, double) point
//
//   Author      : Antoine Gregoire
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TimePoint2D.hpp"

DRLIB_BEGIN_NAMESPACE

/** _TimePoint2DCompareDateC1C2 binary function */
bool TimePoint2D::CompareDateC1C2::operator()(const TimePoint2D& p1, const TimePoint2D& p2) const {
    return p1.lessThanDateC1C2(p2);
}

/** Comparison by 1) date 2) coord1 3) coord2 */
bool TimePoint2D::lessThanDateC1C2(const TimePoint2D& point) const {
    if (date == point.date) {
        if (coord1 == point.coord1) {
            return (coord2 < point.coord2);
        } else {
            return (coord1 < point.coord1);
        }
    } else {
        return (date < point.date);
    }
}

/** Comparison by 1) coord1 2) coord2 3) date */
bool TimePoint2D::lessThanC1C2Date(const TimePoint2D& point) const {
    if (coord1 == point.coord1) {
        if (coord2 == point.coord2) {
            return (date < point.date);
        } else {
            return (coord2 < point.coord2);
        }
    } else {
        return (coord1 < point.coord1);
    }
}

/** Access to date */
DateTime TimePoint2D::getDate() const {
    return date;
}

/** Access to first coordinate */
double TimePoint2D::getCoord1() const {
    return coord1;
}

/** Access to second coordinate */
double TimePoint2D::getCoord2() const {
    return coord2;
}

/** Virtual destructor */
TimePoint2D::~TimePoint2D() {}

/** Constructor */
TimePoint2D::TimePoint2D(const DateTime& date, double coord1, double coord2):
    CObject(TYPE), date(date), coord1(coord1), coord2(coord2) {}

/** Private constructor (for reflection) */
TimePoint2D::TimePoint2D(): CObject(TYPE) {}

/** Default constructor */
IObject* TimePoint2D::defaultTimePoint2D() {
    return new TimePoint2D();
}

/** Invoked when Class is 'loaded' */
void TimePoint2D::load(CClassSP& clazz) {
    clazz->setPrivate(); // don't make visible to EAS/spreadsheet
    REGISTER(TimePoint2D, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultTimePoint2D);
    FIELD_NO_DESC(date);
    FIELD_NO_DESC(coord1);
    FIELD_NO_DESC(coord2);
}

/** TYPE for TimePoint2D */
CClassConstSP const TimePoint2D::TYPE =
    CClass::registerClassLoadMethod(
        "TimePoint2D",
        typeid(TimePoint2D),
        load);

/** TYPE for TimePoint2DArray */
DEFINE_TEMPLATE_TYPE(TimePoint2DArray);

DRLIB_END_NAMESPACE

