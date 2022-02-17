//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCStrata.cpp
//
//   Description : A definition of strata for stratified sampling
//
//   Date        : 22 Nov 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_QMCSTRATA_CPP
#include "edginc/QMCStrata.hpp"

DRLIB_BEGIN_NAMESPACE

CClassConstSP const QMCStrata::TYPE = CClass::registerClassLoadMethod(
    "QMCStrata", typeid(QMCStrata), load);

void QMCStrata::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(QMCStrata, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(QMCStrata::defaultConstructor);

    FIELD(name, "Name for this strata");
    FIELD_MAKE_OPTIONAL(name); // default name empty string
}

bool QMCStrata::notIntersecting(QMCStrataConstSP strata2) const
{
    checkTypes(strata2);
    return false; // undefined strata (i.e. full \Omega) intersects everything else
}


void QMCStrata::checkTypes(QMCStrataConstSP strata2) const
{
    QLIB_VERIFY(this->getClass() == strata2->getClass(), // must be of the same type
        "Presently all the stratas in stratification must be of identically same type "
        "providing stratification of one and only one asset.");
}


//DEFINE_TEMPLATE_TYPE(QMCStrataArray);


CClassConstSP const QMCStrataCIDJumps::TYPE = CClass::registerClassLoadMethod(
    "QMCStrataCIDJumps", typeid(QMCStrataCIDJumps), load);

void QMCStrataCIDJumps::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(QMCStrataCIDJumps, clazz);
    SUPERCLASS(QMCStrata);
    EMPTY_SHELL_METHOD(QMCStrataCIDJumps::defaultConstructor);

	//// SPREADSHEET FILLED-IN FIELDS +++++++++++++++++++++++++++++++++++++++++++
    FIELD(lowestNJumps, "Lowest number of jumps in this strata");
    FIELD(highestNJumps,"Highest number of jumps in this strata (infinity if skipped)");
    FIELD_MAKE_OPTIONAL(highestNJumps); // default name empty string
}

//DEFINE_TEMPLATE_TYPE(QMCStrataCIDJumpsArray);


bool QMCStrataCIDJumps::notIntersecting(QMCStrataConstSP strata2) const
{
    checkTypes(strata2);
    const QMCStrataCIDJumps* p2 = dynamic_cast<const QMCStrataCIDJumps*>(strata2.get());
    QLIB_VERIFY(p2, 
        "Internal... Presently all the stratas in stratification must be of identically same type "
        "providing stratification of one and only one asset.");

    return  (this->lowestNJumps < p2->lowestNJumps && this->highestNJumps <= p2->lowestNJumps) ||
            (this->lowestNJumps >= p2->highestNJumps && this->highestNJumps > p2->highestNJumps);           
}


CClassConstSP const QMCStratificationRecord::TYPE = CClass::registerClassLoadMethod(
    "QMCStratificationRecord", typeid(QMCStratificationRecord), load);

void QMCStratificationRecord::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(QMCStratificationRecord, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(QMCStratificationRecord::defaultConstructor);

    FIELD(name, "Name for this strata");
    FIELD_MAKE_OPTIONAL(name); // default name empty string

	//// SPREADSHEET FILLED-IN FIELDS +++++++++++++++++++++++++++++++++++++++++++
    FIELD(strata, "The attached strata");
    FIELD(weight, "The weight of the attached strata");
}

DEFINE_TEMPLATE_TYPE(QMCStratificationRecordArray);


DRLIB_END_NAMESPACE
