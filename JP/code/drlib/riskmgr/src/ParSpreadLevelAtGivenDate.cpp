//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : ParSpreadLevelAtGivenDate.cpp
//
//   Description : CDS Par Spread level scenario
//                 1. Get spread S corresponding to given date (doing linear interpolation)
//                 2. Flatten spread curve to this value S
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ParSpreadLevelAtGivenDate.hpp"

DRLIB_BEGIN_NAMESPACE

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP ParSpreadLevelAtGivenDate::shiftInterface() const {
    return IShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    VegaParallel.Shift interface */
bool ParSpreadLevelAtGivenDate::nameMatches(const OutputName& name,
                                            IObjectConstSP    obj){
    return name.equals(dynamic_cast<const IShift &>(*obj).sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void ParSpreadLevelAtGivenDate::appendName(OutputNameArray& namesList, 
                                           IObjectConstSP        obj) {
    namesList.push_back(OutputNameSP(new OutputName(
        dynamic_cast<const IShift &>(*obj).sensName(this))));
}

/** Shifts the object (which supports being tweaked
    by this type of sens control) using given shift. The return value
    indicates whether or not components of this object need to be
    tweaked ie true: infrastructure should continue to recurse through
    components tweaking them; false: the infrastructure shouldn't
    touch any components within this object */
bool ParSpreadLevelAtGivenDate::shift(IObjectSP obj) {
    return dynamic_cast<IShift &>(*obj).sensShift(this);
}

/** Returns the level date */
DateTime ParSpreadLevelAtGivenDate::getLevelDate() const {
    return levelDate;
}

void ParSpreadLevelAtGivenDate::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ParSpreadLevelAtGivenDate, clazz);
    SUPERCLASS(Perturbation);
    EMPTY_SHELL_METHOD(defaultConstructor);
    
    FIELD(levelDate,
        "Date used to define the spread level used to flatten the par spreads curve");
}

ParSpreadLevelAtGivenDate::ParSpreadLevelAtGivenDate() : Perturbation(TYPE) {}

IObject* ParSpreadLevelAtGivenDate::defaultConstructor() {
    return new ParSpreadLevelAtGivenDate();
}

CClassConstSP const ParSpreadLevelAtGivenDate::TYPE =
CClass::registerClassLoadMethod("ParSpreadLevelAtGivenDate",
                                typeid(ParSpreadLevelAtGivenDate),
                                ParSpreadLevelAtGivenDate::load);

ParSpreadLevelAtGivenDate::IShift::~IShift() {}

CClassConstSP const ParSpreadLevelAtGivenDate::IShift::TYPE =
CClass::registerInterfaceLoadMethod("ParSpreadLevelAtGivenDate::IShift",
                                    typeid(ParSpreadLevelAtGivenDate::IShift),
                                    0);

DRLIB_END_NAMESPACE

