//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : MeanCreditReversionLevel.cpp
//
//   Description : Mean credit reversion level scenario - set CDSVol mean reversion to supplied value
//
//   Author      : Antoine Gregoire
//
//   Date        : September 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MeanCreditReversionLevel.hpp"

DRLIB_BEGIN_NAMESPACE

MeanCreditReversionLevel::IShift::~IShift(){} // empty

/** constructor with explicit spread level */
MeanCreditReversionLevel::MeanCreditReversionLevel(double spread): 
    ScalarPerturbation(TYPE, spread){}

/** for reflection */
MeanCreditReversionLevel::MeanCreditReversionLevel(): ScalarPerturbation(TYPE){
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP MeanCreditReversionLevel::shiftInterface() const{
    return IShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    IShift interface */
bool MeanCreditReversionLevel::nameMatches(const OutputName&         name,
                                 IObjectConstSP            obj){
    // cast obj to MeanCreditReversionLevel::IShift and then invoke name method
    const IShift& meanCreditReversionLevelObj = 
        dynamic_cast<const IShift&>(*obj);
    return name.equals(meanCreditReversionLevelObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's IShift interface */
void MeanCreditReversionLevel::appendName(OutputNameArray&          namesList,
                                IObjectConstSP            obj){
    // cast obj to MeanCreditReversionLevel::IShift and then invoke name method
    const IShift& meanCreditReversionLevelObj = 
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(
        new OutputName(meanCreditReversionLevelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool MeanCreditReversionLevel::shift(IObjectSP obj) {
    // cast obj to MeanCreditReversionLevel::IShift and then invoke shift method
    IShift& meanCreditReversionLevelObj = 
        dynamic_cast<IShift&>(*obj);
    return meanCreditReversionLevelObj.sensShift(this);
}

class MeanCreditReversionLevelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MeanCreditReversionLevel, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultMeanCreditReversionLevel);
        // no fields
    }

    static IObject* defaultMeanCreditReversionLevel(){
        return new MeanCreditReversionLevel();
    }
};

CClassConstSP const MeanCreditReversionLevel::TYPE = CClass::registerClassLoadMethod(
    "MeanCreditReversionLevel", typeid(MeanCreditReversionLevel), MeanCreditReversionLevelHelper::load);

CClassConstSP const MeanCreditReversionLevel::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "MeanCreditReversionLevel::IShift", typeid(MeanCreditReversionLevel::IShift), 0);

DRLIB_END_NAMESPACE
