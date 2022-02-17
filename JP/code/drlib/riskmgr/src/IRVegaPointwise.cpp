//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRVegaPointwise.cpp
//
//   Description : Pointwise IR vega sensitivity
//
//   Author      : Andrew J Swain
//
//   Date        : 25 February 2002
//
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//        This class is now deprecated - Use IRVegaMatrix instead.
//        See IRVegaMatrix.hpp for more information
//   WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/SensitiveIRVolPoints.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Model.hpp"
#include "edginc/IRVegaPointwise.hpp"

DRLIB_BEGIN_NAMESPACE
IRVegaPointwise::IShift::~IShift(){} // empty
IRVegaPointwise::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega pointwise */
const string IRVegaPointwise::NAME = "IRVEGA_POINTWISE";
const double IRVegaPointwise::DEFAULT_SHIFT = 0.001;

/** constructor with explicit shift size */
IRVegaPointwise::IRVegaPointwise(double shiftSize):
    IRGridShift(TYPE, NAME, shiftSize) {}

/** for reflection */
IRVegaPointwise::IRVegaPointwise():
    IRGridShift(TYPE, NAME) {}

IRVegaPointwise::IRVegaPointwise(CClassConstSP clazz, 
                                 const string name):
    IRGridShift(clazz, name) {}

IRVegaPointwise::IRVegaPointwise(CClassConstSP clazz, 
                                 const string name, 
                                 const double shiftSize):
    IRGridShift(clazz, name, shiftSize) {}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double IRVegaPointwise::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP IRVegaPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP IRVegaPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool IRVegaPointwise::nameMatches(const OutputName&         name,
                                  IObjectConstSP          obj){
    // cast obj to IRVegaPointwise::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(vegaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void IRVegaPointwise::appendName(OutputNameArray&          namesList,
                                 IObjectConstSP          obj){
    // cast obj to IRVegaPointwise::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(vegaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool IRVegaPointwise::shift(IObjectSP obj) {
    // cast obj to IRVegaPointwise::Shift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

void IRVegaPointwise::restore(IObjectSP obj) {
    // cast obj to IRVegaPointwise::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

IRGridPointAbsArrayConstSP IRVegaPointwise::getGridPoints(
    CInstrument*      instrument, 
    OutputNameConstSP outputName) 
{
    static const string method("IRVegaPointwise::getGridPoints");
    try{
        const IModel* model = getModel();
        //  ask model
        const ISensitivePoints* vpModel =
        dynamic_cast<const ISensitivePoints*>(model);
        if (!vpModel) {
            string m("Error: attempted to calculate sensitive grid points for a"
                     "model which is not derived from "
                     "VegaPointwise::ISensitivePoints");
            throw ModelException(method, m);
        }
        IRGridPointAbsArraySP gridPoints(
            vpModel->getSensitiveIRVolPoints(outputName, instrument));
        if (!gridPoints){
            throw ModelException(method, "No IR grid points returned for "+
                                 outputName->toString());
        }
        if (!gridPoints->empty())
        {
            // then trim excess points if we can determine when to stop tweaking
            const LastSensDate* lsd = dynamic_cast<const LastSensDate*>(instrument);
            if (lsd){
                const DateTime& endDate = lsd->endDate(this);
                IRGridPointAbs::trim(*gridPoints, endDate);
            }
        }
        return gridPoints;
    } catch (exception& e){
        throw ModelException(e, method," For model "+
                             getModel()->getClass()->getName()+" with "
                             "instrument "+instrument->getClass()->getName());
    }
}

// should instrument bother doing this greek or not?
bool IRVegaPointwise::avoid(CInstrument* instrument,
                            IModel*      model) {
    ISensitivePoints* vp = dynamic_cast<ISensitivePoints*>(model);
    return (vp == 0);
}

class IRVegaPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new IRVegaPointwise(IRVegaPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new IRVegaPointwise(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IRVegaPointwise, clazz);
        SUPERCLASS(IRGridShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultIRVegaPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(IRVegaPointwise::NAME, 
                                    new Factory(), 
                                    new IRVegaPointwise(IRVegaPointwise::DEFAULT_SHIFT),
                                    IRVegaPointwise::IShift::TYPE);
    }

    static IObject* defaultIRVegaPointwise(){
        return new IRVegaPointwise();
    }
};

CClassConstSP const IRVegaPointwise::TYPE = CClass::registerClassLoadMethod(
    "IRVegaPointwise", typeid(IRVegaPointwise), IRVegaPointwiseHelper::load);

CClassConstSP const IRVegaPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "IRVegaPointwise::IShift", typeid(IRVegaPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(IRVegaPointwise::IRestorableShift, clazz);
    EXTENDS(IRVegaPointwise::IShift);
}
    
CClassConstSP const IRVegaPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "IRVegaPointwise::IRestorableShift",
    typeid(IRVegaPointwise::IRestorableShift), 
    restorableShiftLoad);

static void myLoad(CClassSP& clazz){
    REGISTER_INTERFACE(ISensitiveIRVolPoints, clazz);
    EXTENDS(IObject);
}

CClassConstSP const ISensitiveIRVolPoints::TYPE = 
CClass::registerInterfaceLoadMethod("ISensitiveIRVolPoints",
                                    typeid(ISensitiveIRVolPoints), myLoad);

IRVegaPointwise::ISensitivePoints::~ISensitivePoints(){}

DRLIB_END_NAMESPACE
