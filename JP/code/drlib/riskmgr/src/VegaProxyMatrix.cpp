//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaProxyMatrix.cpp
//
//   Description : Controls calculation of fund proxy vega matrix
//
//   Author      : Andrew J Swain
//
//   Date        : 12 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/VegaProxyPointwise.hpp"
#include "edginc/VegaProxyMatrix.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"

DRLIB_BEGIN_NAMESPACE
VegaProxyMatrix::IShift::~IShift(){} // empty
VegaProxyMatrix::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega pointwise */
const string VegaProxyMatrix::NAME = "VEGA_PROXY_MATRIX";
const double VegaProxyMatrix::DEFAULT_SHIFT = VegaMatrix::DEFAULT_SHIFT;

/** constructor with explicit shift size */
VegaProxyMatrix::VegaProxyMatrix(double     shiftSize):
    MatrixShift(TYPE, NAME, shiftSize) {
}

/** for reflection */
VegaProxyMatrix::VegaProxyMatrix():
    MatrixShift(TYPE, NAME) {
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VegaProxyMatrix::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VegaProxyMatrix::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VegaProxyMatrix::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VegaProxyMatrix::nameMatches(const OutputName&         name,
                                  IObjectConstSP            obj){
    // cast obj to VegaProxyMatrix::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaObj.sensNameMatches(this, name);
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VegaProxyMatrix::appendName(OutputNameArray&          namesList,
                                 IObjectConstSP            obj){
    // cast obj to VegaProxyMatrix::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    vegaObj.sensAppendName(this, namesList);
}

bool VegaProxyMatrix::shift(IObjectSP obj) {
    // cast obj to VegaProxyMatrix::Shift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

    
void VegaProxyMatrix::restore(IObjectSP obj) {
    // cast obj to VegaProxyMatrix::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

/** implements a one sided scalar derivative for each instance of the
    market data which is sensitive to this SensControl */
void VegaProxyMatrix::calculate(TweakGroup*      tweakGroup,
                                CResults*        results) 
{
    try{
        CDoubleArraySP sensitiveStrikes;
        ISensitiveStrikes* VegaProxyMatrixImnt = 
            dynamic_cast<ISensitiveStrikes*>(tweakGroup->getInstrument());
        if (VegaProxyMatrixImnt && 
            !VegaProxyMatrixImnt->avoidVegaMatrix(getModel())) {
            // calculate vega matrix
            MatrixShift::calculate(tweakGroup, results);
        } else {
            // do proxy vega pointwise instead
            VegaProxyPointwiseSP vegaPointwise(
                new VegaProxyPointwise(NAME, this->getShiftSize()));

            vegaPointwise->calculateSens(tweakGroup->getModel(),
                                         tweakGroup->getInstrument(),
                                         getControl(), 
                                         results);
        }
    } catch (exception& e){
        // don't kill everything else
        results->storeGreek(IObjectSP(new Untweakable(e)), 
                            getSensOutputName(),
                            OutputNameSP(new OutputName("")));
    }
}

DoubleArraySP VegaProxyMatrix::getXAxis(CInstrument*      instrument, 
                                        OutputNameConstSP outputName)
{
    CDoubleArraySP sensStrikes;
    ISensitiveStrikes* VegaProxyMatrixImnt =
        dynamic_cast<ISensitiveStrikes*>(instrument);

    if ( !VegaProxyMatrixImnt ) {
        string m("Error: attempted to calculate sensitive strikes"
                 "for a product which is not derived from ISensitiveStrikes: "
                 + instrument->getClass()->getName());
        throw ModelException("VegaProxyMatrix::getXAxis", m);
    }
    sensStrikes = VegaProxyMatrixImnt->getSensitiveStrikes(outputName, getModel());
    return sensStrikes;
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the VegaProxyMatrix.Shift interface */
IObjectConstSP VegaProxyMatrix::qualifier(IObjectConstSP obj){
    // cast obj to VegaProxyMatrix::Shift and then invoke qualifier method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
   return vegaObj.sensExpiries(this);
}

/** returns a new zero shift control */
MatrixShift* VegaProxyMatrix::getZeroShift() const
{
    VegaProxyMatrixSP vega(new VegaProxyMatrix(0.0));
    vega->algorithm = getModel();
    vega->control   = getControl();
    return (vega.release());
}

 
class VegaProxyMatrixHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VegaProxyMatrix(VegaProxyMatrix::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VegaProxyMatrix(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VegaProxyMatrix, clazz);
        SUPERCLASS(MatrixShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVegaProxyMatrix);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(VegaProxyMatrix::NAME, 
                                    new Factory(), 
                                    new VegaProxyMatrix(VegaProxyMatrix::DEFAULT_SHIFT),
                                    VegaProxyMatrix::IShift::TYPE);
    }

    static IObject* defaultVegaProxyMatrix(){
        return new VegaProxyMatrix();
    }
};

CClassConstSP const VegaProxyMatrix::TYPE = CClass::registerClassLoadMethod(
    "VegaProxyMatrix", typeid(VegaProxyMatrix), VegaProxyMatrixHelper::load);

CClassConstSP const VegaProxyMatrix::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VegaProxyMatrix::IShift", typeid(VegaProxyMatrix::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VegaProxyMatrix::IRestorableShift, clazz);
    EXTENDS(VegaProxyMatrix::IShift);
}
    
CClassConstSP const VegaProxyMatrix::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VegaProxyMatrix::IRestorableShift",
    typeid(VegaProxyMatrix::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
