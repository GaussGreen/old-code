//


#include "edginc/config.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/VegaPointwise.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/VegaProxyMatrix.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE
VegaMatrix::IShift::~IShift(){} // empty
VegaMatrix::IShift::IShift(){} // empty
VegaMatrix::IRestorableShift::IRestorableShift(){} // empty
VegaMatrix::IRestorableShift::~IRestorableShift(){} // empty

/** Sens Control for vega pointwise */
const string VegaMatrix::NAME = "VEGA_MATRIX";
const double VegaMatrix::DEFAULT_SHIFT = 0.001;

VegaMatrix::~VegaMatrix(){}

/** constructor with explicit shift size */
VegaMatrix::VegaMatrix(double     shiftSize):
    MatrixShift(TYPE, NAME, shiftSize) {
}

/** constructor with extra chilli sauce */
VegaMatrix::VegaMatrix(double     shiftSize,
                       IModel*    model,
                       Control*   control):
          MatrixShift(TYPE, NAME, shiftSize) {
    this->algorithm = model;
    this->control   = control;
}

/** the full monty */
VegaMatrix::VegaMatrix(
    double                     shiftSize,
    int                        expiryIdx,
    const ExpiryArrayConstSP&  allExpiries,
    int                        strikeIdx,
    const DoubleArraySP&       strikes) : MatrixShift(TYPE, NAME, shiftSize){
    this->expiryIdx = expiryIdx;
    this->expiries  = allExpiries;
    this->xLowerIdx = strikeIdx;
    this->xUpperIdx = strikeIdx;
    this->xAxis     = strikes;
}


/** for reflection */
VegaMatrix::VegaMatrix():
    MatrixShift(TYPE, NAME) {
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double VegaMatrix::divisor() const{
    // add check for 0 shift size - somewhere 
    return 100.0 * getShiftSize();
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP VegaMatrix::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP VegaMatrix::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool VegaMatrix::nameMatches(const OutputName&         name,
                             IObjectConstSP          obj){
    // cast obj to VegaMatrix::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(vegaObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void VegaMatrix::appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj){
    // cast obj to VegaMatrix::Shift and then invoke name method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(vegaObj.sensName(this)));
    namesList.push_back(outputName);
}

bool VegaMatrix::shift(IObjectSP obj) {
    // cast obj to VegaMatrix::Shift and then invoke shift method
    IShift& vegaObj =
        dynamic_cast<IShift&>(*obj);
    return vegaObj.sensShift(this);
}

void VegaMatrix::restore(IObjectSP obj) {
    // cast obj to VegaMatrix::Shift and then invoke restore method
    IRestorableShift& vegaObj =
        dynamic_cast<IRestorableShift&>(*obj);
    vegaObj.sensRestore(this);
}

/** implements a one sided scalar derivative for each instance of the
    market data which is sensitive to this SensControl */
void VegaMatrix::calculate(TweakGroup*      tweakGroup,
                           CResults*        results) 
{
    try{
        CDoubleArraySP sensitiveStrikes;
        ISensitiveStrikes* vegaMatrixImnt = 
            dynamic_cast<ISensitiveStrikes*>(tweakGroup->getInstrument());
        if (vegaMatrixImnt && !vegaMatrixImnt->avoidVegaMatrix(getModel())) {
            // calculate vega matrix
            MatrixShift::calculate(tweakGroup, results);
        } else {
            // do vega pointwise instead - save vegaPointwise object in field
            // to allow getMarketDataName to work correctly
            vegaPointwise.reset(new VegaPointwise(NAME, this->getShiftSize()));

            vegaPointwise->calculateSens(tweakGroup->getModel(),
                                         tweakGroup->getInstrument(),
                                         getControl(), 
                                         results);
            vegaPointwise.reset(); // tidy up
        }
    } catch (exception& e){
        vegaPointwise.reset();
        // don't kill everything else
        results->storeGreek(IObjectSP(new Untweakable(e)), 
                            getSensOutputName(),
                            OutputNameSP(new OutputName("")));
    }
}

/** returns the name identifying the market data to be shifted. Returns
    null if not set */
OutputNameConstSP VegaMatrix::getMarketDataName() const{
    // This value is used by models for optimisation purposes so when we do
    // vega pointwise rather than vega matrix make sure we return the right
    // information
    if (!vegaPointwise){
        return SensControlPerName::getMarketDataName();
    }
    return vegaPointwise->getMarketDataName();
}

DoubleArraySP VegaMatrix::getXAxis(CInstrument*      instrument, 
                                   OutputNameConstSP outputName)
{
    CDoubleArraySP sensStrikes;
    ISensitiveStrikes* vegaMatrixImnt =
        dynamic_cast<ISensitiveStrikes*>(instrument);

    if ( !vegaMatrixImnt ) {
        string m("Error: attempted to calculate sensitive strikes"
                 "for a product which is not derived from ISensitiveStrikes: "
                 + instrument->getClass()->getName());
        throw ModelException("VegaMatrix::getXAxis", m);
    }

    // Dummy try/catch to stop nt.vc6.opt wetting its pants

	try{
        sensStrikes = vegaMatrixImnt->getSensitiveStrikes(outputName, getModel());
	}
	catch (...) { throw; }

    // now strip duplicates
    sort(sensStrikes->begin(), sensStrikes->end());

    CDoubleArraySP unique(new DoubleArray(0));

    if (!sensStrikes->empty()) {
        unique->push_back((*sensStrikes)[0]);
    }

    for (int i = 1; i < sensStrikes->size(); i++) {
        if ((*sensStrikes)[i] > (*sensStrikes)[i-1] + DBL_EPSILON) {
            unique->push_back((*sensStrikes)[i]);
        }
    }

    return unique;
}

/** The supplied object is queried for the expiries array needed
    for doing vega pointwise and this array is returned. The supplied
    object must implement the VegaMatrix.Shift interface */
IObjectConstSP VegaMatrix::qualifier(IObjectConstSP obj){
    // cast obj to VegaMatrix::Shift and then invoke qualifier method
    const IShift& vegaObj =
        dynamic_cast<const IShift&>(*obj);
    return vegaObj.sensExpiries(this);
}

/** returns a new zero shift control */
MatrixShift* VegaMatrix::getZeroShift() const
{
    return (new VegaMatrix(0.0, getModel(), getControl()));
}

/** nasty */
VegaMatrix* VegaMatrix::fromProxy(VegaProxyMatrix* proxy) {
    VegaMatrixSP vm(new VegaMatrix(proxy->getShiftSize()));
    vm->expiries  = proxy->getExpiries();
    vm->xAxis     = CDoubleArraySP::constCast(proxy->getXAxisValues());
    vm->xLowerIdx = proxy->getLowerIdx();
    vm->xUpperIdx = proxy->getUpperIdx();
    vm->expiryIdx = proxy->getExpiryIdx();

    vm->setMarketDataName(proxy->getMarketDataName());   
    return vm.release();
}
 
class VegaMatrixHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new VegaMatrix(VegaMatrix::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new VegaMatrix(shiftSize);
        }
    };
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VegaMatrix, clazz);
        SUPERCLASS(MatrixShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultVegaMatrix);
        // one field
        FIELD_NO_DESC(vegaPointwise);
        FIELD_MAKE_TRANSIENT(vegaPointwise);
        // register how to build our sensitivity
        SensitivityFactory::addSens(VegaMatrix::NAME, 
                                    new Factory(), 
                                    new VegaMatrix(VegaMatrix::DEFAULT_SHIFT),
                                    VegaMatrix::IShift::TYPE);
    }

    static IObject* defaultVegaMatrix(){
        return new VegaMatrix();
    }
};

CClassConstSP const VegaMatrix::TYPE = CClass::registerClassLoadMethod(
    "VegaMatrix", typeid(VegaMatrix), VegaMatrixHelper::load);

CClassConstSP const VegaMatrix::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "VegaMatrix::IShift", typeid(VegaMatrix::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(VegaMatrix::IRestorableShift, clazz);
    EXTENDS(VegaMatrix::IShift);
}
    
CClassConstSP const VegaMatrix::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "VegaMatrix::IRestorableShift",
    typeid(VegaMatrix::IRestorableShift), 
    restorableShiftLoad);

CClassConstSP const ISensitiveStrikes::TYPE = 
CClass::registerInterfaceLoadMethod("ISensitiveStrikes",
                                    typeid(ISensitiveStrikes), 0);
ISensitiveStrikes::ISensitiveStrikes(){}
ISensitiveStrikes::~ISensitiveStrikes(){}

bool VegaMatrixLoad() {
    return VegaMatrix::TYPE != NULL && ISensitiveStrikes::TYPE != NULL;
}

DRLIB_END_NAMESPACE

