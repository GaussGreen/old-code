//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : ParSpreadMaxTweakSize.cpp
//
//   Description : Pseudo sensitivity that uses the qualifier mechanism
//                 to return maximum tweak sizes that may be applied to
//                 the par spread curve.
//
//   Author      : Gordon Stephens
//
//   Date        : 22 June 2005
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

ParSpreadMaxTweakSize::IShift::~IShift(){}

/** Returns this */
ITweakNameResolver* ParSpreadMaxTweakSize::nameResolver(){
    return this;
}

/** returns the name identifying the market data to be shifted. */
OutputNameConstSP ParSpreadMaxTweakSize::getMarketDataName() const{
    return name;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP ParSpreadMaxTweakSize::shiftInterface() const{
    return IShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    ParSpreadMaxTweakSize.Shift interface */
bool ParSpreadMaxTweakSize::nameMatches(const OutputName&  name,
                                        IObjectConstSP     obj)
{
    // cast obj to ParSpreadMaxTweakSize::IShift and then invoke name method
    const IShift& parSpreadMaxTweakSizeObj = dynamic_cast<const IShift&>(*obj);
    return name.equals(parSpreadMaxTweakSizeObj.sensName(this));
}

/** Casts supplied obj to IShift and calls sensMaxTweakSize() */
IObjectConstSP ParSpreadMaxTweakSize::qualifier(IObjectConstSP obj){
    const IShift& parSpreadMaxTweakSizeObj = dynamic_cast<const IShift&>(*obj);
    return parSpreadMaxTweakSizeObj.sensMaxTweakSize(this);
}

//-----------------------------
//ParSpreadMaxTweakSize methods
//-----------------------------
IExpiryRiskPropertyConstSP ParSpreadMaxTweakSize::wrt(){
    return withRespectTo;
}

ParSpreadMaxTweakSize::ParSpreadMaxTweakSize(
        OutputNameConstSP name,
        IExpiryRiskPropertyConstSP withRespectTo):
    name(name), withRespectTo(withRespectTo)
{}

ParSpreadMaxTweakSize::~ParSpreadMaxTweakSize(){}

/** Calculates appropriate shift sizes for each expiry given original
 * shift size */
DoubleArrayConstSP ParSpreadMaxTweakSize::calculateTweakSizes(
    const IObject*     tweakGroup,
    DoubleArrayConstSP origShiftSizes,
    OutputNameConstSP  name,
    const ExpiryArray& expiries) {
    static const string method = "ParSpreadMaxTweakSize::calculateTweakSizes";
    try{
        TRACE_METHOD;

        ASSERT(origShiftSizes->size() == expiries.size());

        int i;
        DoubleArraySP tweakSizes; // return value

        //get all implementing objects
        ObjectArrayConstSP maxTwkSizesObjs(SensMgrConst(tweakGroup).all(IShift::TYPE,this));
        if (!maxTwkSizesObjs) {
            // just use supplied shift
            tweakSizes = DoubleArraySP::constCast(origShiftSizes);
        } else {

            //will result in array of mktobj.expiries in length
            DoubleArrayConstSP maxTwks;

            //get the max tweak size for the first match
            //implemented via qualifier mechanism
            maxTwks = DoubleArrayConstSP::dynamicCast(qualifier((*maxTwkSizesObjs)[0]));
            //copy as the smallest so far
            DoubleArraySP maxTwkSizes(DoubleArraySP::constCast(maxTwks));

            //now compare the other matches
            for (i=1; i<maxTwkSizesObjs->size(); i++)
            {
                maxTwks = DoubleArrayConstSP::dynamicCast(qualifier((*maxTwkSizesObjs)[i]));

                //ensure array length consistency
                if (maxTwks->size() != maxTwkSizes->size())
                {
                    throw ModelException(method, "unexpected array length mismatch");
                }

                for (int j=0; j<maxTwks->size(); j++)
                {
                    //update maxTwkSizes where appropriate
                    if ((*maxTwks)[j] < (*maxTwkSizes)[j])
                    {
                        (*maxTwkSizes)[j] = (*maxTwks)[j];
                    }
                }
            }

            // rescale to be expiries.size
            ExpiryWindowArrayConstSP allExpiries = 
                withRespectTo->subjectQualifiers(
                    IObjectConstSP::attachToRef(tweakGroup),
                    name);
            tweakSizes.reset(new DoubleArray(expiries.size()));
            for (i=0; i<expiries.size(); i++) {
                int idx = expiries[i]->search(ExpiryWindow::expiries(allExpiries).get());
                // take smallest absolute value
                (*tweakSizes)[i] = Maths::min(fabs((*origShiftSizes)[i]),
                                              fabs((*maxTwkSizes)[idx]));
                // then preserve sign
                if ((*origShiftSizes)[i] < 0.0){
                    (*tweakSizes)[i] = -(*tweakSizes)[i]; 
                }
            }
        }
        // Validate all tweaksizes are large enough to prevent
        // numerical error swamping results
        for (i=0; i<tweakSizes->size(); i++) {
            if (fabs((*tweakSizes)[i]) < 1e-9) {
                throw ModelException(method, "Pointwise tweaking " +
                                     expiries[i]->toString() +
                                     " expiry results in bad CDS curve");
            }
        }
        return tweakSizes;
    } catch (ModelException& e) {
        throw ModelException(&e, method);
    }
}

CClassConstSP const ParSpreadMaxTweakSize::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "ParSpreadMaxTweakSize::IShift", typeid(ParSpreadMaxTweakSize::IShift), 0);

class ParSpreadMaxTweakSize_AdaptedRiskProperty:
       public CObject,
       public virtual IExpiryRiskProperty {

    static IObject* emptyShell() {
        return new ParSpreadMaxTweakSize_AdaptedRiskProperty(
                       IExpiryRiskPropertySP());
    }

    static void load(CClassSP& clazz) {
        REGISTER(ParSpreadMaxTweakSize_AdaptedRiskProperty, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IExpiryRiskProperty);
        EMPTY_SHELL_METHOD(emptyShell);
        FIELD(underlying, "underlying");
    }

    IExpiryRiskPropertyConstSP underlying;

public:

    static CClassConstSP const TYPE;

    ParSpreadMaxTweakSize_AdaptedRiskProperty(
            IExpiryRiskPropertyConstSP underlying):
        CObject(TYPE),
        underlying(underlying)
    {}

    bool discrete() const {
        return underlying->discrete();
    }

    CClassConstSP subjectInterface() const {
        return underlying->subjectInterface();
    }

    OutputNameArrayConstSP subjectNames(IObjectConstSP world) const {
        return underlying->subjectNames(world);
    }

    IRiskAxisConstSP axisFor(OutputNameConstSP subjectName,
                             QualifierConstSP qualifier) const {
        return underlying->axisFor(subjectName, qualifier);
    }

    QualifierArrayConstSP subjectQualifiers(IObjectConstSP world,
                                            OutputNameConstSP name) const {
        return underlying->subjectQualifiers(world, name);
    }

    DoubleArrayConstSP adaptiveCoefficients(
            MultiTweakGroupConstSP world,
            OutputNameConstSP name,
            ExpiryWindowArrayConstSP exps,
            double targetCoefficient) const {
        TRACE_METHOD;

        DoubleArrayConstSP targets = underlying->adaptiveCoefficients(
            world, name, exps, targetCoefficient);

        return ParSpreadMaxTweakSize(name, underlying).
                calculateTweakSizes(world.get(),
                                    targets,
                                    name,
                                    *ExpiryWindow::expiries(exps));
    }

    string toString() const {
        return underlying->toString() + " (with adaptive tweak size)";
    }
};

CClassConstSP const ParSpreadMaxTweakSize_AdaptedRiskProperty::TYPE = CClass::registerClassLoadMethod(
    "ParSpreadMaxTweakSize_AdaptedRiskProperty", typeid(ParSpreadMaxTweakSize_AdaptedRiskProperty), load);

IExpiryRiskPropertyConstSP ParSpreadMaxTweakSize::adapted(
        IExpiryRiskPropertyConstSP property) {
    return IExpiryRiskPropertyConstSP(
        new ParSpreadMaxTweakSize_AdaptedRiskProperty(property));
}

DRLIB_END_NAMESPACE
