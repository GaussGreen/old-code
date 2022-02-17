//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParamSens.hpp
//
//   Description : Handles sensitivities for individual fields
//
//   Author      : Mark A Robson
//
//   Date        : 22 Nov 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/ScalarShift.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/TweakFunction.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/Results.hpp"
#include "edginc/FieldPath.hpp"
#include "edginc/fieldRiskProperty.hpp"
#include "edginc/IFieldTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/** Handles sensitivities for individual fields. Typically these are
    fields in parameterised vols. Some way of controlling which fields
    can be tweaked is needed and this class uses the
    Calibrator::IAdjustable interface to do this.
    It assumes that the tweaks are 'additive' can be summed for composites */
class ParamSens: public Sensitivity,
                 public virtual Additive{
    /// fields ///
    string       typeToTweak;   // identifies type of object to tweak
    double       shiftSize;     // what to shift by (same for all fields)
    StringArray  fieldsToTweak; // [optional] which fields to tweak
public:
    static CClassConstSP const TYPE;
    static const string NAME; //PARAM_SENS

    /** identifies the name used for storing associated results in the output*/
    const string& getSensOutputName() const{
        return NAME;
    }
    
    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one (return value: false) */
    bool discreteShift() const{
        return false;
    }
    //// constructor for the case of one field
    ParamSens(const string& typeToTweak, double shiftSize, 
              const string& fieldToTweak): 
        Sensitivity(TYPE), typeToTweak(typeToTweak), shiftSize(shiftSize),
        fieldsToTweak(1, fieldToTweak){}

    /** calculates given sensitivity and returns the names calculated.
        Implementation below */
    OutputNameArrayConstSP calculateAndReturnNames(Control*         control,
                                                   TweakGroup*      tweakGroup,
                                                   Results*         results);

private:
    ParamSens(): Sensitivity(TYPE){}

    static IObject* defaultMaker(){
        return new ParamSens();
    }
    /** Can't really offer much default functionality */
    class Factory: public virtual SensitivityFactory::ICreation{
    public:
        Factory(){}
    };
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ParamSens, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultMaker);
        FIELD(typeToTweak, "Type of object whose fields "
                     "are to be tweaked");
        FIELD(shiftSize, "Amount to shift parameters by");
        FIELD(fieldsToTweak, "Fields of object to be tweaked");
        FIELD_MAKE_OPTIONAL(fieldsToTweak);
        // register how to build our sensitivity
        SensitivityFactory::addSens(ParamSens::NAME, 
                                    new Factory(), 
                                    new ParamSens(),
                                    Calibrator::IAdjustable::TYPE);
    }

    /** Return the fields (as CFields) that we should be tweaking  */
    CFieldArray getFields(const CClassConstSP& clazz){
        static const string method("ParamSens::getFields");
        CFieldArray fields;
        if (!fieldsToTweak.empty()){
            // to do: move block somewhere central
            fields.reserve(fieldsToTweak.size());
            for (int i = 0; i < fieldsToTweak.size(); i++){
                const string& field = fieldsToTweak[i];
                CClassConstSP c = clazz;
                CFieldConstSP realField;
                do {
                    realField = c->hasDeclaredField(field);
                } while (!realField && (c = c->getSuperClass()) != 0);
                if (!realField){
                    throw ModelException(method, "Field '"+
                                         field + "' not "
                                         "found in "+ typeToTweak);
                }
                fields.push_back(realField);
            }
        } else {
            fields = Calibrator::IAdjustable::getFields(clazz);
        }
        if (fields.empty()){
            throw ModelException(method, "No fields to tweak!");
        }
        return fields;
    }

protected:
    /** calculates given sensitivity - invoked by calculateSens. */
    virtual void calculate(TweakGroup*      tweakGroup,
                           Results*         results){
        calculateAndReturnNames(getControl(), tweakGroup, results);
    }
public:
    class Ctrl;
    class Output;
};


/** Class used to return outputs. For clarity I think this class is worth
    having rather than return arrays of objects etc. It also makes the 
    handling needed for composites easier */
class ParamSens::Output: public CObject,
                 public virtual CombinableResult{
public:
    string       type; // identifies type of object which was tweaked
    StringArray  fieldNames; // identifies what was tweaked
    ObjectArray  results;   // Result for the tweaking of each field
    static CClassConstSP const TYPE;
    Output(const string&      type): CObject(TYPE), type(type){}

    /** Appends result to internal array */
    void addResult(const string& field, const IObjectSP& result){
        fieldNames.push_back(field);
        results.push_back(result);
    }

    /** scale by factor x */
    virtual void scale(double x){
        for (int i = 0; i < results.size(); i++){
            const IObjectSP& obj = results[i];
            CombinableResult* cr = 
                dynamic_cast<CombinableResult*>(obj.get());
            if (!cr){
                throw ModelException("ParamSens::Output::scale",
                                     "Could not scale object of type "+
                                     obj->getClass()->getName());
            }
            cr->scale(x);
        }
    }

    /** add an object (scaled by scaleFactor) to this
        result. Implementations should modify this result. If the x is
        not the same type as this then a [class cast] exception will
        be thrown */
    virtual void add(const CombinableResult& x, double scaleFactor){
        static const string method("ParamSens::Output::add");
        const Output& output = dynamic_cast<const Output&>(static_cast<const IObject&>(x));
        // we fail here if we try and combine things that don't match.
        // May have to review this before it goes into Pyramid
        if (output.type != type){
            throw ModelException(method, "Cannot combine results of different"
                                 " type ("+output.type+" and "+type+")");
        }
        if (fieldNames.size() != output.fieldNames.size()){
            throw ModelException(method, "Cannot combine results for different"
                                 " sets of fields");
        }
        for (int i = 0; i < results.size(); i++){
            const IObjectSP& obj1 = results[i];
            const IObjectConstSP& obj2 = output.results[i];
            CombinableResult* cr1 = 
                dynamic_cast<CombinableResult*>(obj1.get());
            const CombinableResult* cr2 = 
                dynamic_cast<const CombinableResult*>(obj2.get());
            if (!cr1 || !cr2 || !obj1->getClass()->isInstance(obj2)){
                throw ModelException(method,
                                     "Could not combine objects of type "+
                                     obj1->getClass()->getName()+ " and "+
                                     obj1->getClass()->getName());
            }
            cr1->add(*cr2, scaleFactor);
        }
    }
private:
    Output(): CObject(TYPE){}
    static IObject* defaultMaker(){
        return new Output();
    }
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(Output, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultMaker);
        IMPLEMENTS(CombinableResult);
        FIELD(type, "identifies type of object which was tweaked");
        FIELD(fieldNames, "identifies what was tweaked");
        FIELD(results, "Result for the tweaking of each field");
    }
};

CClassConstSP const ParamSens::Output::TYPE = CClass::registerClassLoadMethod(
    "ParamSens::Output", typeid(ParamSens::Output), load);

/* Onto the interesting part of the class. This is what does the
   actual work */
class ParamSens::Ctrl: public ScalarShift{
private:
    CFieldConstSP field; // what we're tweaking $unregistered
    Range         validRange; // permitted values of the field $unregistered
public:
    static CClassConstSP const TYPE;

    /** Constructor */
    Ctrl(double               shiftSize, 
         const CFieldConstSP& field):
        ScalarShift(TYPE, ParamSens::NAME, shiftSize), field(field), 
        validRange(Calibrator::IAdjustable::getRange(field)){}
    
    /** Override as fields don't support clone */
    IObject* clone(){
        Ctrl* ctrl = new Ctrl(getShiftSize(), field);
        return ctrl;
    }

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    virtual double divisor() const{
        // avoid arbitrary scaling as the fields will represent different 
        // things
        return getShiftSize();
    }

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    virtual CClassConstSP shiftInterface() const{
        return field->getDeclaringClass();
    }

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    virtual CClassConstSP restorableShiftInterface() const{
        return field->getDeclaringClass();
    }

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity.  The object must implement the
        shiftInterface() interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP          obj){
        const Calibrator::IAdjustable& object = 
            dynamic_cast<const Calibrator::IAdjustable&>(*obj);
        return name.equals(object.getName());
    }

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj){
        const Calibrator::IAdjustable& object = 
            dynamic_cast<const Calibrator::IAdjustable&>(*obj);
        OutputNameSP outputName(new OutputName(object.getName()));
        namesList.push_back(outputName);
    }

    /** Shifts the object (which supports being tweaked
        by this type of sens control) using given shift. The return value
        indicates whether or not components of this object need to be
        tweaked ie true: infrastructure should continue to recurse through
        components tweaking them; false: the infrastructure shouldn't
        touch any components within this object */
    virtual bool shift(IObjectSP obj){
        // save the current value
        double val = field->getDouble(obj);
        setInitialValue(val);
        val += getShiftSize();
        // check value is still in range
        Range::checkVariableIsInRange(validRange, val, field->getName());
        // then set it
        field->setDouble(obj, val);
        // and tell object that we've changed something
        obj->fieldsUpdated(CFieldArray(1, field)); 
        return true;
    }

    /** Restores the object (which supports being tweaked
        by this type of sens control) to its original form */
    virtual void restore(IObjectSP obj){
        // get the current value
        double val = getInitialValue();
        // then set it
        field->setDouble(obj, val);
        // and tell object that we've changed something
        obj->fieldsUpdated(CFieldArray(1, field)); 
    }

private:
    static void load(CClassSP& clazz){
        REGISTER(Ctrl, clazz);
        SUPERCLASS(ScalarShift);
        // don't support external creation
    }
};

CClassConstSP const ParamSens::Ctrl::TYPE = CClass::registerClassLoadMethod(
    "ParamSens::Ctrl", typeid(ParamSens::Ctrl), load);

/** calculates given sensitivity - invoked by calculateSens */
OutputNameArrayConstSP ParamSens::calculateAndReturnNames(
    Control*         control,
    TweakGroup*      tweakGroup,
    Results*         results){
    static const string method("ParamSens::calculateAndReturnNames");
    try{
        // check specified type supports Calibrator::IAdjustable
        CClassConstSP clazz = CClass::forName(typeToTweak);
        if (!Calibrator::IAdjustable::TYPE->isAssignableFrom(clazz)){
            throw ModelException(method, "Object of type "+typeToTweak+" do "
                                 "not support parameterised tweaking");
        }
        // then get hold of fields to tweak
        CFieldArray fields(getFields(clazz));
        // look up price up front
        double price;
        CInstrument* inst = tweakGroup->getInstrument();
        price = getSensPrice(results, inst, tweakGroup->getModel(), control);

        const string& ccy = results->getCcyName();
        OutputNameArrayConstSP namesToDo; // get first time through
        if (hasOverrideNames()){
            // unless they've been specified
            namesToDo = overrideNames();
        }

        int iName = 0;
        // loop across names to tweak (except we don't them yet)...
        do {
            // set up space for results
            smartPtr<Output> allResults(new Output(typeToTweak));
            for (unsigned int i = 0; i < fields.size(); i++){
                try{
                    // check that this field supports being 'calibrated'
                    Calibrator::IAdjustable::getRange(fields[i]);
                } catch (exception&){
                    // ignore exception - would probably confuse matters
                    throw ModelException(method, "Field "+fields[i]->getName()+
                                         " does not support being tweaked");
                }
                //now build up the SensControl which is going to do the work
                smartPtr<Ctrl> sensCtrl(new Ctrl(shiftSize, fields[i]));
                if (!namesToDo || namesToDo->empty()){
                    if (!namesToDo){
                        // Get hold of the names of what we'll be tweaking
                        namesToDo = sensCtrl->names(tweakGroup);
                    }
                    if (namesToDo->empty()) {
                        results->storeNotApplicable(this);
                        return namesToDo;
                    }
                }
                const OutputNameSP& name = (*namesToDo)[iName];
                // do one name at a time (easier for result handling)
                OutputNameArraySP oneName(new OutputNameArray(1, name));
                sensCtrl->storeOverrideNames(oneName);
                // create new Results object to allow us to store in our form
                Results tmpResults;
                tmpResults.storePrice(price, ccy);
                // go through calculateSens as it sets up SensCtrl properly
                sensCtrl->calculateSens(tweakGroup->getModel(),
                                        tweakGroup->getInstrument(),
                                        control,
                                        &tmpResults);
                // then get hold of and store the result
                IObjectConstSP greek(tmpResults.retrieveGreek(NAME, name));
                allResults->addResult(fields[i]->getName(), 
                                      IObjectSP(greek->clone()));
            }
            // then store allResults in main results object
            results->storeGreek(allResults, NAME, (*namesToDo)[iName]);
            iName++;
        } while (iName < namesToDo->size());
        return namesToDo;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

const string ParamSens::NAME = "PARAM_SENS";

CClassConstSP const ParamSens::TYPE = CClass::registerClassLoadMethod(
    "ParamSens", typeid(ParamSens), ParamSens::load);

/** Class offering alternative view onto ParamSens - in particular allows
    definitions of fields etc to be passed in via constructor ie they
    can be hard coded */
class ParamSensSingleField: public Sensitivity,
                            public virtual Additive,
                            public virtual IPerturbation {
    double shiftSize;
    string sensOutputName; // transient 
    string typeToTweak; // transient 
    string fieldToTweak; // transient 
public:
    static CClassConstSP const TYPE;

    /** identifies the name used for storing associated results in the output*/
    virtual const string& getSensOutputName() const{
        return sensOutputName;
    }

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one (return value: false) */
    virtual bool discreteShift() const{
        return false;
    }

    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*      tweakGroup,
                           Results*         results){
        // build relevant instance of ParamSens
        ParamSens paramSens(typeToTweak, shiftSize, fieldToTweak);
        // create empty results except for price
        Results tmpResults;
        tmpResults.storePrice(results->retrievePrice(), results->getCcyName());
        OutputNameArrayConstSP names(paramSens.
                                     calculateAndReturnNames(getControl(),
                                                             tweakGroup,
                                                             &tmpResults));
        if (names->empty()){
            results->storeNotApplicable(this);
        } else {
            for (int i = 0; i < names->size(); i++){
                const OutputNameSP& name = (*names)[i]; // for ease
                IObjectConstSP singleResult = 
                    tmpResults.retrieveGreek(paramSens.getSensOutputName(),
                                             name);
                const ParamSens::Output& output = 
                    dynamic_cast<const ParamSens::Output&>(*singleResult);
                if (output.results.size() != 1){
                    throw ModelException("ParamSensSingleField::calculate",
                                         "Internal error");
                }
                results->storeGreek(output.results.front(),
                                    sensOutputName, name);
            }
        }
    }

    /** IPerturbation implementation */
    bool findAndShift(IObjectSP objectToShift,
                      OutputNameConstSP name) {
        try {
            IScalarRiskPropertyConstSP property = fieldRiskProperty::scalar(
                CClass::forName(typeToTweak),
                IFieldTweak::bulk(
                    FieldPath::SP(fieldToTweak),
                    FieldPathConstSP(),
                    CDouble::SP(shiftSize),
                    IFieldTweak::IOperator::numeric(
                        TweakFunction::additive(),
                        InfiniteRange(),
                        true,    // useCalibratorRange, overriding InfiniteRange
                        false)), // clip
                true);           // absoluteDistance

            IHypothesisConstSP hyp = property->axisFor(name)->hypothesis(1.);
            // applyScenario should be a const method
            return IHypothesisSP::constCast(hyp)->applyScenario(objectToShift);
        }
        catch (exception& e) {
            throw ModelException(e, "ParamSensSingleField::findAndShift()");
        }
    }

protected:
    ParamSensSingleField(CClassConstSP clazz, 
                         const char* _sensOutputName,
                         const char* _typeToTweak,
                         const char* _fieldToTweak,
                         double defaultShiftSize): 
        Sensitivity(clazz), sensOutputName(_sensOutputName), 
        typeToTweak(_typeToTweak), fieldToTweak(_fieldToTweak), shiftSize(defaultShiftSize) {}
private:
    static void load(CClassSP& clazz){
        REGISTER(ParamSensSingleField, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        IMPLEMENTS(IPerturbation);
        FIELD(shiftSize, "Amount to shift parameters by");
        FIELD_NO_DESC(sensOutputName);
        FIELD_MAKE_TRANSIENT(sensOutputName);
        FIELD_NO_DESC(typeToTweak);
        FIELD_MAKE_TRANSIENT(typeToTweak);
        FIELD_NO_DESC(fieldToTweak);
        FIELD_MAKE_TRANSIENT(fieldToTweak);
    }
};                                      
                
CClassConstSP const ParamSensSingleField::TYPE = 
CClass::registerClassLoadMethod(
    "ParamSensSingleField", typeid(ParamSensSingleField), load);
#if 0           
class VHTRefSkew: public ParamSensSingleField{
public:
    static CClassConstSP const TYPE;

    VHTRefSkew():
        ParamSensSingleField(TYPE, "VHT_REF_SKEW", "VolHyperTrig", "skew"){}
private:
    static IObject* defaultConstructor(){
        return new VHTRefSkew();
    }
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VHTRefSkew, clazz);
        SUPERCLASS(ParamSensSingleField);
        EMPTY_SHELL_METHOD(defaultConstructor);
        // no fields
    }
};
CClassConstSP const VHTRefSkew::TYPE = CClass::registerClassLoadMethod(
    "VHTRefSkew", typeid(VHTRefSkew), load);
#endif

//// macro definition to allow us to create identical classes but with
//// different strings in them. Not sure how to do it with templates

#define GENERIC_PARAM_SENS_SINGLE_FIELD(className, sensOutputName, \
                                        typeToTweak, fieldToTweak, defaultShift) \
class className: public ParamSensSingleField{ \
public: \
    static CClassConstSP const TYPE; \
    const static double DEFAULT_SHIFT; \
 \
    className(double _shiftSize): \
        ParamSensSingleField(TYPE, sensOutputName,  \
                             typeToTweak, fieldToTweak, _shiftSize){} \
    \
    class Factory: public SensitivityFactory::IDefault, \
        public SensitivityFactory::IScalar {    \
    public: \
        virtual Sensitivity* createDefault(){   \
            return new className(className::DEFAULT_SHIFT);   \
        }   \
        virtual Sensitivity* createScalar(double shiftSize){    \
            return new className(shiftSize);       \
        }   \
    };  \
private: \
        className(): \
        ParamSensSingleField(TYPE, sensOutputName,  \
                typeToTweak, fieldToTweak, className::DEFAULT_SHIFT){}  \
    static IObject* defaultConstructor(){ \
        return new className(); \
    } \
    static void load(CClassSP& clazz){ \
        clazz->setPublic(); /* make visible to EAS/spreadsheet */ \
        REGISTER(className, clazz); \
        SUPERCLASS(ParamSensSingleField); \
        EMPTY_SHELL_METHOD(defaultConstructor); \
        SensitivityFactory::addSens(sensOutputName, \
            new Factory(),  \
            new className,  \
            NULL);  \
        /* no fields */ \
    } \
}; \
 \
const double className::DEFAULT_SHIFT = defaultShift;   \
CClassConstSP const className::TYPE = CClass::registerClassLoadMethod( \
#className, typeid(className), load);

GENERIC_PARAM_SENS_SINGLE_FIELD(VHTRefSkew, "VHT_REF_SKEW", 
                                "VolHyperTrig", "skew", -0.0005);
GENERIC_PARAM_SENS_SINGLE_FIELD(VHTRefConv, "VHT_REF_CONV", 
                                "VolHyperTrig", "convexity", 0.0001);
GENERIC_PARAM_SENS_SINGLE_FIELD(VHTRefSkewCelerity, "VHT_REF_SKEW_CELERITY", 
                                "VolHyperTrig", "skewCelerity", 0.05);
GENERIC_PARAM_SENS_SINGLE_FIELD(VHTRefConvCelerity, "VHT_REF_CONV_CELERITY", 
                                "VolHyperTrig", "convexityCelerity", 0.05);

// VolSV
GENERIC_PARAM_SENS_SINGLE_FIELD(SVInitialVol, "SV_INITIAL_VOL", 
                                "VolSV", "initialVol", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVMeanVol, "SV_MEAN_VOL", 
                                "VolSV", "meanVol", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVMeanReversRate, "SV_MEAN_REVERS_RATE", 
                                "VolSV", "meanReversRate", 0.05);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVVolVol, "SV_VOL_VOL", 
                                "VolSV", "volVol", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCorrelation, "SV_CORRELATION", 
                                "VolSV", "correlation", -0.01);

//VolSVJ
GENERIC_PARAM_SENS_SINGLE_FIELD(SVJInitialVol, "SVJ_INITIAL_VOL", 
                                "VolSVJ", "initialVol", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVJMeanVol, "SVJ_MEAN_VOL", 
                                "VolSVJ", "meanVol", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVJMeanReversRate, "SVJ_MEAN_REVERS_RATE", 
                                "VolSVJ", "meanReversRate", 0.05);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVJVolVol, "SVJ_VOL_VOL", 
                                "VolSVJ", "volVol", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVJCorrelation, "SVJ_CORRELATION", 
                                "VolSVJ", "correlation", -0.01);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVJCrashRate, "SVJ_CRASH_RATE", 
                                "VolSVJ", "crashRate", 0.01);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVJCrashSizeMean, "SVJ_CRASH_SIZE_MEAN", 
                                "VolSVJ", "crashSizeMean", -0.01);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVJCrashSizeUncertainty, "SVJ_CRASH_SIZE_UNCERTAINTY", 
                                "VolSVJ", "crashSizeUncertainty", 0.01);

//VolSVCJ
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJInitialVol, "SVCJ_INITIAL_VOL", 
                                "VolSVCJ", "initialVol", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJMeanVol, "SVCJ_MEAN_VOL", 
                                "VolSVCJ", "meanVol", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJMeanReversRate, "SVCJ_MEAN_REVERS_RATE", 
                                "VolSVCJ", "meanReversRate", 0.05);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJVolVol, "SVCJ_VOL_VOL", 
                                "VolSVCJ", "volVol", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJCorrelation, "SVCJ_CORRELATION", 
                                "VolSVCJ", "correlation", -0.01);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJCrashRate, "SVCJ_CRASH_RATE", 
                                "VolSVCJ", "commonCrashRate", 0.01);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJStockCrashSizeMean, "SVCJ_STOCK_CRASH_SIZE_MEAN", 
                                "VolSVCJ", "commonStockCrashSizeMean", -0.01);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJStockCrashSizeUncertainty, "SVCJ_STOCK_CRASH_SIZE_UNCERTAINTY", 
                                "VolSVCJ", "commonStockCrashSizeUncertainty", 0.01);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJVolCrashSizeMean, "SVCJ_VOL_CRASH_SIZE_MEAN", 
                                "VolSVCJ", "commonVolCrashSizeMean", 0.001);
GENERIC_PARAM_SENS_SINGLE_FIELD(SVCJStockVolCrashSizeCorrelation, "SVCJ_STOCK_VOL_CRASH_SIZE_CORRELATION", 
                                "VolSVCJ", "stockVolCrashSizeCorrelation", -0.01);


bool ParamSensLinkIn(){
    return (ParamSens::TYPE != 0);
}


DRLIB_END_NAMESPACE
