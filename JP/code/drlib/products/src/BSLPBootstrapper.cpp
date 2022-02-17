//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BSLPBootstrapper.cpp
//
//   Description : Class to bootstrap a SkewSurface for BSLP.
//
//   Author      : Matthias Arnsdorf
//
//   Date        : December 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/Results.hpp"
#include "edginc/BSLPBootstrapper.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// AUXILLIARY OBJECTIVE FUNCTION
///////////////////////////////////////////////////////////////////////////////
/** objective function for calibrating first point at each maturity */

class PRODUCTS_DLL AuxObjFunc :
    public Calibrator::ObjFuncLeastSquare
{
public:	
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~AuxObjFunc();

    /** Public constructor  */
    AuxObjFunc(
        OptimizerNDSP optimizer,
        TrancheIndexLeastSquareFitSP trancheObjFunc,
        CInstrumentArraySP quotes,
        Calibrator::InstanceIDArraySP ids
        );


    /**
    * Give a chance to do s'thing with the market
    * [Implements Calibrator::ObjFunc]
    * */
    virtual void getMarket(MarketData* market);

    /**
    * Returns an IObect that contains all the IAdjustable objects 
    * that the calibrator is to operate upon
    * [Implements Calibrator::ObjFunc]
    * */
    virtual IObjectSP getAdjustableGroup();

    /**
    * Returns the number of functions
    * [Implements Calibrator::ObjFuncLeastSquare]
    * */
    virtual int getNbFuncs() const;

    /**
    * Calculates the values of the objective functions
    * [Implements Calibrator::ObjFuncLeastSquare]
    * */
    virtual void calcValue(CDoubleArray&   funcvals) const;

    virtual double calcValue() const;


    // append objective function values
    void addObjFuncVals(DoubleArray & objFuncVals);

    // append calibrated function values
    void addCalibVals(DoubleArray & calibVals);


private:

    /** Private constructor (only build instances of that class using reflection) */
    AuxObjFunc();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Private constructor (only build instances of that class using reflection) */
    //AuxObjFunc();

    /** Default constructor */
    static IObject* defaultConstructor();

    // FIELDS
    /** Calibrator */
    CalibratorSP calibrator;

    /** objective function */
    TrancheIndexLeastSquareFitSP trancheObjFunc;

    /** aray of corresponding instruments (tranche quotes) */
    CInstrumentArraySP quotes;

    /** ids to bootstrap */
    Calibrator::InstanceIDArraySP ids;


    /** output array of calibrate values */
    mutable DoubleArraySP calibratedValues;

    /** output array of objective functions */
    mutable DoubleArraySP objFuncValues;


};
/** Destructor */
AuxObjFunc::~AuxObjFunc() {}

/** Invoked when Class is 'loaded' */
void AuxObjFunc::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(AuxObjFunc, clazz);
    SUPERCLASS(Calibrator::ObjFuncLeastSquare);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(calibrator);
    FIELD_NO_DESC(trancheObjFunc);
    FIELD_NO_DESC(quotes);
    FIELD_NO_DESC(ids);
    FIELD_NO_DESC(calibratedValues);
    FIELD_NO_DESC(objFuncValues);
}

/** Private constructor (only build instances of that class using reflection) */
AuxObjFunc::AuxObjFunc():
Calibrator::ObjFuncLeastSquare(TYPE),
calibrator(),
trancheObjFunc(),
quotes(),
ids()
{}

/** Default constructor */
IObject* AuxObjFunc::defaultConstructor() {
    return new AuxObjFunc();
}

/** TYPE for AuxObjFunc */
CClassConstSP const AuxObjFunc::TYPE =
CClass::registerClassLoadMethod(
                                "AuxObjFunc",
                                typeid(AuxObjFunc),
                                AuxObjFunc::load);



/** Public constructor  */
AuxObjFunc::AuxObjFunc(
    OptimizerNDSP optimizer,
    TrancheIndexLeastSquareFitSP trancheObjFunc,
    CInstrumentArraySP quotes,
    Calibrator::InstanceIDArraySP ids
    ):
Calibrator::ObjFuncLeastSquare(TYPE),
calibrator(new Calibrator(optimizer)),
trancheObjFunc(trancheObjFunc),
quotes(quotes),
ids(ids)
{
    // check that have at least two ids
    QLIB_VERIFY(ids->size()>1,"Need at least 2 id's");

    //ids and quotes have to have same length
    QLIB_VERIFY(ids->size()==quotes->size(), "Quotes and ids not same size");

    // intialise objfuncs and calibvals
    objFuncValues = DoubleArraySP(new DoubleArray(ids->size()-1));
    calibratedValues = DoubleArraySP(new DoubleArray(ids->size()-1));

}
/**
* Give a chance to do s'thing with the market
* [Implements Calibrator::ObjFunc]
* */
void AuxObjFunc::getMarket(MarketData* market)
{
    // not needed for AuxObjFunc
    // getMarket on input objFunc (TrancheIndexLeastSquareFit) should already have been called
    // at initialisation
}

/**
* Returns an IObect that contains all the IAdjustable objects 
* that the calibrator is to operate upon
* [Implements Calibrator::ObjFunc]
* */
IObjectSP AuxObjFunc::getAdjustableGroup()
{
    return trancheObjFunc->getAdjustableGroup();
}

/**
* Returns the number of functions
* [Implements Calibrator::ObjFuncLeastSquare]
* */
int AuxObjFunc::getNbFuncs() const
{
    // just calculating one function value for each choice of first contagion point
    return 1;
}

// append objective function values
void AuxObjFunc::addObjFuncVals(DoubleArray & objFuncVals)
{
    int i;
    for(i=0 ; i< objFuncValues->size(); i++)
    {
        objFuncVals.push_back((*objFuncValues)[i]);
    }

}

// append calibrated function values
void AuxObjFunc::addCalibVals(DoubleArray & calibVals)
{
    int i;
    for(i=0 ; i< calibratedValues->size(); i++)
    {
        calibVals.push_back((*calibratedValues)[i]);
    }

}



/** objective function value */
double AuxObjFunc::calcValue() const
{
    DoubleArray val(1);
    calcValue(val);
    return val[0];
}

/**
* Calculates the values of the objective functions
* [Implements Calibrator::ObjFuncLeastSquare]
* Here this is the objective function obtained by bootstrapping contagion facs
* contingent on first point
* */
void AuxObjFunc::calcValue(CDoubleArray& funcvals) const
{
    static const string method = "AuxObjFunc::calcValue";
    try
    {
        /** output value */
        double outVal = 0;
        // loop over points to bootstrap, i.e. points after first
        int i;
        for(i = 1; i< ids->size(); i++)
        {
            // set corresponding instrument in TrancheIndexLeastSquareFit
            // this is point with high strike corresponding to contagion fac point
            CInstrumentArraySP inst(new CInstrumentArray(1));
            (*inst)[0] = (*quotes)[i-1];
            trancheObjFunc->setCurrentInstruments(inst);

            Calibrator::InstanceIDArray theseIds(1);
            theseIds[0] = (*ids)[i];
            // Run Calibrator  to bootstrap this point
            CResultsSP currentResults =
                calibrator->run(*trancheObjFunc, theseIds);


            // calculate obj func and export obj func value
            // TODO this value is already available in currentResults
            // so should'nt need to recompute...
            (*objFuncValues)[i-1] = trancheObjFunc->calcValue(); 
            outVal+= (*objFuncValues)[i-1];


            // get calibrated value. //TODO check if can use getVAlue
            DoubleArray vals(0);
            Calibrator::getResultValues(currentResults, vals);

            (*calibratedValues)[i-1] = vals[0];

            // apply calibrated values to obj func 
            Calibrator::InstanceID::applyAdjustment(
                theseIds,
                trancheObjFunc->getAdjustableGroup(),
                vals);

        }

        // calculate value of final quote. Note that this point is not bootstrapped her
        // It is calculated in the outer loop which has this as obj function
        CInstrumentArraySP inst(new CInstrumentArray(1));
        (*inst)[0] = quotes->back();
        trancheObjFunc->setCurrentInstruments(inst);

        outVal+= trancheObjFunc->calcValue();

        // set output. Normalise by number of functions. Note that this is already a sum of squares
        funcvals[0] =  outVal/((double)ids->size());
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////
// BSLPBootstrapper
///////////////////////////////////////////////////////////////////////////////////////////

/** botstrap type corresponding to this algorithm */
const string BSLPBootstrapper::BOOTSTRAP_BSLP = "BSLP";


/** Invoked when Class is 'loaded' */
void BSLPBootstrapper::load(CClassSP& clazz) {
    clazz->setPrivate(); // don't make visible to EAS/spreadsheet
    REGISTER(BSLPBootstrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IBootstrapper);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(bootstrapMaturities);
    FIELD_NO_DESC(timePoints);
    FIELD_NO_DESC(currBootstrapIdx);
    FIELD_NO_DESC(objFunc);
    FIELD_NO_DESC(pointsToCalibrate);
    FIELD_NO_DESC(cdoQuotesIterator);
    FIELD_NO_DESC(name);
    FIELD_NO_DESC(fieldToCalibrate);
    FIELD_NO_DESC(rangeOverride);

}

/** Default constructor */
IObject* BSLPBootstrapper::defaultConstructor() {
    return new BSLPBootstrapper();
}


/** TYPE (for reflection) */
CClassConstSP const BSLPBootstrapper::TYPE =
CClass::registerClassLoadMethod(
                                "BSLPBootstrapper",
                                typeid(BSLPBootstrapper),
                                load);

/** Virtual destructor */
BSLPBootstrapper::~BSLPBootstrapper() {}

/** Constructor (internal - used by reflection) */
BSLPBootstrapper::BSLPBootstrapper():
CObject(TYPE),
bootstrapMaturities(0),
timePoints(0),
currBootstrapIdx(0),
objFunc(0),
pointsToCalibrate(0),
cdoQuotesIterator(0),
name(""),
fieldToCalibrate(""),
rangeOverride(0)
{}

/** Constructor (external) */
BSLPBootstrapper::BSLPBootstrapper(
    TrancheIndexLeastSquareFitSP objFunc,
    const Calibrator::InstanceIDArray & ids):
CObject(TYPE),
bootstrapMaturities(0),
timePoints(0),
currBootstrapIdx(0),
objFunc(objFunc),
pointsToCalibrate(0),
cdoQuotesIterator(0),
name(""),
fieldToCalibrate(""),
rangeOverride(0)
{
    static const string method = "BSLPBootstrapper::~BSLPBootstrapper";
    int i;

    QLIB_VERIFY(ids.size()>=1,"need at least one InstanceID");

    // assume that all ids correspond to skewSurface. This can be generalised
    SkewSurfaceConstSP  skewSurface = SkewSurfaceConstSP::dynamicCast(
        ids[0]->getObject(objFunc->getAdjustableGroup()));

    // assume that only have one field to calibrate
    fieldToCalibrate = ids[0]->getFieldName();

 

    // init name
    name = skewSurface->getName();

    // init range Override
    // assume for  now that ranges of all instanceID's are the same
    RangeArray rangeArray(ids[0]->numVariables());
    ids[0]->getRanges(rangeArray, 0);

    // assume for now that all ranges are the same
    rangeOverride = rangeArray[0];

    // construct bootstrap maturities
    bootstrapMaturities = DateTimeArraySP(new DateTimeArray(*skewSurface->maturities));
    DateTime::removeDuplicates(*bootstrapMaturities, true);

    //construct pointsToCalibrate
    // this is two dimensional in case we want to support calibrating several skew points
    // at each bootstrap step, e.g. the BOOTSTRAP_TIME case

    pointsToCalibrate.resize(0);

    for(i=0; i < bootstrapMaturities->size(); i++)
    {
        pointsToCalibrate.push_back(IntArraySP(new IntArray()));
    }
    for(i=0; i < (int)skewSurface->maturities->size(); i++)
    {
        DateTime date = (*skewSurface->maturities)[i];
        pointsToCalibrate[date.find(*bootstrapMaturities)]->push_back(i);
    }

    // WE have assumed that all points skews are to be calibrated
    // SHOULD use INstanceID idx for filter here instead // TODO


    // set up cdoQuoteBootstrappers
    // Ideally we should extract quote directly given skewSurface points
    // i.e. internalise mapping
    cdoQuotesIterator = objFunc->getQuotesIterator();

}



/** Main method that runs bootstrap calibration and returns results */
void BSLPBootstrapper::bootstrap(
                OptimizerNDSP optimizer,
                int & nbVars,                                 // (O) number of variables calibrated
                Calibrator::InstanceIDArray & aggregIds,    // (0)
                DoubleArray & aggregCalibValues,                   // (O) calibrated values
                DoubleArray & objFuncValues                  // (O) objective function values
                )
{
    static const string method = "BSLPBootstrapper::bootstrap";
    try
    {
        nbVars = 0;
        aggregIds.resize(0);
        aggregCalibValues.resize(0);
        objFuncValues.resize(0);

        // Bootstrapping loop in time
        for(init(); !end(); next())
        {
            // get instance ids for this time step
            Calibrator::InstanceIDArraySP theseIds = getCurrentInstanceIDs();

            //set up auxilliary objective function
            AuxObjFunc  auxObjFunc(
                optimizer,
                objFunc,
                cdoQuotesIterator->buildCurrentInstruments(),
                theseIds);

            // define calibrator
            Calibrator calibrator(optimizer);
            // Run Calibrator to calibrate first point
            Calibrator::InstanceIDArraySP firstPoint(new Calibrator::InstanceIDArray(1));
            (*firstPoint)[0] = theseIds->front();
            
            
            CResultsSP currentResults =
                calibrator.run(auxObjFunc, *firstPoint);

            // calculate obj func for first point and export obj func value
            // TODO this value is already available in currentResults
            // so should'nt need to recompute...
            objFuncValues.push_back(auxObjFunc.calcValue());
            // add objective function vals for other points at this maturity
            auxObjFunc.addObjFuncVals(objFuncValues);

            // get calibrated value for first point
            DoubleArray vals(0);
            Calibrator::getResultValues(currentResults, vals);

            // apply calibrated values to obj func
            Calibrator::InstanceID::applyAdjustment(
                *theseIds,
                objFunc->getAdjustableGroup(),
                vals);
   

            // append value of first point to aggregCalibValues
            int i;
            for(i = 0; i < vals.size(); i++)
                aggregCalibValues.push_back(vals[i]);

            // append remaining values at this maturity
            auxObjFunc.addCalibVals(aggregCalibValues);

            // add ids for values we are calibrating and ammend numVariables
            for (int iIds = 0 ; iIds < theseIds->size() ; iIds++){
                aggregIds.push_back((*theseIds)[iIds]);
                nbVars += (*theseIds)[iIds]->numVariables();
            }

        }

    }
    catch( exception & e)
    {
        throw ModelException(e, method);
    }


}

/**
* Method called before first step of the loop
* [Implements IBootstrapperTime]
* */
void BSLPBootstrapper::init() 
{
    currBootstrapIdx = 0;
    cdoQuotesIterator->init();
}

/**
* Method called after each step of the loop
* [Implements IBootstrapperTime]
* */
void BSLPBootstrapper::next() 
{
    // check that SkewSurface iteration has not ended already
    QLIB_VERIFY(!end(), "Next called past end of skew surface iteration");
    // iterate skewSurface
    currBootstrapIdx++;

    // check that cdoQuoteIterator has not ended already
    QLIB_VERIFY(!cdoQuotesIterator->end(), 
        "Not enough quotes to bootstrap. CdoQuotesIterator ends to early");
    // iterate quotes   
    cdoQuotesIterator->next();     
}

/**
* Method called to test the end of the loop
* [Implements IBootstrapperTime]
* */
bool BSLPBootstrapper::end() const 
{
    bool isEnd = (currBootstrapIdx >= pointsToCalibrate.size());
   
    return isEnd;
}

/**
* Returns "state" corresponding to current step of the loop
* In time only case this just corresponds to bootstrap date
* */
IObjectSP BSLPBootstrapper::getCurrentState() const 
{
    return DateTimeSP(new DateTime((*bootstrapMaturities)[currBootstrapIdx]));  
}   

/**
* Produces an InstanceID corresponding to the current bootstrapping step
* [Implements IInstanceIDBootstrapperTime]
* */
Calibrator::InstanceIDArraySP BSLPBootstrapper::getCurrentInstanceIDs() const 
{

    Calibrator::InstanceIDArraySP outArray(new Calibrator::InstanceIDArray(0));

    int numPoints = pointsToCalibrate[currBootstrapIdx]->size();
    int t;
    for(t=0; t < numPoints; t++)
    {
        outArray->push_back(Calibrator::InstanceIDDbSP(new Calibrator::InstanceIDDb(
            SkewSurface::TYPE->getName(),
            name, 
            fieldToCalibrate,             // field name
            false,                                   // do not override
            0.0,                                     // dummy override value
            (*pointsToCalibrate[currBootstrapIdx])[t],  // index
            rangeOverride							 // override for range	
            ))); 

    }


    return outArray;
}




DRLIB_END_NAMESPACE
