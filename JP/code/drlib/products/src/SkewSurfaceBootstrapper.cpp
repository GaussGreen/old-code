//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : SkewSurfaceBootstrapper.cpp
//
//   Description : Class to bootstrap a SkewSurface.
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
#include "edginc/SkewSurfaceBootstrapper.hpp"

DRLIB_BEGIN_NAMESPACE

/** bootstrap in time and strike dimension */
const string SkewSurfaceBootstrapper::BOOTSTRAP_TIME_STRIKE = "TIME_STRIKE";

/** bootstrap in time dimension only (global calibration in strike dim) */
const string SkewSurfaceBootstrapper::BOOTSTRAP_TIME = "TIME";

/** TYPE (for reflection) */
CClassConstSP const SkewSurfaceBootstrapper::TYPE =
CClass::registerClassLoadMethod(
                                "SkewSurfaceBootstrapper",
                                typeid(SkewSurfaceBootstrapper),
                                load);

/** Virtual destructor */
SkewSurfaceBootstrapper::~SkewSurfaceBootstrapper() {}

/** Constructor (internal - used by reflection) */
SkewSurfaceBootstrapper::SkewSurfaceBootstrapper():
CObject(TYPE),
bootstrapType(""),
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
SkewSurfaceBootstrapper::SkewSurfaceBootstrapper(
    TrancheIndexLeastSquareFitSP objFunc,
    const Calibrator::InstanceIDArray & ids):
CObject(TYPE),
bootstrapType(""),
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
    static const string method = "SkewSurfaceBootstrapper::~SkewSurfaceBootstrapper";
    int i;

    QLIB_VERIFY(ids.size()>=1,"need at least one InstanceID");

    // assume that all ids correspond to skewSurface. This can be generalised
    SkewSurfaceConstSP  skewSurface = SkewSurfaceConstSP::dynamicCast(
        ids[0]->getObject(objFunc->getAdjustableGroup()));

    // assume that only have one field to calibrate
    fieldToCalibrate = ids[0]->getFieldName();

    bootstrapType = objFunc->getBootstrapType();

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
    if(bootstrapType==BOOTSTRAP_TIME_STRIKE)
    {
        for (i = 0; i < skewSurface->strikes->size(); ++i) {
            // ignore "extrapolated" betas i.e. "0%" and "100%" betas
            // Should really use index in INstanceID for filter // TODO
            if ((*skewSurface->strikes)[i] > 0.0 && (*skewSurface->strikes)[i] < 1.0) 
            {
                pointsToCalibrate.push_back(IntArraySP(new IntArray(1)));
                (*pointsToCalibrate.back())[0] = i;
            }
        }   
    }
    else if(bootstrapType==BOOTSTRAP_TIME)
    {
        for(i=0; i < bootstrapMaturities->size(); i++)
        {
            pointsToCalibrate.push_back(IntArraySP(new IntArray()));
        }
        for(i=0; i < (int)skewSurface->maturities->size(); i++)
        {
            DateTime date = (*skewSurface->maturities)[i];
            pointsToCalibrate[date.find(*bootstrapMaturities)]->push_back(i);
        }

    }
    else
    {
        throw ModelException("bootstrapType = "+ bootstrapType +" is not supported", method);
    }
    
   
    // SHOULD use INstanceID idx for filter here // TODO


    if(bootstrapType == BOOTSTRAP_TIME_STRIKE)
    {
        // This is just needed for validation in TIME_STRIKE case
        // init timePoints:
        // We map the (maturity, strike) surface (excluding strikes <= 0) into
        // a (maturity, 0%, strike) cube. Eg:
        // 1Y,0% => ignored
        // 1Y,3% => 1Y,0%,3%
        // 1Y,6% => 1Y,3%,6%
        // 2Y,3% => 2Y,0%,3%
        // 2Y,6% => 2Y,3%,6%
        int nbPoints = pointsToCalibrate.size();
        timePoints = TimePoint2DArraySP(new TimePoint2DArray(nbPoints));
        for (i = 0; i < nbPoints; ++i)
        {
            if (i == 0 || (*skewSurface->strikes)[(*pointsToCalibrate[i])[0]] 
            < (*skewSurface->strikes)[(*pointsToCalibrate[i-1])[0]])

            {
                (*timePoints)[i] = TimePoint2DSP(new TimePoint2D(
                    (*skewSurface->maturities)[(*pointsToCalibrate[i])[0]],
                    0.0,
                    (*skewSurface->strikes)[(*pointsToCalibrate[i])[0]]));
            }
            else
            {
                (*timePoints)[i] = TimePoint2DSP(new TimePoint2D(
                    (*skewSurface->maturities)[(*pointsToCalibrate[i])[0]],
                    (*skewSurface->strikes)[(*pointsToCalibrate[i-1])[0]],
                    (*skewSurface->strikes)[(*pointsToCalibrate[i])[0]]));            
            }
        }

    }

    // set up cdoQuoteBootstrappers
    // Ideally we should extract quote directly given skewSurface points
   // i.e. internalise mapping
    cdoQuotesIterator = objFunc->getQuotesIterator();

}



/** Main method that runs bootstrap calibration and returns results */
void SkewSurfaceBootstrapper::bootstrap(
    OptimizerNDSP optimizer,
    int & nbVars,                                 // (O) number of variables calibrated
    Calibrator::InstanceIDArray & aggregIds,    // (0)
    DoubleArray & aggregCalibValues,                   // (O) calibrated values
    DoubleArray & objFuncValues                  // (O) objective function values
    )
{
    static const string method = "SkewSurfaceBootstrapper::bootstrap";
    try
    {
        nbVars = 0;
        aggregIds.resize(0);
        aggregCalibValues.resize(0);
        objFuncValues.resize(0);

        // Bootstrapping loop
        for(init(); !end(); next())
        {

            // check that bootstrappers are in the same state
            QLIB_VERIFY(
                getCurrentState()->equalTo(cdoQuotesIterator->getCurrentState().get()),
                "Bad mapping between tranche quotes and skewSurface. Iterators not in same state"
                );


            Calibrator::InstanceIDArray theseIds = *getCurrentInstanceIDs();


            for (int iIds = 0 ; iIds < theseIds.size() ; iIds++){
                aggregIds.push_back(theseIds[iIds]);
                nbVars += theseIds[iIds]->numVariables();
            }

            // need to set the current instruments to be valued in objective function
            objFunc->setCurrentInstruments(cdoQuotesIterator->buildCurrentInstruments());

            Calibrator calibrator(optimizer);
            // Run Calibrator !
            CResultsSP currentResults =
                calibrator.run(*objFunc, theseIds);

            // calculate obj func and export obj func value
            // TODO this value is already available in currentResults
            // so should'nt need to recompute...
            objFuncValues.push_back(objFunc->calcValue());

            DoubleArray vals(0);
            Calibrator::getResultValues(currentResults, vals);

            
            // apply calibrated values to obj func
            Calibrator::InstanceID::applyAdjustment(theseIds,
                objFunc->getAdjustableGroup(),
                vals);

            // append vals to aggregCalibValues
            int i;
            for(i = 0; i < vals.size(); i++)
                aggregCalibValues.push_back(vals[i]);

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
void SkewSurfaceBootstrapper::init() 
{
    currBootstrapIdx = 0;
    cdoQuotesIterator->init();
}

/**
* Method called after each step of the loop
* [Implements IBootstrapperTime]
* */
void SkewSurfaceBootstrapper::next() 
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
bool SkewSurfaceBootstrapper::end() const 
{
    bool isEnd = (currBootstrapIdx >= pointsToCalibrate.size());
    //check codquotesIterator and skewSurface have ended at same time
    // Don' think we need this but copied over from old implementation (MA)
    if(bootstrapType==BOOTSTRAP_TIME_STRIKE)
    {
        if(isEnd && !cdoQuotesIterator->end())
        {
            throw ModelException("cdoQuotes has different length to SkewSurface. Iterators end at different points",
                "SkewSurfaceBootstrapper::end()");
        }
    }
    return isEnd;
}

/**
* Returns "state" corresponding to current step of the loop
* In time only case this just corresponds to bootstrap date
* */
IObjectSP SkewSurfaceBootstrapper::getCurrentState() const 
{
    if(bootstrapType==BOOTSTRAP_TIME_STRIKE)
    {
        return (*timePoints)[currBootstrapIdx];
    }
    else
    {
        return DateTimeSP(new DateTime((*bootstrapMaturities)[currBootstrapIdx]));
    }
}   

/**
* Produces an InstanceID corresponding to the current bootstrapping step
* [Implements IInstanceIDBootstrapperTime]
* */
Calibrator::InstanceIDArraySP SkewSurfaceBootstrapper::getCurrentInstanceIDs() const 
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



/** Invoked when Class is 'loaded' */
void SkewSurfaceBootstrapper::load(CClassSP& clazz) {
    clazz->setPrivate(); // don't make visible to EAS/spreadsheet
    REGISTER(SkewSurfaceBootstrapper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IBootstrapper);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(bootstrapType);
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
IObject* SkewSurfaceBootstrapper::defaultConstructor() {
    return new SkewSurfaceBootstrapper();
}


DRLIB_END_NAMESPACE
