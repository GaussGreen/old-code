//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : TrancheIndexLeastSquareFit.cpp
//
//   Description : Objective function used to calibrate CDO models
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TrancheIndexLeastSquareFit.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/SkewSurfaceBootstrapper.hpp"
#include "edginc/BSLPBootstrapper.hpp"
#include "edginc/FDModel.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/Timer.hpp"
#include "edginc/CDO.hpp"
#include "edginc/RFLMixtureModelCalibrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** flag to indicate whether we want to output the number of calls to calcValue */
//#define COUNT_VALUE_CALLS
#ifdef COUNT_VALUE_CALLS
	static int COUNTER = 0;
#endif

// MTM^2
const string TrancheIndexLeastSquareFit::OBJ_FUNC_MTM = "MTM";
// (MTM / FeeLeg)^2
const string TrancheIndexLeastSquareFit::OBJ_FUNC_MTM_FEE_LEG = "MTM/FEE_LEG";

const double TrancheIndexLeastSquareFit::DEFAULT_OBJ_FUNC_FLOOR = 0.0;

/**
 * Give a chance to do s'thing with the market
 * [Implements Calibrator::ObjFunc]
 * */
void TrancheIndexLeastSquareFit::getMarket(MarketData* market)
{
    cdoQuotes.getData(model.get(), market);
    CDOQuotesConstSP processedQuotes = cdoQuotes->buildCDOQuotes();

    // Use fine grid (if any)
    CDOQuotesConstSP extendedQuotes;
    if (fineGrid.get() != 0) 
    {
		  fineGrid->getMarket(market);
      extendedQuotes = fineGrid->extendQuotes(processedQuotes);
    } 
    else 
    {
      extendedQuotes = processedQuotes;
    }
    
    // Init instruments & controls
    CInstrumentArraySP inst;
    instruments = CInstrumentArraySP(new CInstrumentArray());
    controls = CControlArraySP(new CControlArray());
    
	if(bootstrapType == SkewSurfaceBootstrapper::BOOTSTRAP_TIME_STRIKE ||
        bootstrapType == BSLPBootstrapper::BOOTSTRAP_BSLP)
	{
		cdoQuotesIterator =
			CDOQuotesBootstrapperSP(new CDOQuotesBootstrapper(extendedQuotes));
	}
	else if(bootstrapType == SkewSurfaceBootstrapper::BOOTSTRAP_TIME)
	{
		cdoQuotesIterator =
			CDOQuotesBootstrapperTimeOnlySP(new CDOQuotesBootstrapperTimeOnly(extendedQuotes));
	}
	else
	{
		throw ModelException("bootstrapType ("+ bootstrapType +") is not supported");
	}

	
    OutputRequestArraySP outputRequest(new OutputRequestArray(0));
    if(objFuncType == OBJ_FUNC_MTM_FEE_LEG)
    {
        outputRequest->push_back(OutputRequestSP(new OutputRequest(OutputRequest::TRANCHE_FEE_LEG_PRICE)));
    }
    
    CControlSP control = CControlSP(new Control(
        SensitivityArraySP(new SensitivityArray(0)),
        outputRequest,
        false,
        ""));

	for (cdoQuotesIterator->init(); !cdoQuotesIterator->end(); cdoQuotesIterator->next()) 
  {

        // get instruments to value at this bootstrap state
        inst = cdoQuotesIterator->buildCurrentInstruments();
        
        // loop through instruments for this bootstrap state
		int i = 0;
		for(; i < (int)inst->size(); i++)
		{

            CDO* cdoInst = dynamic_cast<CDO*>((*inst)[i].get());
            CDOPortfolioSP ptf = cdoInst->getCDOPortfolio();
            if (initialGuess.get())
            ptf->setEngineParameters(initialGuess);

        // disabled here for now since it causes BSLP failure (MA)
        currentGuess = (dynamic_cast<const PortfolioName&>(*(ptf->getInnerLossConfig(0)))).getCreditAsset()->getEngineParams();
            
            Calibrator::ObjFunc::Helper::getMarket(
				(*inst)[i].get(), 
				model.get(),
				control.get(),
				market);
			
			// in boostrapping case instruments field should be reset at each
            // bootstrap iteration by IBootstrapper
			instruments->push_back((*inst)[i]);
				
			controls->push_back(control);
		}
  }
}

/**
 * Returns an IObect that contains all the IAdjustable objects 
 * that the calibrator is to operate upon
 * [Implements Calibrator::ObjFunc]
 * */
IObjectSP TrancheIndexLeastSquareFit::getAdjustableGroup() {

	// temporary: Have special case if model is of type SpreadLossTree
	// want to search for field to calibrate in Model and not in instruments

	if(FDModel::TYPE->isInstance(model.get()))
	{
		return model;
	}
	else 
	{
        return instruments;
    }
}

/**
 * Returns the number of functions
 * [Implements Calibrator::ObjFuncLeastSquare]
 * */
int TrancheIndexLeastSquareFit::getNbFuncs() const 
{
     return instruments->size();
}

/**
 * Calculates the values of the objective functions
 * [Implements Calibrator::ObjFuncLeastSquare]
 * */

void TrancheIndexLeastSquareFit::calcValue(CDoubleArray& funcvals) const {
    const string method = "TrancheIndexLeastSquareFit::calcValue(...)";
    try 
    {
        CResultsSP results;
        
		// sanity check
		if (funcvals.size() != instruments->size()) 
		{
			throw ModelException(method,
				"Internal error: inconsistent number of least square values");
		}
        

		// loop over instruments
		for (int i=0; i< (int)instruments->size(); i++) 
		{
            double weight = 1.0;
            if (cdoQuotesWeights.get())
            {
                double      k1,k2;
                DateTime    t;
                // casts the instrument into a CDO
                CDO* cdoInst = dynamic_cast<CDO*>((*instruments)[i].get());

                cdoInst->getMaturityAndStrikesPercent(t, k1, k2);
                weight = (*cdoQuotesWeights.get()).getWeight(t,k1,k2);
            }

#ifdef COUNT_VALUE_CALLS
			COUNTER++; 
			ErrorHandler::writeMsg("Num value calls: " +Format::toString(COUNTER) );
			cout << "Num value calls: " << COUNTER << endl;

			Timer timer;
#endif

            if (weight != 0.0)
            {
    			results = CResultsSP(model->Run(
				    (*instruments)[i].get(),
				    (*controls)[i].get()));
    			
                if(objFuncType == OBJ_FUNC_MTM)
                {
                    funcvals[i] = results->retrievePrice();
                }
                else if (objFuncType == OBJ_FUNC_MTM_FEE_LEG)
                {
                    const CDOLegOutput& out = 
                        dynamic_cast<const CDOLegOutput&>(*(results->retrieveRequestResult("TRANCHE_FEE_LEG_PRICE")));
                    double feeLegVal = out.Price();         
/*                        CDoubleConstSP::dynamicCast(
                        results->retrieveRequestResult("TRANCHE_FEE_LEG_PRICE")
                        )->doubleValue();
*/                    
                    funcvals[i] = results->retrievePrice() / feeLegVal;
                }
                else
                {
                    throw ModelException("objFuncType ("+objFuncType+") is not supported.");
                }

                // multiply by the weight computed above
                funcvals[i] *= weight;
            }
            else
            {
                // the instrument is not used at all in that case
                funcvals[i] = 0.0;
            }

#ifdef COUNT_VALUE_CALLS
			double valTime = timer.calcTime();
			ErrorHandler::writeMsg("Price time: " +Format::toString(valTime) );
			cout << "Price time: " << valTime << endl;
#endif
		} // end instruments loop
		
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** override for calcValue() so that can put a floor on the objective function */
double TrancheIndexLeastSquareFit::calcValue() const{
	try{
		int nbFuncs = getNbFuncs();

		if (!funcs){
			funcs = CDoubleArraySP(new CDoubleArray(nbFuncs));
		} else {
			if (nbFuncs != funcs->size()) {
				// nbFuncs may vary for the same objective function
				funcs = CDoubleArraySP(new CDoubleArray(nbFuncs));
			}
		}
		DoubleArray& funcvals = *funcs;
		calcValue(funcvals);
		double res = 0.0;
		int iFunc = 0;
		for (; iFunc < nbFuncs; ++iFunc)
			res += Maths::square(funcvals[iFunc]);
    res += currentGuess->ParameterBadness();
		if (CString::equalsIgnoreCase(objFuncNormalized,TrancheIndexLeastSquareFit::DO_NOT_NORMALIZE)){
			return res;
		}
		else{
			return Maths::max((res / nbFuncs), objFuncFloor);
		}
	}
	catch(exception& e){
		throw ModelException(e, "Calibrator::ObjFuncLeastSquare::calcValue");
	}    
}


/** main method used to retrieve a bootstrapper capable of bootstrapping */
IBootstrapperSP TrancheIndexLeastSquareFit::getBootstrapper(
    const Calibrator::InstanceIDArray & ids)
{

    IBootstrapperSP out;
    // need to check bootstrap type and return a bootstrapper here
    if(bootstrapType == SkewSurfaceBootstrapper::BOOTSTRAP_TIME_STRIKE ||
        bootstrapType == SkewSurfaceBootstrapper::BOOTSTRAP_TIME)
    {
        
        return SkewSurfaceBootstrapperSP(new SkewSurfaceBootstrapper(
            TrancheIndexLeastSquareFitSP(this), ids));
            

    }
    else if (bootstrapType == "RFLMixtureCalibration")
    {
      return RFLMixtureModelCalibratorSP(new RFLMixtureModelCalibrator(
            TrancheIndexLeastSquareFitSP(this), ids));
    }
    else if(bootstrapType == BSLPBootstrapper::BOOTSTRAP_BSLP)
    {
        return BSLPBootstrapperSP(new BSLPBootstrapper(
            TrancheIndexLeastSquareFitSP(this), ids));
    }
    else
    {
        throw ModelException("bootstrapType ("+ bootstrapType +") is not supported");
    }

    return out;
}



/** set the current instruments valued in calcValue. Used fore bootstrapping */
void TrancheIndexLeastSquareFit::setCurrentInstruments(CInstrumentArraySP inst)
{
    instruments = inst;
}

/** Destructor */
TrancheIndexLeastSquareFit::~TrancheIndexLeastSquareFit() {}

/** Invoked when Class is 'loaded' */
void TrancheIndexLeastSquareFit::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(TrancheIndexLeastSquareFit, clazz);
    SUPERCLASS(Calibrator::ObjFuncLeastSquare);
    IMPLEMENTS(Calibrator::ObjFunc::IGenericBootstrappable);
    EMPTY_SHELL_METHOD(defaultConstructor);
	FIELD(model, "CDO model");

	FIELD(cdoQuotes, "Index tranche quotes or generator of quotes");

    FIELD(cdoQuotesWeights, "Index tranche Weights");
    FIELD_MAKE_OPTIONAL(cdoQuotesWeights);

    FIELD(initialGuess, "initial guess override for the calibration");
    FIELD_MAKE_OPTIONAL(initialGuess);

    FIELD(bootstrapType, 
		SkewSurfaceBootstrapper::BOOTSTRAP_TIME_STRIKE + ": bootstrap both in time and strike dim (default). " +
		SkewSurfaceBootstrapper::BOOTSTRAP_TIME + ": bootstrap in time dim only.");
	FIELD_MAKE_OPTIONAL(bootstrapType);
	FIELD(objFuncFloor, "Floor for objective function. [default = "+
		Format::toString(DEFAULT_OBJ_FUNC_FLOOR)+"]");

    FIELD(objFuncType,
        OBJ_FUNC_MTM + ": (MTM)^2 (default). " +
        OBJ_FUNC_MTM_FEE_LEG +": (MTM / Fee_Leg)^2.");
    FIELD_MAKE_OPTIONAL(objFuncType);
	FIELD_MAKE_OPTIONAL(objFuncFloor);
	FIELD(fineGrid, "Defines the 'fine grid'");
	FIELD_MAKE_OPTIONAL(fineGrid);

	FIELD_NO_DESC(cdoQuotesIterator);
	FIELD_MAKE_TRANSIENT(cdoQuotesIterator);
    FIELD_NO_DESC(instruments);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(instruments);
    FIELD_NO_DESC(idBootstrappers);
    FIELD_MAKE_TRANSIENT(idBootstrappers);
    FIELD_NO_DESC(controls);
    FIELD_MAKE_TRANSIENT(controls);
    FIELD_NO_DESC(doBootstrapping);
    FIELD_MAKE_TRANSIENT(doBootstrapping);
	FIELD_NO_DESC(funcs);
	FIELD_MAKE_TRANSIENT(funcs);
    
}

/** Private constructor (only build instances of that class using reflection) */
TrancheIndexLeastSquareFit::TrancheIndexLeastSquareFit() :
	Calibrator::ObjFuncLeastSquare(TYPE),
    cdoQuotes(0),
    cdoQuotesWeights(0),
    initialGuess(0),
	cdoQuotesIterator(0),
    instruments(0),
    idBootstrappers(0),
    controls(0),
    doBootstrapping(false),
	bootstrapType(SkewSurfaceBootstrapper::BOOTSTRAP_TIME_STRIKE),
    objFuncType(OBJ_FUNC_MTM),
	objFuncFloor(DEFAULT_OBJ_FUNC_FLOOR)
	{}

/** Default constructor */
IObject* TrancheIndexLeastSquareFit::defaultConstructor() {
    return new TrancheIndexLeastSquareFit();
}

/** TYPE for TrancheIndexLeastSquareFit */
CClassConstSP const TrancheIndexLeastSquareFit::TYPE =
	CClass::registerClassLoadMethod(
    	"TrancheIndexLeastSquareFit",
    	typeid(TrancheIndexLeastSquareFit),
    	TrancheIndexLeastSquareFit::load);
    	
// for class loading 
bool TrancheIndexLeastSquareFitLoad() {
    return (TrancheIndexLeastSquareFit::TYPE != 0);
}
  	

DRLIB_END_NAMESPACE

