//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : RFLMixtureModelCalibrator.cpp
//
//   Description : Class to calibrate the RFLMixtureModel
//
//   Author      : Jakob Sidenius
//
//   Date        : December 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/Results.hpp"
#include "edginc/RFLMixtureModelCalibrator.hpp"
#include "edginc/ConvolutionEngine.hpp"
#include "edginc/BespokeCDOModelTest.hpp"
#include "edginc/ConvolutionModelConfig.hpp"

DRLIB_BEGIN_NAMESPACE

/** TYPE (for reflection) */
CClassConstSP const RFLMixtureModelCalibrator::TYPE =
CClass::registerClassLoadMethod(
                                "RFLMixtureModelCalibrator",
                                typeid(RFLMixtureModelCalibrator),
                                load);

/** Constructor (internal - used by reflection) */
RFLMixtureModelCalibrator::RFLMixtureModelCalibrator():
CObject(TYPE),
objFunc(0),
theIds(0)
{
  // empty
}

/** Constructor (external) */
RFLMixtureModelCalibrator::RFLMixtureModelCalibrator(
    TrancheIndexLeastSquareFitSP objFunc,
    const Calibrator::InstanceIDArray & ids):
CObject(TYPE),
objFunc(objFunc),
theIds(new Calibrator::InstanceIDArray(0))
{
  theIds->reserve(ids.size());
  for (size_t i = 0; i < ids.size(); ++i)
    theIds->push_back(ids[i]);
}

struct RootSearchHelper
{
  bool xIsW;
  TrancheIndexLeastSquareFitSP  objFunc;
}; // RootSearchHelper

double FunctionToRootSearch(double x, void  *s)
{
  RootSearchHelper* data = static_cast<RootSearchHelper*>(s);
  TrancheIndexLeastSquareFit& objFunc = *(data->objFunc);
  if (data->xIsW)
  {
    // Set parameter
    ConvolutionEngine& ce = dynamic_cast<ConvolutionEngine&>(*objFunc.Model());
    IConvolutionModel& cm = *ce.ConvolutionModel();
    BespokeCDOModelTest& bcmt = dynamic_cast<BespokeCDOModelTest&>(cm);
    IModelConfigMapper& mcm = *bcmt.ModelConfigMapper();
    ConvolutionModelConfig& cmc = dynamic_cast<ConvolutionModelConfig&>(mcm);
    IConditionalDefaultsModel& cdm = *cmc.CondDefaultsModel();
    RFLMixtureDefaultsModel& model = dynamic_cast<RFLMixtureDefaultsModel&>(cdm);

    model.SetW(x);

    // Set instrument and controls for most senior tranche
    const CInstrumentArraySP instruments = objFunc.GetInstruments();
    const size_t nInsts = instruments->size(); 
    const CControlArraySP controls = objFunc.GetControls(); 

    // Calculate function value
  	CResultsSP results = CResultsSP(objFunc.Model()->Run((*instruments)[nInsts-1].get(),(*controls)[nInsts-1].get()));
    return results->retrievePrice();
  }
  // Set parameter
  const_cast<RflMixtureParameters&>(dynamic_cast<const RflMixtureParameters&>(*(objFunc.CurrentGuess()))).SetGlobalBeta(x);

  // Set instrument and control for equity junior tranche
  const CInstrumentArraySP instruments = objFunc.GetInstruments();
  const CControlArraySP controls = objFunc.GetControls(); 

  // Calculate function value
	CResultsSP results = CResultsSP(objFunc.Model()->Run((*instruments)[0].get(),(*controls)[0].get()));
  return results->retrievePrice();
}

/** Main method that runs bootstrap calibration and returns results */
void RFLMixtureModelCalibrator::bootstrap(
    OptimizerNDSP optimizer,
    int & nbVars,                                 // (O) number of variables calibrated
    Calibrator::InstanceIDArray & aggregIds,    // (0)
    DoubleArray & aggregCalibValues,                   // (O) calibrated values
    DoubleArray & objFuncValues)                  // (O) objective function values
{
    static const string functionName = "RFLMixtureModelCalibrator::bootstrap";
    try
    {
      // Preprocessing step
      // Here we alternate between root searches which match the most junior and the most senior input tranches, resp.
      // These are assumed to be first and last in the input instrument array. We stop when both are matched to sufficient
      // accuracy. 

      RootSearchHelper helper;
      helper.xIsW = true;
      helper.objFunc = objFunc;

      const double tol = 1.0e-6;
      double pvSenior = 1.0, pvEquity = 1.0;
      for(; pvSenior*pvSenior + pvEquity*pvEquity < tol; helper.xIsW = !helper.xIsW)
      {
        double rangeLow = 0.0, rangeHigh = 1.0, funcLow, funcHigh;
        ZbracReturn zbrac = zbracUseful(&FunctionToRootSearch, &helper, &rangeLow, &rangeHigh, &funcLow, &funcHigh);
        
        // Call zbrent
        try 
        {
          double root = zbrentUsefulBoundary(&FunctionToRootSearch, &helper, rangeLow, rangeHigh, 0.1 * tol, funcLow, funcHigh);
        }
        catch (NRException&) 
        {
          // This is zbrent failing to solve - hide "zbrent" message,
          // which confuses users, and output a more meaningful error
          throw NRException(functionName, "Zbrent error: failed to find root");
        }
      }

      // Optimizing step
      // Here we let an optimizer loose on the problem of minimizing the pricing error on all given tranches
      CResultsSP currentResults = Calibrator(optimizer).run(*objFunc, *theIds);
      aggregCalibValues.clear();
      Calibrator::getResultValues(currentResults, aggregCalibValues);
      aggregIds = *theIds;
    }
    catch( exception & e)
    {
      throw ModelException(e, functionName);
    }
}


/** Invoked when Class is 'loaded' */
void RFLMixtureModelCalibrator::load(CClassSP& clazz) {
    clazz->setPrivate(); // don't make visible to EAS/spreadsheet
    REGISTER(RFLMixtureModelCalibrator, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IBootstrapper);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(objFunc);
    FIELD_NO_DESC(theIds);   
}

/** Default constructor */
IObject* RFLMixtureModelCalibrator::defaultConstructor() {
    return new RFLMixtureModelCalibrator();
}


DRLIB_END_NAMESPACE
