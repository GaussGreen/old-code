//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VAsset.hpp
//
//   Description : Implement volatility asset such as VIX
//
//
//   Author      : xiaolan zhang
//
//   Date        : 24 Mar 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VAsset.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/AssetUtil.hpp"
//#include "edginc/MarketObservable.hpp"


DRLIB_BEGIN_NAMESPACE

const string VIXFModel::BS_SV = "bs_sv";
const string VIXFModel::SV = "sv";
const string VIXFModel::BS_SV_NO_BASIS = "bs_sv_no_basis";
const string VIXFModel::SV_NO_BASIS = "sv_no_basis";
const string VIXFModel::VSW_BS= "vsw_bs";
const string VIXFModel::VSW_SV= "vsw_sv";

void VIXFModel::validate(){

    //Validate pricingMethod used to priceing VIX Future/ Forward
    if (   !CString::equalsIgnoreCase(pricingMethod, BS_SV)
        && !CString::equalsIgnoreCase(pricingMethod, SV)
        && !CString::equalsIgnoreCase(pricingMethod, BS_SV_NO_BASIS)
        && !CString::equalsIgnoreCase(pricingMethod, SV_NO_BASIS)
        && !CString::equalsIgnoreCase(pricingMethod, VSW_BS)
        && !CString::equalsIgnoreCase(pricingMethod, VSW_SV) ){
        throw ModelException("pricingMethod chosen is wrong");
    }
}

// inherited from VanVSModel
void VIXFModel::Price(CInstrument*  instrument,
                       CControl*     control,
                       CResults*     results)
{
    static const string method = "VIXFModel::Price";
    IProduct* product = 0;
    try {
        validate();

        if (!IIntoProduct::TYPE->isInstance(instrument)){
            throw ModelException("Instrument of type "+
                                 instrument->getClass()->getName() +
                                 " does not support VIXFModel::IntoProduct");
        }
        if (instrument->priceDeadInstrument(control, results)) {
            return; // done for a dead instrument
        }
        // cast to VIXFModel::IIntoProduct
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);

        // create and price the product
        product = intoProd.createProduct(this);
        product->price(this, control, results);
        delete product;
    }
    catch (exception& e) {
        delete product;
        throw ModelException(e, method);
    }
};

CClassConstSP const VIXFModel::TYPE = 
    CClass::registerClassLoadMethod("VIXFModel", typeid(VIXFModel),
    VIXFModel::load);

CClassConstSP const VIXFModel::IIntoProduct::TYPE =
    CClass::registerInterfaceLoadMethod("VIXFModel::IIntoProduct",
                                        typeid(VIXFModel::IIntoProduct), 
                                        VIXFModel::loadIntoProduct);

//////////////////////////////////////////////////////////////////////////////////////////////

// IMS compatible VIXFModel shell for ClosedFormIntegrateLN (varswap)
// and FourierEngine (convexity)  (SVJ)
class BSVFutSVJ: public VIXFModel {
public:
    static CClassConstSP const TYPE;

private:
    BSVFutSVJ():VIXFModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(BSVFutSVJ, clazz);
        SUPERCLASS(VIXFModel);
        EMPTY_SHELL_METHOD(defaultBSVFutSVJ);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultBSVFutSVJ(){
        return new BSVFutSVJ();
    }
};

CClassConstSP const BSVFutSVJ::TYPE =
CClass::registerClassLoadMethod(
    "BSVFutSVJ", typeid(BSVFutSVJ), load);


//////////////////////////////////////////////////////////////////////////////////////////////

// IMS compatible VIXFModel shell for ClosedFormIntegrateLN (varswap)
// and FourierEngine (convexity)  (VSCurve)
class BSVFutVSCurve: public VIXFModel {
public:
    static CClassConstSP const TYPE;

private:
    BSVFutVSCurve():VIXFModel(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(BSVFutVSCurve, clazz);
        SUPERCLASS(VIXFModel);
        EMPTY_SHELL_METHOD(defaultBSVFutVSCurve);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultBSVFutVSCurve(){
        return new BSVFutVSCurve();
    }
};

CClassConstSP const BSVFutVSCurve::TYPE =
CClass::registerClassLoadMethod(
    "BSVFutVSCurve", typeid(BSVFutVSCurve), load);


//////////////////////////////////////////////////////////////////////////////////////////////

/** Fourier payoff: product */
class VAssetFP:  virtual public FourierProduct,
                 virtual public FourierProductIntegrator1D,
                 virtual public FwdStFourierProductQuadVar,  //for close form in fourier for fwd start var swap.
                 virtual public FwdStFourierProductExpQuadVar{
public:
    /** equivalent to InstIntoFourierProduct */
    VAssetFP(const VAsset* inst, const VAssetAlgorithm* data, const DateTime& s, const DateTime& matDate); //s : startDate

    /** Constructs integrands for var option */
    Function1DDoubleArrayConstSP Integrand(const FourierEngine * model,
                                           const Integrator1D*  integrator);

    const DateTime& getStartDate() const {
        return startDate;
    }

    /** Post process method for VAssetFP*/
    virtual void postResults(const FourierEngine* model,
                             const Integrator1D*  integrator,
                             const FourierProductIntegrator1D::IntegralArray& integrals,
                             CControl*            control,
                             CResults*            results);

private:
    const VAsset*   inst;
    DateTime        maturity;
    DateTime        startDate;
};

///////////////////////////// VAsset ////////////////////////////////

/**
* PARVOL_SWAPVAR = sqt(E[Relized Var(T,T+x)]_0)
* PARVOL_SWAPVAR_SV = sqt(E[Relized Var(T,T+x)]_0) in SV
* SQT_SWAPVAR_SV = E{sqt[E(Relized Var(T,T+x)]_0}  in SV
* PARVOL_SWAPVAR_VIX_BASIS0 = difference between parvol_swapvar and VIX at value date
*/
const string VAsset::PARVOL_SWAPVAR = "parvol_swapvar";
const string VAsset::PARVOL_SWAPVAR_SV = "parvol_swapvar_sv";
const string VAsset::SQT_SWAPVAR_SV = "sqt_swapvar_sv";
const string VAsset::PARVOL_SWAPVAR_VIX_BASIS0 = "parvol_swapvar_vix_basis0";


VAsset::~VAsset(){}

string VAsset::sensName(Spot*) const {
    return name;
}
/* AD can't handle it, so comment for now
TweakOutcome VAsset::sensShift(const PropertyTweak<Spot>& shift) {
    static const string method = "VAsset::sensShift";

    double origSpot = spot;
    try {
        if (!Maths::isZero(shift.coefficient)){ // avoid altering if zero shift size
            spot = origSpot * (1 + shift.coefficient);       // shifted, as of valueDate
        }
    } catch(exception& e){
        spot = origSpot;
        throw ModelException(e, method);
    }

    // Arg1 is the initial spot to be stored on the delta shift
    // Arg2-Arg1 is the required divisor
    return TweakOutcome(origSpot, origSpot*(1+shift.coefficient), true);   // our components has a delta type sensitivity
}
*/
/** Restores the object to its original form */
/*
void VAsset::sensRestore(const PropertyTweak<Spot>& shift) {
    static const string method = "VAsset::sensRestore";

    double origSpot = spot;
    try {
        if (!Maths::isZero(shift.coefficient)) {
            spot = origSpot / (1 + shift.coefficient);
        }
    }
    catch(exception& e){
        spot = origSpot;
        throw ModelException(e, method);
    }
}
*/
/** Pull out the asset from the market data */
void VAsset::getMarket(const IModel* model, const MarketData* market)  {
    static const string method("VAsset::getMarket");
    try {
        const VIXFModel* vsModel = dynamic_cast<const VIXFModel*>(model);

        if (vsModel) {
            // retrieve asset data with BS vol, and history etc.
            IndexSpecEQ::setup(vsModel->varSwapModel.get(), market);

            // retrieve asset data with sv vol, just for vol
            undSV = CAssetWrapper(asset->getName());
            undSV.getData(vsModel->capModel.get(), market);
		}
		else {
            throw ModelException("VAsset::getMarket:" "Only VIXFModel supported");
		}

        // store algorithm data
        CControlSP ctrl(Control::makeFromFlags("P", 0.0)); // price only
        data = VAssetAlgorithmSP(new VAssetAlgorithm(VIXFModelConstSP::attachToRef(vsModel), ctrl, PARVOL_SWAPVAR));

        //basis = VIXSpot - par vol of Var Swap BS
        basisBS = spot - getBSSpot();
        basisSV = spot - getSVSpot();
    }
    catch (exception& e){
        throw ModelException(e, method, "Failed for VAsset ");
    }
}

/** returns the spot price */
/** if MTM = true, Market spot price; else, calculated spot price*/
double VAsset::getSpot() const{
    if (useMTM) {
        return spot;   //input market VIX spot
    }
    else{
        //calculate the VIX spot by Var swap BS
        return getBSSpot();
    }
}

/** returns the market spot price */
double VAsset::getMarketSpot() const{
    return spot;   //input market VIX spot
}


double VAsset::getBSSpot()const{
    //calculate the VIX spot by Var swap BS
    data->setType(VAsset::PARVOL_SWAPVAR);
    return priceBS(data->model->getValueDate())*100.0;
}

double VAsset::getSVSpot()const{
    //calculate the VIX spot by Var swap SV 
    data->setType(VAsset::PARVOL_SWAPVAR_SV);
    return priceSV(data->model->getValueDate())*100.0;
}

/** Basis defined as the diff between Market value and model value for VAsset spot*/
double VAsset::getBasis(const string basisType) const {

    if (CString::equalsIgnoreCase(basisType, VIXFModel::BS_SV)){
        return basisBS;
    }else if(CString::equalsIgnoreCase(basisType, VIXFModel::SV)){
        return basisSV;
    }else{
        return 0.0;
        //throw ModelException("VAsset::getBasis:" "basis Type must be BS_SV or SV.");
    }    
}

double VAsset::getFixing (DateTime& t, string source) const {

    AssetHistoryConstSP assetHistory = history->getAssetHistory(source) ;
    CashFlowArraySP samples = CashFlowArraySP(new CashFlowArray(0));
    CashFlow dateOnly(t, 0);
    samples->push_back(dateOnly);
    assetHistory.get()->getSamples(samples, t);
    double fixing = (*samples)[0].amount;

    return fixing;
}

/** vector "out" contains BS, SV components, only has values after calling vAsset->fwdValue*/
/** record "out" into results based on the control */
const void VAsset::recordExtraOutput(Control*   control,  
                          CResults*  results) const{

    if (control && control->isPricing()) {
        OutputRequest* request ;
        request = control->requestsOutput(OutputRequest::PARVOL_SWAPVAR);
        if (request) {
            //calc fwd Var Swap BS
            results->storeRequestResult(request, out[0]); //parVolSwapVar
        }
        request =control->requestsOutput(OutputRequest::PARVOL_SWAPVAR_SV);
        if (request) {
            //calc fwd Var Swap SV
            results->storeRequestResult(request, out[1]); //parVolSwapVarSV
        }
        request = control->requestsOutput(OutputRequest::SQT_SWAPVAR_SV);
        if (request) {
            //calc fwd Var E(sqrt(K^2(t,x)] using SV
            results->storeRequestResult(request, out[2]); //sqtSwapVarSV
        }
        request = control->requestsOutput(OutputRequest::PARVOL_SWAPVAR_VIX_BASIS0);
        if (request) {
            results->storeRequestResult(request, out[3]); //parVolSwapvarVixBasis0
        }
    }
}

/** Calculates the expected spot price of the asset at the given date */
double VAsset::fwdValue(const DateTime& date) const{
    if (useMTM) {
        return IndexSpecEQ::fwdValue(date);
    }else{
        return fwdValue(date,  const_cast<VAsset*>(this)->out);
    }
}

/** This is the method responsible for pricing the Fwd VarSwap 
    and the convexty adj using SVJ and aggregating the results
*/
/** compute fwd price according to model flag. all component values 
    are output to prices */
double VAsset::fwdValue(const DateTime& maturity,  DoubleArray& out) const{
    out.resize(4);
    //calc fwd Var Swap BS
    getAlgorithm()->setType(VAsset::PARVOL_SWAPVAR);
    out[0] = priceFwd(maturity);  //parVolSwapVar

    //calc fwd Var Swap SV
    getAlgorithm()->setType(VAsset::PARVOL_SWAPVAR_SV);
    out[1] = priceFwd(maturity);  //parVolSwapVarSV
    
    //calc fwd Var E(sqrt(K^2(t,x)] using SV
    getAlgorithm()->setType(VAsset::SQT_SWAPVAR_SV);
    out[2] = priceFwd(maturity);   //sqtSwapVarSV

    const string pM = getAlgorithm()->model->getPricingMethod();
    out[3] = getBasis(pM)/100.0;  //parVolSwapvarVixBasis0

    /** 
    BS_SV : BS + convexity + basis
    SV: E_t[sqrt(Fwd Var(T,T+ tau))] + basis
    BS_SV_NO_BASIS: BS + convexity
    SV_NO_BASIS: E_t[sqrt(Fwd Var(T,T+ tau))]
    VSW_BS: Var Swap using BS
    VSW_SV: Var Swap using SVJ
    */

    double totalValue;
    if (CString::equalsIgnoreCase(pM, VIXFModel::BS_SV)){        
        totalValue = out[0] + out[2] -out[1] + out[3];
    }else if(CString::equalsIgnoreCase(pM, VIXFModel::SV)){
        totalValue = out[2] + out[3];
    }else if(CString::equalsIgnoreCase(pM, VIXFModel::BS_SV_NO_BASIS)){
        totalValue = out[0] + out[2] -out[1] ;
    }else if(CString::equalsIgnoreCase(pM, VIXFModel::SV_NO_BASIS)){
        totalValue = out[2] ;
    }else if(CString::equalsIgnoreCase(pM, VIXFModel::VSW_BS)){
        totalValue = out[0] ;
    }else if(CString::equalsIgnoreCase(pM, VIXFModel::VSW_SV)){
        totalValue = out[1] ;
    }else{
        throw ModelException("pricingMethod chosen is wrong");
    }
    return totalValue;
}

/** price each components based on data type, need to set the data->type before calling*/
double VAsset::priceFwd(const DateTime& date) const{
    if (data->type == PARVOL_SWAPVAR)
        return priceBS(date);
    else if ((data->type == PARVOL_SWAPVAR_SV) || (data->type == SQT_SWAPVAR_SV))
        return priceSV(date);
    else{
        throw ModelException("VAsset::fwdValue:" "This type of fwdValue: " +
                                data->type + " isn't supported.");
    }
}

/** return BS model price */
double VAsset::priceBS(const DateTime& startDate) const{

    // transfer tenorPeriod  = nb of calandar days, should be 30d
    DateTime vswapMat =   tenor->toDate(startDate);
    
    HolidayConstSP hol = AssetUtil::getHoliday(asset.get());

    double parVolSwapVar;

    //if startDate + 30D is a holiday, then do interpolation on var between the two near bus day var
    if (hol->isHoliday(vswapMat)) {
        //look for the first date before t which isn't holiday
        DateTime noHolBefore;
        //look for the first date after t which isn't holiday
        DateTime noHolAfter;
        
        int i = 0;
        DateTime t = vswapMat;
        while (hol->isHoliday(t)){
            t = vswapMat.rollDate(i-1);
            i--;
        }
        noHolBefore = t;

        t = vswapMat;
        i = 0;
        while (hol->isHoliday(t)){
            t = vswapMat.rollDate(i+1);
            i++;
        }
        noHolAfter = t;

        double var1 =calcPriceBS(startDate, noHolBefore);
        double var2 =calcPriceBS(startDate, noHolAfter);
        double daysdiff1 = noHolBefore.yearFrac(vswapMat);
        double daysdiff2 = vswapMat.yearFrac(noHolAfter);

        if (Maths::isZero(daysdiff2+daysdiff1)){
            throw ModelException("VAsset::priceBS:" "shouldn't be called here! " );
            
        }else{
            parVolSwapVar = var1 + (var2 - var1) *daysdiff1 / (daysdiff2+daysdiff1);
        }
    }else{
        parVolSwapVar = calcPriceBS(startDate, vswapMat);
    }

    double nbDaysInVIX = vswapMat.daysDiff(startDate);

    //VIX is calculated using seconds
    double totalYears = startDate.yearFrac(vswapMat);

    parVolSwapVar = sqrt(parVolSwapVar/totalYears);

    return parVolSwapVar;
}

/** return BS model variance sigma^2*T */
double VAsset::calcPriceBS(const DateTime& startDate, const DateTime& matDate) const{

    VarianceSwapSP vswap = createVarSwap(asset, yc, data.get(), startDate, matDate);
	vswap->notional = 0.01;  //notional of VIX future

    Results results;
    // Price the var swap
    data->model->varSwapModel->Price(vswap.get(), data->ctrl.get(), &results);

    //value for variance swap based on variance swap model (close form....)
    double pv = vswap->instSettle->pv(vswap->valueDate,matDate,vswap->discount.get(),vswap->asset.get());
//    double parVolSwapVar = sqrt(results.retrievePrice()/(totalYears  * pv));   //need to remove pv from price to get vol
    double parVolSwapVar = results.retrievePrice()/pv;   //need to remove pv from price to get vol

    return parVolSwapVar;
}

/** create a var swap instrument*/
VarianceSwapSP VAsset::createVarSwap(const CAssetWrapper& asset, const YieldCurveWrapper& yc,
                                     const VAssetAlgorithm* data, 
                                     const DateTime& startDate, const DateTime& matDate) const{
	// Create VarianceSwap inst with no strike and insuring no mean and no scale.
    // set ppy = 1 and numTotalReturns = 1 and do scaling with timeMetric yearFrac later
    CashFlowArray samples(2);
    samples[0].date = startDate;
    //samples[1].date = tenor->toDate(startDate);
    samples[1].date = matDate;
    samples[0].amount = 0.0;
    samples[1].amount = 0.0;
    bool fwdStarting = data->model->getValueDate() < startDate;

    VarianceSwapSP vswap(new VarianceSwap(
                        data->model->getValueDate(),
                        startDate,
                        fwdStarting,
                        0.0,                        //inst->initialSpot
                        "V",
                        settle,
                        InstrumentSettlementSP(   ),  // premium settlement
                        asset,
                        yc,
                        samples,
                        0.0,                        //strikeVol = 0.0
                        1,                          //observationsPerYear =252
                        false,                      //subtractMeanVol,
                        "FORWARD",                  //payoffType = "FORWARD";
                        true,                       //dontScaleByStrike,
                        2.5,                        //cap
                        true,                       // noDivAdj = True
                        false,                      // const bool divAdjOnExDate
                        true,                       //isVanilla = true, numTotalReturns will be taken as input
                        0,                          // numPastReturns = 0
                        1));                        // numTotalReturns = 1
    vswap->Validate();
    return vswap;
}

// reprot sensitive strikes through var swap instrument
DoubleArraySP VAsset::sensitiveStrikes(OutputNameConstSP outputName, const VIXFModel* model) const{
    VarianceSwapSP vswap = createVarSwap(asset, yc, data.get(), model->getValueDate(), tenor->toDate(data->model->getValueDate()));
    return vswap->getSensitiveStrikes(outputName, model->varSwapModel.get());
}

/** return stochastic vol model price*/
double VAsset::priceSV(const DateTime& startDate) const{

    double parVolSwapVar;

    // transfer tenorPeriod  = nb of calandar days, should be 30d
    DateTime vswapMat =   tenor->toDate(startDate);
    HolidayConstSP hol = AssetUtil::getHoliday(asset.get());

    if (hol->isHoliday(vswapMat)) {
        //look for the first date before t which isn't holiday
        DateTime noHolBefore;
        //look for the first date after t which isn't holiday
        DateTime noHolAfter;
        
        int i = 0;
        DateTime t = vswapMat;
        while (hol->isHoliday(t)){
            t = vswapMat.rollDate(i-1);
            i--;
        }
        noHolBefore = t;

        t = vswapMat;
        i = 0;
        while (hol->isHoliday(t)){
            t = vswapMat.rollDate(i+1);
            i++;
        }
        noHolAfter = t;

        double daysdiff1 = noHolBefore.yearFrac(vswapMat);
        double daysdiff2 = vswapMat.yearFrac(noHolAfter);

        double var1 =calcPriceSV(startDate, noHolBefore);
        double var2 =calcPriceSV(startDate, noHolAfter);

        if (Maths::isZero(daysdiff2+daysdiff1)){
            throw ModelException("VAsset::priceBS:" "shouldn't be called here! " );
            
        }else{
            parVolSwapVar = var1 + (var2 - var1) *daysdiff1 / (daysdiff2+daysdiff1);
        }
    }else{
        parVolSwapVar = calcPriceSV(startDate, vswapMat);
    }

    double nbDaysInVIX = vswapMat.daysDiff(startDate);
    double totalYears = startDate.yearFrac(vswapMat);

    parVolSwapVar = sqrt(parVolSwapVar/totalYears);

    return parVolSwapVar  ;
}

/** return stochastic vol model variance sigma^2*T*/
double VAsset::calcPriceSV(const DateTime& startDate, const DateTime& matDate) const{

    const FourierEngine* m = dynamic_cast<const FourierEngine*>(data->model->capModel.get());
    FourierProcess& process = const_cast<FourierProcess&>(m->getProcess());

    VAssetFP prod(this, data.get(), startDate, matDate);
    prod.validate(&process);

    Results results;
    prod.price(m, data->ctrl.get(), &results);

    return results.retrievePrice();
}

DateTime VAsset::endDate(const Sensitivity* sensControl) const {
    return tenor->toDate(data->model->getValueDate());
}

/** record forwards at maturity*/
void VAsset::recordFwdAtMat(OutputRequest*  request,
                         CResults*       results,
                         const DateTime& maturityDate) const
{
    //assuming that the fwd value for VIX is the 
    //E_0(sqrt(K^2(t,t+tau)] using pricingMethod in VIXFModel
    //t is the last date of the instrument, which = VIX Fut + 30 Days
  //put it into comments to speed 
    double fwd = fwdValue(maturityDate);  

    // wrap forward price into object
    CDoubleSP fwdObj(CDouble::create(fwd));

    results->storeRequestResult(request,
                                fwdObj,
                                OutputNameSP(new OutputName(getName())));


    // record forward at maturity for asset (ex: S&P500) BS
    asset->recordFwdAtMat(request, results, maturityDate);
}


/* for reflection */
VAsset::VAsset() :
    IndexSpecEQ(TYPE), 
    spot( 0. ),
    tenor(new MaturityPeriod("30D")),
    useMTM(false){

    settle = InstrumentSettlementSP(new CashSettlePeriod(0));
    out.resize(4, 0.0);
}

class VAssetHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VAsset, clazz);
        SUPERCLASS(IndexSpecEQ);
        //IMPLEMENTS(IRestorableWithRespectTo<Spot>);
        EMPTY_SHELL_METHOD(defaultVAsset);
        //optional
        FIELD(spot, "current index spot level");
        FIELD(tenor, "VIX tenor, 30 days as default value");
        FIELD_MAKE_OPTIONAL(tenor);
        FIELD(useMTM, "");
        FIELD_MAKE_OPTIONAL(useMTM);

        // transient fields
        FIELD(basisBS,"spot - varSwap price");
	    FIELD_MAKE_TRANSIENT(basisBS);
        FIELD(basisSV,"spot - varSwap price");
	    FIELD_MAKE_TRANSIENT(basisSV);
        FIELD(undSV,"asset with stochastic vol parameters");
	    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(undSV);
        FIELD(data,"");
	    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(data);
        FIELD(out,"");
	    FIELD_MAKE_TRANSIENT(out);
    }

    static IObject* defaultVAsset(){
        return new VAsset();
    }
};

CClassConstSP const VAsset::TYPE = CClass::registerClassLoadMethod(
    "VAsset", typeid(VAsset), VAssetHelper::load);

CClassConstSP const VAssetAlgorithm::TYPE = CClass::registerClassLoadMethod(
    "VAssetAlgorithm", typeid(VAssetAlgorithm), VAssetAlgorithm::load);

//////////////////////////////////////////////////////////////////////////////////

//**********************************************************************//
//****************   PRODUCT CLASS AND METHODS    **********************//
//**********************************************************************//
// ** using Fourier method for SVJ **
/** E[sqrt(K^2(t,x))]*/
class VAssetIntegrand:public Function1DDouble {
public:
    VAssetIntegrand(const FwdStFourierProcessExpQuadVar& process,
					   const FwdStFourierProductExpQuadVar& product,
					   const DateTime& maturity):
        Function1DDouble(Range(OpenBoundary(0.0), Infinity(Infinity::Plus))),
        process(process), product(product), matdate(maturity) {}

    double operator()(double  x) const {
			// based on Arnaud and Manos's presentation for VAsset
			return (1.0 - exp(process.cumulant(product, -x, matdate) ).real()) / (x * sqrt(x));
    }

private:
    const FwdStFourierProcessExpQuadVar& process;
    const FwdStFourierProductExpQuadVar& product;
    DateTime							 matdate;
};

//////////////////////////////////////////////////////////////////////////////////

/** equivalent to InstIntoFourierProduct
*/
VAssetFP::VAssetFP(const VAsset* inst, const VAssetAlgorithm* data, const DateTime& s, const DateTime& matDate):
        FourierProduct(inst->undSV.get(),
                        data->model->getValueDate(),
                        inst->yc.get(),
                        inst->settle.get() /* instrument settlement*/),
                        inst(inst){
    startDate = s;
    //maturity = inst->tenor->toDate(startDate);
    maturity = matDate;
}

/** Constructs integrands for VIX Future */
Function1DDoubleArrayConstSP VAssetFP::Integrand(const FourierEngine * model,
													 const Integrator1D*  integrator) {
    static const string method = "VAssetFP::Integrand";

    try{
        const FourierProcess& process = model->getProcess();
        const FwdStFourierProductExpQuadVar& thisProd = *this;
        const FwdStFourierProcessExpQuadVar* thisProc
                        = dynamic_cast<const FwdStFourierProcessExpQuadVar*>(&process);
        if(!thisProc) {
            throw ModelException(method, "Process does not support FwdStFourierProcessExpQuadVar interface." );
        }

		Function1DDoubleArraySP functions(new Function1DDoubleArray(1));
        (*functions)[0] = Function1DDoubleSP(
                        new VAssetIntegrand(*thisProc, thisProd, maturity));
        return functions;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Post process method for VAssetFP */
void VAssetFP::postResults(const FourierEngine* model,
                               const Integrator1D*  integrator,
                               const FourierProductIntegrator1D::IntegralArray& integrals,
                               CControl*            control,
                               CResults*            results) {
    static const string method = "VAssetFP::postResults";

    try {
        double value = 0.0;
        const FourierProcess& process = model->getProcess();
        //double nbDaysInVIX = maturity.daysDiff(startDate);
        //double totalYears = nbDaysInVIX / 365.0;
        double expectedVar;
        double years ;

        if (inst->data->type == VAsset::SQT_SWAPVAR_SV){
            double sqtSwapVarSV = integrals[0] / (2 * sqrt(Maths::PI));
            //sqtSwapVarSV /= sqrt(totalYears);
            //value = sqtSwapVarSV;
            //total var
            value = sqtSwapVarSV*sqtSwapVarSV;

            // See if future price is positive
            if(sqtSwapVarSV <= 0.0){
                throw ModelException("The value of sqt Swap Var SV is negative " +
                                    Format::toString(sqtSwapVarSV) +
                                    ". Check integrator parameters e.g. tolerance, nbSteps.");
            }
        }else if(inst->data->type == VAsset::PARVOL_SWAPVAR_SV){
			const FwdStFourierProductQuadVar& thisProd = *this;
			const FwdStFourierProcessQuadVar* thisProc = dynamic_cast<const FwdStFourierProcessQuadVar*>(&process);
			if(!thisProc) {
				throw ModelException(method, "Process does not support FwdStFourierProcessQuadVar interface" );
			}
			// Expectation prices annualized var from startDate to startDate + tenor
			expectedVar =  thisProc->expectation(thisProd, maturity);
			years = process.getTimeMetric().yearFrac(startDate, maturity);
			expectedVar *= years;

            // Vix is defined in terms of 30 / 365
            //double parVolSwapVarSV = sqrt(expectedVar) / sqrt(totalYears);
            double parVolSwapVarSV = expectedVar ;
            value =  parVolSwapVarSV;
        }else{
            throw ModelException(method, "This type: " +
                                 inst->data->type + " of calculation isn't supported.");
        }

        //no pv is needed for future
        results->storePrice(value, discount->getCcy());
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

//////////////////////////////////////////////////////////////////////////////////

bool VAssetLoad()
{
    return (VAsset::TYPE != 0);
}

DRLIB_END_NAMESPACE

