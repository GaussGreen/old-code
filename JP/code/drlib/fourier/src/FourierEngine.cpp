//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierEngine.cpp
//
//   Description : Fourier Engine Integrator
//
//   Date        : February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/Addin.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/FourierProcessSV.hpp"
#include "edginc/FourierProcessSVJ.hpp"
#include "edginc/FourierProcessSVCJ.hpp"
#include "edginc/FourierProcessVSCurve.hpp"


DRLIB_BEGIN_NAMESPACE

// INTEGRATOR TYPE HELPER
typedef FourierEngine::IntegratorType FEIT;
template<> string nameForType<FEIT>(FEIT*){
    return "FourierEngine::IntegratorType";
}

typedef Enum2StringListHelper<FEIT> FEITHelper;
template<> string FEITHelper::names[FEITHelper::EnumList::NB_ENUMS] = {
    "INTEGRAL",
    "FFT"
};

// ISAP
void FourierEngine::ISAP::load(CClassSP &clazz){
    REGISTER_INTERFACE(ISAP, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}

CClassConstSP const FourierEngine::ISAP::TYPE = CClass::registerInterfaceLoadMethod(
    "FourierEngine::ISAP", typeid(FourierEngine::ISAP), FourierEngine::ISAP::load);


MarketDataFetcherSP FourierEngine::createMDF() const {
    return process->marketDataFetcher();
}

void FourierEngine::Price(CInstrument*  instrument,
                          CControl*     control,
                          CResults*     results)
{
    static const string method = "FourierEngine::Price";
    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support FourierEngine::IntoProduct");
    }
    try{
        // valueDate >= matDate is taken care of here
        //if (instrument->priceDeadInstrument(control, results)) {
        //    return; // done for a dead instrument
        //}

        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        auto_ptr<FourierProduct> prodAP(intoProd.createProduct(this));
        prodAP->validate(process.get());
        prodAP->price(this, control, results);
    } catch (exception& e){
        throw ModelException(e, method, "Failed to createProduct for "
                             "instrument of type "+
                             instrument->getClass()->getName());
    }
}


void FourierEngine::validatePop2Object(){
    static const string method = "FourierEngine::validatePop2Object";
    integratorTypeIndex = FEITHelper::getIndex(integratorType);
}

const FourierProcess& FourierEngine::getProcess() const {
    return *process;
}

IModel::WantsRiskMapping FourierEngine::wantsRiskMapping() const {
    return process->wantsRiskMapping();
}

const FourierEngine::ISAP& FourierEngine::getISAP() const {
    if (!isap){
        throw ModelException("FourierEngine::getISAP", "isap is missing");
    }
    return *isap;
}

FourierEngine::FourierEngine():
CModel(TYPE),
integrator(new IntFuncInf()),
fftIntegrator(new FFTIntegrator1D()),
integratorType(FEITHelper::getDefaultName()){}

FourierEngine::FourierEngine(CClassConstSP clazz):
CModel(clazz),
integrator(new IntFuncInf()),
fftIntegrator(new FFTIntegrator1D()),
integratorType(FEITHelper::getDefaultName()){}

FourierEngine::FourierEngine(const FourierProcessSP&  p,
                             const string&            integratorType,
                             const Integrator1DSP&    integrator,
                             const FFTIntegrator1DSP& fftIntegrator):
CModel(TYPE),
process(p),
integrator(integrator),
fftIntegrator(fftIntegrator),
integratorType(integratorType){
    validatePop2Object();
}

FourierEngine::FourierEngine(const FourierProcessSP&  p):
CModel(TYPE),
process(p),
integrator(IntFuncInf::create()),
fftIntegrator(FFTIntegrator1D::create()),
integratorType(FEITHelper::getDefaultName()){
    validatePop2Object();
}

class FourierEngineHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FourierEngine, clazz);
        SUPERCLASS(CModel);
        EMPTY_SHELL_METHOD(defaultFourierEngine);
        FIELD(process, "Fourier Process");
        FIELD(isap, "Instrument Specific Algorithm Parameters");
        FIELD_MAKE_OPTIONAL(isap);
        FIELD(integrator, "Integrator");
        FIELD_MAKE_OPTIONAL(integrator);
        FIELD(fftIntegrator, "FFT Integrator");
        FIELD_MAKE_OPTIONAL(fftIntegrator);
        FIELD(integratorType, "Integrator Type ("+ FEITHelper::getNameList() +")");
        FIELD_MAKE_OPTIONAL(integratorType);
        FIELD(integratorTypeIndex, "");
        FIELD_MAKE_TRANSIENT(integratorTypeIndex);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultFourierEngine(){
        return new FourierEngine();
    }
};

CClassConstSP const FourierEngine::TYPE = CClass::registerClassLoadMethod(
    "FourierEngine", typeid(FourierEngine), FourierEngineHelper::load);

CClassConstSP const FourierEngine::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("FourierEngine::IIntoProduct",
                                    typeid(FourierEngine::IIntoProduct), 0);

string FourierProductIntegrator1D::errorHandling(int                  iIntegrand,
                                                 const Integrator1D*  integrator) const {
    string errormsg = Format::toString("Failed to compute integrand number %ld.", iIntegrand + 1);
    return errormsg;
}

string FourierProductFFTIntegrator1D::errorHandling(int                    iIntegrand,
                                                    const FFTIntegrator1D* integrator) const {
    string errormsg = Format::toString("Failed to compute integrand number %ld.", iIntegrand + 1);
    return errormsg;
}

void FourierProduct::validate(FourierProcess* process) const{
    const static string method("FourierProduct::validate");
    int i = 0;
    for(; i < mAsset->NbAssets(); i++) {
        string CcyTtreatment(mAsset->assetGetCcyTreatment(i));
        if(CcyTtreatment == "S" || CcyTtreatment == "P") {
            throw ModelException(method,
                                 "Currency Protected or Struck assets not supported.");
        }
    }
    // give the process the chance to do some extra validation
    // and to initialize its transient fields
    process->validate(this);
}

/** Gives a point in the interior of the intersection of two regions */
double FourierProduct::intersectRanges(double l1,
                                       double u1,
                                       double l2,
                                       double u2,
                                       double p) {
    const static string method("FourierProduct::intersectRanges");

    if(!Maths::isPositive(p) || !Maths::isPositive(1.0 - p)){
        throw ModelException(method,
                             "p (" + Format::toString(p) + ") must be strictly between 0.0 and 1.0.");
    }

    if(l1 > u1 || l2 > u2) {
        throw ModelException(method,
                             Format::toString("Order must be l1 <= u1, l2 <= u2; got %f, %f, %f and %f respectively.",
                                              l1, u1, l2, u2));
    }

    double l = Maths::max(l1, l2);
    double u = Maths::min(u1, u2);

    if(!Maths::isNegative(l - u)) {
        throw ModelException(method,
                             Format::toString("Intervals [%f, %f] and [%f, %f] do not intersect or just touch.",
                                              l1, u1, l2, u2));
    }

    return (p * l + (1.0 - p) * u);
}

/** price method */
void FourierProduct::recordFwdAtMat(Control*             control,
                                    const DateTime&      matDate,
                                    Results*             results){
    const static string method("FourierProduct::recordFwdAtMat");
    try{
        OutputRequest* request = 0;
        if (control && control->isPricing() &&
            (request = control->requestsOutput(OutputRequest::FWD_AT_MAT))) {
                if (matDate.isGreaterOrEqual(today)){
                    mAsset->recordFwdAtMat(request, results, matDate);
                }
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}


const IMultiFactors& FourierProduct::getMultiAsset() const{
    return *mAsset;
}

/** price method */
void FourierProduct::price(const FourierEngine* model,
                           Control*             control,
                           Results*             results){
    const static string method("FourierProduct::price");

    try{
        if (model->integratorTypeIndex == FEITHelper::EnumList::INTEGRAL){
            FourierProductIntegrator1D* thisProd = dynamic_cast<FourierProductIntegrator1D*>(this);
            if (!thisProd){
                throw ModelException(method, Format::toString(typeid(this).name()) + " does not support Integrator1D's");
            }

            Integrator1D* thisIntegrator = model->integrator.get();
            if (!thisIntegrator){
                throw ModelException(method, "A Integrator1D-type integrator is required");
            }

            Function1DDoubleArrayConstSP myfuncs(thisProd->Integrand(model, thisIntegrator));

            int NbIntegrals = myfuncs->size();
            FourierProductIntegrator1D::IntegralArray integrals(NbIntegrals);

            int iIntegral = 0;
            for(; iIntegral < NbIntegrals; iIntegral++) {
                try {
                    integrals[iIntegral] = thisIntegrator->integrate(*(*myfuncs)[iIntegral]);
                } catch(exception& e) {
                    try {
                        // Be careful as this might fail
                        string payoffMessage  = thisProd->errorHandling(iIntegral, thisIntegrator);
                        string processMessage = model->process->getParameters();
                        string message = payoffMessage + "\n" + processMessage;
                        throw ModelException::addTextToException(e, message);
                    } catch(exception &ee) {
                        // Just throw the original exception
                        throw ModelException(e, method);
                    }
                }
            }

            thisProd->postResults(model,
                                  thisIntegrator,
                                  integrals,
                                  control,
                                  results);
        }
        else {  // FFT
            FourierProductFFTIntegrator1D* thisProd = dynamic_cast<FourierProductFFTIntegrator1D*>(this);
            if (!thisProd){
                throw ModelException(method, Format::toString(typeid(this).name()) + " does not support FFTIntegrator1D's");
            }

            FFTIntegrator1D* thisIntegrator = model->fftIntegrator.get();
            if (!thisIntegrator){
                throw ModelException(method, "A FFTIntegrator1D-type integrator is required");
            }

            Function1DComplexArrayConstSP myfuncs(thisProd->Integrand(model, thisIntegrator));

            int NbIntegrals = myfuncs->size();
            FourierProductFFTIntegrator1D::IntegralArray integrals(NbIntegrals);

            int iIntegral = 0;
            for(; iIntegral < NbIntegrals; iIntegral++) {
                try {
                    integrals[iIntegral] = thisIntegrator->integrate(*(*myfuncs)[iIntegral]);
                } catch(exception& e) {
                    try {
                        // Be careful as this might fail
                        string payoffMessage  = thisProd->errorHandling(iIntegral, thisIntegrator);
                        string processMessage = model->process->getParameters();
                        string message = payoffMessage + "\n" + processMessage;
                        throw ModelException::addTextToException(e, message);
                    } catch(exception &ee) {
                        // Just throw the original exception
                        throw ModelException(e, method);
                    }
                }
            }

            thisProd->postResults(model,
                                  thisIntegrator,
                                  integrals,
                                  control,
                                  results);
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// turns Asset into MultiFactors
IMultiFactorsConstSP FourierProduct::convertToMultiFactors(const Asset* asset){
    IMultiFactorsConstSP mAsset(asset->asMultiFactors());
    return mAsset;
}

/** 'Full' constructor for single factor payoffs */
FourierProduct::FourierProduct(
    const CAsset*               asset,
    const DateTime&             today,
    const YieldCurve*           discount,
    const InstrumentSettlement* instSettle):
mAsset(convertToMultiFactors(asset)),
today(today),
discount(discount),
instSettle(instSettle),
NbAssets(mAsset->NbAssets()){}

/** 'Full' constructor for N factor payoffs */
FourierProduct::FourierProduct(
    const IMultiFactors*        mFactors,
    const DateTime&             today,
    const YieldCurve*           discount,
    const InstrumentSettlement* instSettle):
mAsset(IMultiFactorsConstSP::attachToRef(mFactors)),
today(today),
discount(discount),
instSettle(instSettle),
NbAssets(mAsset->NbAssets()){}

// Trivial FwdStarting product
class TrivialFwdStFourierProductLogRtn: public FwdStFourierProductLogRtn {
public:
    TrivialFwdStFourierProductLogRtn(const DateTime& startDate): startDate(startDate) {};
    virtual const DateTime& getStartDate() const {
        return startDate;
    }

private:
    DateTime startDate;
};

// Trivial FwdStarting product
class TrivialFwdStFourierProductIntVar: public FwdStFourierProductIntVar {
public:
    TrivialFwdStFourierProductIntVar(const DateTime& startDate): startDate(startDate) {};
    virtual const DateTime& getStartDate() const {
        return startDate;
    }

private:
    DateTime startDate;
};

// Trivial FwdStarting product
class TrivialFwdStFourierProductQuadVar: public FwdStFourierProductQuadVar {
public:
    TrivialFwdStFourierProductQuadVar(const DateTime& startDate): startDate(startDate) {};
    virtual const DateTime& getStartDate() const {
        return startDate;
    }

private:
    DateTime startDate;
};

// Addin that compute the Laplace transform of log returns
class LaplaceTransformAddin: public CObject{
    static CClassConstSP const TYPE;

    class DummyProduct: public FourierProduct{
    public:
        DummyProduct(const CAsset*               asset,        // single factor
                     const DateTime&             today):       // value date
        FourierProduct(asset, today, 0, 0){}
    };

    CMarketDataSP    market;
    CAssetSP         asset;
    FourierProcessSP process;
    DoubleArray      realFrequency;
    DoubleArray      imagFrequency;
    DateTimeArray    maturities;
    DateTime         startDate;

    void validatePop2Object() {
        static const string routine = "LaplaceTransformAddin::validatePopToObject";
        try {
            if (startDate.empty()){
                startDate = market->GetReferenceDate();
            }

            DateTime::ensureIncreasing(maturities,
                                       "Maturities",
                                       true);

            if(maturities[0] <= startDate) {
                throw ModelException("First date " +
                                     maturities[0].toString() +
                                     " is before start date " +
                                     startDate.toString());
            }
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

    static IObjectSP computeLaplaceReturn(LaplaceTransformAddin* params){
        static const string routine = "LaplaceTransformAddin::computeLaplaceReturn";
        try {
            FourierEngine model(params->process);
            params->asset->getMarket(&model, params->market.get());
            // need to create a dummy model to grab vol object
            DummyProduct product(params->asset.get(), params->startDate);
            product.validate(params->process.get());

            int n = params->realFrequency.size();
            int nbMaturities = params->maturities.size();
            CDoubleMatrixSP realLaplace(new CDoubleMatrix(nbMaturities, n));
            CDoubleMatrixSP imagLaplace(new CDoubleMatrix(nbMaturities, n));
            CDoubleMatrixSP modulusLaplace(new CDoubleMatrix(nbMaturities, n));

            double* realFrequency = &params->realFrequency[0];
            double* imagFrequency = &params->imagFrequency[0];

            if(params->startDate <= params->market->GetReferenceDate()) {
                const StFourierProcessLogRtn* stReturn =
                    dynamic_cast<const StFourierProcessLogRtn*>(params->process.get());
                if(!stReturn) {
                    throw ModelException("Process does not support StFourierProcessLogRtn interface");
                }

                StFourierProductLogRtn trivialStProduct;
                for(int iMat = 0; iMat < nbMaturities; iMat++) {
                    double* realLaplaceMat = (*realLaplace)[iMat];
                    double* imagLaplaceMat = (*imagLaplace)[iMat];
                    double* modulusLaplaceMat = (*modulusLaplace)[iMat];
                    for(int i = 0; i < n; i++) {
                        Complex z(realFrequency[i], imagFrequency[i]);
                        Complex cumulant = stReturn->scalelessCumulant(trivialStProduct,
                                                                       z,
                                                                       params->maturities[iMat]);
                        Complex laplace = exp(cumulant);
                        realLaplaceMat[i] = laplace.real();
                        imagLaplaceMat[i] = laplace.imag();
                        modulusLaplaceMat[i] = sqrt(Maths::square(realLaplaceMat[i]) +
                                                    Maths::square(imagLaplaceMat[i]));
                    }
                }
            } else {
                const FwdStFourierProcessLogRtn* fwdReturn =
                    dynamic_cast<const FwdStFourierProcessLogRtn*>(params->process.get());
                if(!fwdReturn) {
                    throw ModelException("Process does not support FwdStFourierProcessLogRtn interface");
                }

                TrivialFwdStFourierProductLogRtn trivialFwdStProduct(params->startDate);

                for(int iMat = 0; iMat < nbMaturities; iMat++) {
                    double* realLaplaceMat = (*realLaplace)[iMat];
                    double* imagLaplaceMat = (*imagLaplace)[iMat];
                    double* modulusLaplaceMat = (*modulusLaplace)[iMat];
                    for(int i = 0; i < n; i++) {
                        Complex z(realFrequency[i], imagFrequency[i]);
                        Complex cumulant = fwdReturn->scalelessCumulant(trivialFwdStProduct,
                                                                        z,
                                                                        params->maturities[iMat]);
                        Complex laplace = exp(cumulant);
                        realLaplaceMat[i] = laplace.real();
                        imagLaplaceMat[i] = laplace.imag();
                        modulusLaplaceMat[i] = sqrt(Maths::square(realLaplaceMat[i]) +
                                                    Maths::square(imagLaplaceMat[i]));
                    }
                }
            }

            realLaplace->transpose();
            imagLaplace->transpose();
            modulusLaplace->transpose();
            ObjectArraySP output(new ObjectArray(3));
            (*output)[0] = realLaplace;
            (*output)[1] = imagLaplace;
            (*output)[2] = modulusLaplace;

            return IObjectSP(output.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    //computes LaplaceTransform IntVar
    static IObjectSP computeLaplaceIntVar(LaplaceTransformAddin* params){
        static const string routine = "LaplaceTransformAddin::computeLaplaceIntVar";
        try {
            FourierEngine model(params->process);
            params->asset->getMarket(&model, params->market.get());
            // need to create a dummy model to grab vol object
            DummyProduct product(params->asset.get(), params->startDate);
            product.validate(params->process.get());

            int n = params->realFrequency.size();
            int nbMaturities = params->maturities.size();
            CDoubleMatrixSP realLaplace(new CDoubleMatrix(nbMaturities, n));
            CDoubleMatrixSP imagLaplace(new CDoubleMatrix(nbMaturities, n));
            CDoubleMatrixSP modulusLaplace(new CDoubleMatrix(nbMaturities, n));

            double* realFrequency = &params->realFrequency[0];
            double* imagFrequency = &params->imagFrequency[0];

            if(params->startDate <= params->market->GetReferenceDate()) {
                const StFourierProcessIntVar* stVar =
                    dynamic_cast<const StFourierProcessIntVar*>(params->process.get());
                if(!stVar) {
                    throw ModelException("Process does not support StFourierProcessIntVar interface");
                }

                StFourierProductIntVar trivialStProduct;
                for(int iMat = 0; iMat < nbMaturities; iMat++) {
                    double* realLaplaceMat = (*realLaplace)[iMat];
                    double* imagLaplaceMat = (*imagLaplace)[iMat];
                    double* modulusLaplaceMat = (*modulusLaplace)[iMat];
                    for(int i = 0; i < n; i++) {
                        Complex z(realFrequency[i], imagFrequency[i]);
                        Complex cumulant = stVar->cumulant(trivialStProduct,
                                                           z,
                                                           params->maturities[iMat]);
                        Complex laplace = exp(cumulant);
                        realLaplaceMat[i] = laplace.real();
                        imagLaplaceMat[i] = laplace.imag();
                        modulusLaplaceMat[i] = sqrt(Maths::square(realLaplaceMat[i]) +
                                                    Maths::square(imagLaplaceMat[i]));
                    }
                }
            } else {
                const FwdStFourierProcessIntVar* fwdVar =
                    dynamic_cast<const FwdStFourierProcessIntVar*>(params->process.get());
                if(!fwdVar) {
                    throw ModelException("Process does not support FwdStFourierProcessIntVar interface");
                }

                TrivialFwdStFourierProductIntVar trivialFwdStProduct(params->startDate);

                for(int iMat = 0; iMat < nbMaturities; iMat++) {
                    double* realLaplaceMat = (*realLaplace)[iMat];
                    double* imagLaplaceMat = (*imagLaplace)[iMat];
                    double* modulusLaplaceMat = (*modulusLaplace)[iMat];
                    for(int i = 0; i < n; i++) {
                        Complex z(realFrequency[i], imagFrequency[i]);
                        Complex cumulant = fwdVar->cumulant(trivialFwdStProduct,
                                                            z,
                                                            params->maturities[iMat]);
                        Complex laplace = exp(cumulant);
                        realLaplaceMat[i] = laplace.real();
                        imagLaplaceMat[i] = laplace.imag();
                        modulusLaplaceMat[i] = sqrt(Maths::square(realLaplaceMat[i]) +
                                                    Maths::square(imagLaplaceMat[i]));
                    }
                }
            }

            realLaplace->transpose();
            imagLaplace->transpose();
            modulusLaplace->transpose();
            ObjectArraySP output(new ObjectArray(3));
            (*output)[0] = realLaplace;
            (*output)[1] = imagLaplace;
            (*output)[2] = modulusLaplace;

            return IObjectSP(output.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    //computes Laplace transform Quad Var
    static IObjectSP computeLaplaceQuadVar(LaplaceTransformAddin* params){
        static const string routine = "LaplaceTransformAddin::computeLaplaceQuadVar";
        try {
            FourierEngine model(params->process);
            params->asset->getMarket(&model, params->market.get());
            // need to create a dummy model to grab vol object
            DummyProduct product(params->asset.get(), params->startDate);
            product.validate(params->process.get());

            int n = params->realFrequency.size();
            int nbMaturities = params->maturities.size();
            CDoubleMatrixSP realLaplace(new CDoubleMatrix(nbMaturities, n));
            CDoubleMatrixSP imagLaplace(new CDoubleMatrix(nbMaturities, n));
            CDoubleMatrixSP modulusLaplace(new CDoubleMatrix(nbMaturities, n));

            double* realFrequency = &params->realFrequency[0];
            double* imagFrequency = &params->imagFrequency[0];

            if(params->startDate <= params->market->GetReferenceDate()) {
                const StFourierProcessQuadVar* stVar =
                    dynamic_cast<const StFourierProcessQuadVar*>(params->process.get());
                if(!stVar) {
                    throw ModelException("Process does not support StFourierProcessQuadVar interface");
                }

                StFourierProductQuadVar trivialStProduct;
                for(int iMat = 0; iMat < nbMaturities; iMat++) {
                    double* realLaplaceMat = (*realLaplace)[iMat];
                    double* imagLaplaceMat = (*imagLaplace)[iMat];
                    double* modulusLaplaceMat = (*modulusLaplace)[iMat];
                    for(int i = 0; i < n; i++) {
                        Complex z(realFrequency[i], imagFrequency[i]);
                        Complex cumulant = stVar->cumulant(trivialStProduct,
                                                           z,
                                                           params->maturities[iMat]);
                        Complex laplace = exp(cumulant);
                        realLaplaceMat[i] = laplace.real();
                        imagLaplaceMat[i] = laplace.imag();
                        modulusLaplaceMat[i] = sqrt(Maths::square(realLaplaceMat[i]) +
                                                    Maths::square(imagLaplaceMat[i]));
                    }
                }
            } else {
                const FwdStFourierProcessQuadVar* fwdVar =
                    dynamic_cast<const FwdStFourierProcessQuadVar*>(params->process.get());
                if(!fwdVar) {
                    throw ModelException("Process does not support FwdStFourierProcessQuadVar interface");
                }

                TrivialFwdStFourierProductQuadVar trivialFwdStProduct(params->startDate);

                for(int iMat = 0; iMat < nbMaturities; iMat++) {
                    double* realLaplaceMat = (*realLaplace)[iMat];
                    double* imagLaplaceMat = (*imagLaplace)[iMat];
                    double* modulusLaplaceMat = (*modulusLaplace)[iMat];
                    for(int i = 0; i < n; i++) {
                        Complex z(realFrequency[i], imagFrequency[i]);
                        Complex cumulant = fwdVar->cumulant(trivialFwdStProduct,
                                                            z,
                                                            params->maturities[iMat]);
                        Complex laplace = exp(cumulant);
                        realLaplaceMat[i] = laplace.real();
                        imagLaplaceMat[i] = laplace.imag();
                        modulusLaplaceMat[i] = sqrt(Maths::square(realLaplaceMat[i]) +
                                                    Maths::square(imagLaplaceMat[i]));
                    }
                }
            }

            realLaplace->transpose();
            imagLaplace->transpose();
            modulusLaplace->transpose();
            ObjectArraySP output(new ObjectArray(3));
            (*output)[0] = realLaplace;
            (*output)[1] = imagLaplace;
            (*output)[2] = modulusLaplace;

            return IObjectSP(output.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    LaplaceTransformAddin(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LaplaceTransformAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultLaplaceTransformAddin);
        FIELD(market, "Market");
        FIELD(asset, "Asset");
        FIELD(process, "Fourier process");
        FIELD(realFrequency, "Real part of frequency");
        FIELD(imagFrequency, "Im part of frequency");
        FIELD(maturities, "Maturities");
        FIELD(startDate, "Start Date");
        FIELD_MAKE_OPTIONAL(startDate);

        Addin::registerClassObjectMethod("LAPLACE_TRANSFORM",
                                         Addin::RISK,
                                         "The Laplace transform of the log-return",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)computeLaplaceReturn);

        Addin::registerClassObjectMethod("LAPLACE_TRANSFORM_INTVAR",
                                         Addin::RISK,
                                         "The Laplace transform of the integrated variance",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)computeLaplaceIntVar);

        Addin::registerClassObjectMethod("LAPLACE_TRANSFORM_QUADVAR",
                                         Addin::RISK,
                                         "The Laplace transform of the quadratic variation",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)computeLaplaceQuadVar);
    }

    static IObject* defaultLaplaceTransformAddin(){
        return new LaplaceTransformAddin();
    }
};

CClassConstSP const LaplaceTransformAddin::TYPE = CClass::registerClassLoadMethod(
    "LaplaceTransformAddin", typeid(LaplaceTransformAddin), load);

/*********************************************************************/

/* The ultimate wrapping of FourierProcess, mainly for use in Pyramid
 */
#define FOURIERPROCESS_TYPE_SV          "SV"
#define FOURIERPROCESS_TYPE_SVJ         "SVJ"
#define FOURIERPROCESS_TYPE_SVCJ        "SVCJ"
#define FOURIERPROCESS_TYPE_VSCURVE     "VSCURVE"

class FourierProcessWrapper : public CObject,
                              virtual public ITypeConvert {
public: // how can I have this protected or private?
    string                      processType;     // SV, SVJ, SVCJ or VSCurve
    FourierProcessSVSP          processSV;
    FourierProcessSVJSP         processSVJ;
    FourierProcessSVCJSP        processSVCJ;
    FourierProcessVSCurveSP     processVSCurve;

private:
    FourierProcessSP            realProcess;

public:
    static CClassConstSP const TYPE;

    // validation
    void validatePop2Object(){
        static const string routine =
            "FourierProcessWrapper::validatePop2Object";
        try{
            if (processType.empty()){
                throw ModelException(routine,
                                     "Blank Fourier Process specified!");
            }
            if (CString::equalsIgnoreCase(processType,FOURIERPROCESS_TYPE_SV)) {
                if (processSV.get()) {
                    realProcess = processSV;
                } else {
                    throw ModelException(routine, "Expected Fourier Process SV "
                                         "but none supplied!");
                }
            } else if (CString::equalsIgnoreCase(processType,FOURIERPROCESS_TYPE_SVJ)) {
                if (processSVJ.get()) {
                    realProcess = processSVJ;
                } else {
                    throw ModelException(routine, "Expected Fourier Process SVJ "
                                         "but none supplied!");
                }
            } else if (CString::equalsIgnoreCase(processType,FOURIERPROCESS_TYPE_SVCJ)) {
                if (processSVCJ.get()) {
                    realProcess = processSVCJ;
                } else {
                    throw ModelException(routine, "Expected Fourier Process SVCJ "
                                         "but none supplied!");
                }
            } else if (CString::equalsIgnoreCase(processType,FOURIERPROCESS_TYPE_VSCURVE)) {
                if (processVSCurve.get()) {
                    realProcess = processVSCurve;
                } else {
                    throw ModelException(routine, "Expected Fourier Process VSCURVE "
                                         "but none supplied!");
                }
            } else {
                throw ModelException(routine, "Unrecognised Fourier Process "
                                     + processType + ". Expected "
                                     + FOURIERPROCESS_TYPE_SV + ", "
                                     + FOURIERPROCESS_TYPE_SVJ + ", "
                                     + FOURIERPROCESS_TYPE_SVCJ + " or "
                                     + FOURIERPROCESS_TYPE_VSCURVE);
            }
        }  catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** create a proper ISAP */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const {
        static const string method = "FourierProcessWrapper::convert";
        try {
            if (requiredType != FourierProcess::TYPE) {
                throw ModelException(method,
                                     "Cannot convert a FourierProcessWrapper into "
                                     "object of type "+requiredType->getName());
            }
            object = realProcess;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FourierProcessWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ITypeConvert);
        EMPTY_SHELL_METHOD(defaultFourierProcessWrapper);
        FIELD(processType, "SV, SVJ or SVCJ");
        FIELD(processSV,  "fourier process SV");
        FIELD_MAKE_OPTIONAL(processSV);
        FIELD(processSVJ,  "fourier process SVJ");
        FIELD_MAKE_OPTIONAL(processSVJ);
        FIELD(processSVCJ,  "fourier process SVCJ");
        FIELD_MAKE_OPTIONAL(processSVCJ);
        FIELD(processVSCurve,  "fourier process VSCurve");
        FIELD_MAKE_OPTIONAL(processVSCurve);
        FIELD(realProcess, "real fourier process");
        FIELD_MAKE_TRANSIENT(realProcess);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    // for reflection
    FourierProcessWrapper(): CObject(TYPE){}

    static IObject* defaultFourierProcessWrapper(){
        return new FourierProcessWrapper();
    }
};

typedef smartPtr<FourierProcessWrapper> FourierProcessWrapperSP;

CClassConstSP const FourierProcessWrapper::TYPE =
CClass::registerClassLoadMethod("FourierProcessWrapper",
                                typeid(FourierProcessWrapper), load);


/** "bread & butter" FourierEngine that can be captured in Pyramid using
     current IMS */
class FourierEngineDefault: public FourierEngine {
public:
    static CClassConstSP const TYPE;

private:
    FourierEngineDefault():FourierEngine(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(FourierEngineDefault, clazz);
        SUPERCLASS(FourierEngine);
        EMPTY_SHELL_METHOD(defaultFourierEngineDefault);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultFourierEngineDefault(){
        return new FourierEngineDefault();
    }
};


typedef smartPtr<FourierEngineDefault> FourierEngineDefaultSP;

CClassConstSP const FourierEngineDefault::TYPE =
CClass::registerClassLoadMethod(
    "FourierEngineDefault", typeid(FourierEngineDefault), load);


/////////////////////////////////////////////////////////////////////////


class FourierEngineVSCurve: public FourierEngine {
public:
    static CClassConstSP const TYPE;

private:
    FourierEngineVSCurve():FourierEngine(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(FourierEngineVSCurve, clazz);
        SUPERCLASS(FourierEngine);
        EMPTY_SHELL_METHOD(defaultFourierEngineVSCurve);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultFourierEngineVSCurve(){
        return new FourierEngineVSCurve();
    }
};

typedef smartPtr<FourierEngineVSCurve> FourierEngineVSCurveSP;

CClassConstSP const FourierEngineVSCurve::TYPE = 
CClass::registerClassLoadMethod(
    "FourierEngineVSCurve", typeid(FourierEngineVSCurve), load);


/////////////////////////////////////////////////////////////////////////


class CalibrationEngineSV: public FourierEngine {
public:
    static CClassConstSP const TYPE;

private:
    CalibrationEngineSV():FourierEngine(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(CalibrationEngineSV, clazz);
        SUPERCLASS(FourierEngine);
        EMPTY_SHELL_METHOD(defaultCalibrationEngineSV);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCalibrationEngineSV(){
        return new CalibrationEngineSV();
    }
};

typedef smartPtr<CalibrationEngineSV> CalibrationEngineSVSP;

CClassConstSP const CalibrationEngineSV::TYPE = 
CClass::registerClassLoadMethod(
    "CalibrationEngineSV", typeid(CalibrationEngineSV), load);


/////////////////////////////////////////////////////////////////////////


/* The ultimate wrapping of model, mainly for use in Pyramid
 */
#define MODEL_TYPE_FOURIERENGINE          "FourierEngine"

class ModelWrapper : public CObject,
                     virtual public ITypeConvert {
public: // how can I have this protected or private?
    string              modelType;     // FourierEngine
    FourierEngineSP     fourierEngine;

private:
    IModelSP            realModel;

public:
    static CClassConstSP const TYPE;

    // validation
    void validatePop2Object(){
        static const string routine =
            "ModelWrapper::validatePop2Object";
        try{
            if (modelType.empty()){
                throw ModelException(routine,
                                     "Blank model specified!");
            }
            if (CString::equalsIgnoreCase(modelType,MODEL_TYPE_FOURIERENGINE)) {
                if (fourierEngine.get()) {
                    realModel = fourierEngine;
                } else {
                    throw ModelException(routine, "Expected Fourier Engine "
                                         "but none supplied!");
                }
            } else {
                throw ModelException(routine, "Unrecognised model "
                                     + modelType + ". Expected "
                                     + MODEL_TYPE_FOURIERENGINE);
            }
        }  catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** create a proper model */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const {
        static const string method = "ModelWrapper::convert";
        try {
            if (requiredType != IModel::TYPE) {
                throw ModelException(method,
                                     "Cannot convert a ModelWrapper into "
                                     "object of type "+requiredType->getName());
            }
            object = realModel;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ModelWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ITypeConvert);
        EMPTY_SHELL_METHOD(defaultModelWrapper);
        FIELD(modelType, "FourierEngine");
        FIELD(fourierEngine,  "fourier engine");
        FIELD_MAKE_OPTIONAL(fourierEngine);
        FIELD(realModel, "real model");
        FIELD_MAKE_TRANSIENT(realModel);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    // for reflection
    ModelWrapper(): CObject(TYPE){}

    static IObject* defaultModelWrapper(){
        return new ModelWrapper();
    }
};

typedef smartPtr<ModelWrapper> ModelWrapperSP;

CClassConstSP const ModelWrapper::TYPE =
CClass::registerClassLoadMethod("ModelWrapper",
                                typeid(ModelWrapper), load);

class EmptyISAPHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(EmptyISAP, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FourierEngine::ISAP);
        EMPTY_SHELL_METHOD(defaultEmptyISAP);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultEmptyISAP(){
        return new EmptyISAP();
    }
};

CClassConstSP const EmptyISAP::TYPE = CClass::registerClassLoadMethod(
    "EmptyISAP", typeid(EmptyISAP), EmptyISAPHelper::load);

DRLIB_END_NAMESPACE
