/**
 * @file VanillaGridMultiInstrumentCollection.cpp
 *
 * See VanillaGridMultiInstrumentCollection.hpp for comments
 */

#include "edginc/config.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/VanillaGridMulti.hpp" 
#include "edginc/VanillaGridMultiInstrumentCollection.hpp" 

DRLIB_BEGIN_NAMESPACE

VanillaGridMultiInstrumentCollection::VanillaGridMultiInstrumentCollection():
    CInstrumentCollection(TYPE), grid(0)
{}

VanillaGridMultiInstrumentCollection::VanillaGridMultiInstrumentCollection(
        VanillaGridMultiSP grid):
    CInstrumentCollection(TYPE),
    grid(grid.clone())
{
    validatePop2Object();
}

VanillaGridMultiConstSP VanillaGridMultiInstrumentCollection::getGrid() const {
    return grid;
}

const vector<pair<int, int> >& VanillaGridMultiInstrumentCollection::indices()
        const {
    if (_indices.empty()) {
                
        _indices.resize(grid->avgOutDates.size() * grid->nbPerfPerMat);

        for (int i = 0; i < grid->nbPerfPerMat; ++i) {
            for (int j = 0; j < grid->avgOutDates.size(); ++j) {
                _indices[j + grid->avgOutDates.size() * i] =
                    make_pair(i, j);
            }
        }
    }    

    return _indices;
}

void VanillaGridMultiInstrumentCollection::validatePop2Object() {
    try {
        ASSERT(!!grid);
        indices();
    }
    catch (exception& e) {
        throw ModelException(
            e, "VanillaGridMultiInstrumentCollection::validatePop2Object()");
    }
}                

int VanillaGridMultiInstrumentCollection::size() const {
    return indices().size();
}

InstrumentSP VanillaGridMultiInstrumentCollection::operator [](int v) {
    ASSERT(0 <= v && v < size());    
    int perfIdx = indices()[v].first;
    int dateIdx = indices()[v].second;
        
    DateTimeArray avgOutDatesTmp(1);
    avgOutDatesTmp[0] = grid->avgOutDates[dateIdx];
    IAggregateMakerArraySP basketsTmp(new IAggregateMakerArray(1));
    if ( grid->baskets->size() > 1 ) {
        (*basketsTmp)[0] = (*grid->baskets)[perfIdx];
    } else {
        (*basketsTmp)[0] = (*grid->baskets)[0];
    }
    IDoubleArrayModifierMakerArraySP perfTypesTmp(new IDoubleArrayModifierMakerArray(1));
    if ( grid->perfTypes->size() > 1 ) {
        (*perfTypesTmp)[0] = (*grid->perfTypes)[perfIdx];
    } else {
        (*perfTypesTmp)[0] = (*grid->perfTypes)[0];
    }

    return InstrumentSP(new VanillaGridMulti(grid.get(),
                                             avgOutDatesTmp,                                               
                                             basketsTmp,
                                             perfTypesTmp));
}

void VanillaGridMultiInstrumentCollection::Validate() {
    grid->Validate();
}

void VanillaGridMultiInstrumentCollection::GetMarket(const IModel* model,
                                                const CMarketDataSP market) {
    grid->GetMarket(model, market);
}

DateTime VanillaGridMultiInstrumentCollection::getValueDate() const {
    return grid->getValueDate();
}

string VanillaGridMultiInstrumentCollection::discountYieldCurveName() const {
    return grid->discountYieldCurveName();
}

DateTime VanillaGridMultiInstrumentCollection::endDate(
        const Sensitivity* sensControl) const {    
    return grid->avgOutDates.back();
}

IMCProduct* VanillaGridMultiInstrumentCollection::createProduct(
        const MonteCarlo* model) const {
    return grid->createProduct(model);
}

void VanillaGridMultiInstrumentCollection::Price(IModel *model,
                                            CControl *control,
                                            CResultsArraySP resultss) {
    try {
        ASSERT(resultss->size() == size());

        ResultsSP gridResults(new Results());

        CControlSP gridControl(copy(control));
        gridControl->addRequest(OutputRequestSP(
            new OutputRequest(OutputRequest::OPTION_PRICE)));
        
        model->Price(grid.get(), gridControl.get(), gridResults.get());

        VanillaGrid::OutputConstSP prices = VanillaGrid::OutputConstSP::dynamicCast(
            gridResults->retrieveRequestResult(OutputRequest::OPTION_PRICE));

        VanillaGrid::OutputConstSP stderrs =
            control->requestsOutput(OutputRequest::VALUE_STE) ?
                VanillaGrid::OutputConstSP::dynamicCast(
                    gridResults->retrieveRequestResult(OutputRequest::VALUE_STE)) :
                VanillaGrid::OutputConstSP(   );

        for (int v = 0; v < size(); ++v) {
            int i = indices()[v].first;     // perfIdx
            int j = indices()[v].second;    // dateIdx

            ResultsSP results = (*resultss)[v];

            results->storePrice(prices->getValue(i, j),
                                grid->discountYieldCurveName());

            if (!!stderrs && control->isPricing()) {                

                OutputRequest* request = 
                    control->requestsOutput(OutputRequest::VALUE_STE);
                if (request) {
                    results->storeRequestResult(request,
                                                stderrs->getValue(i, j));
                }
            }
        }            
    }
    catch (exception& e) {
        throw ModelException(e, "VanillaGridMultiInstrumentCollection::Price()");
    }
}

void VanillaGridMultiInstrumentCollection::scaleOutputs(
    CControlSP control, CResultsArraySP results)
{}

IObject* VanillaGridMultiInstrumentCollection::
        defaultVanillaGridMultiInstrumentCollection() {
    return new VanillaGridMultiInstrumentCollection();
}

void VanillaGridMultiInstrumentCollection::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VanillaGridMultiInstrumentCollection, clazz);
    SUPERCLASS(CInstrumentCollection);
    EMPTY_SHELL_METHOD(defaultVanillaGridMultiInstrumentCollection);
    IMPLEMENTS(IMCIntoProduct);
    FIELD(grid, "grid");
}

CClassConstSP VanillaGridMultiInstrumentCollection::TYPE =
    CClass::registerClassLoadMethod(
        "VanillaGridMultiInstrumentCollection", typeid(VanillaGridMultiInstrumentCollection),
        VanillaGridMultiInstrumentCollection::load);

bool VanillaGridMultiInstrumentCollectionLoad() {
    return VanillaGridMultiInstrumentCollection::TYPE != NULL;
}

DRLIB_END_NAMESPACE
