/**
 * @file VanillaGridInstrumentCollection.cpp
 *
 * See VanillaGridInstrumentCollection.hpp for comments
 */

#include "edginc/config.hpp"
#include "edginc/VanillaGrid.hpp" 
#include "edginc/VanillaGridInstrumentCollection.hpp" 

DRLIB_BEGIN_NAMESPACE

VanillaGridInstrumentCollection::VanillaGridInstrumentCollection():
    CInstrumentCollection(TYPE), grid(0)
{}

VanillaGridInstrumentCollection::VanillaGridInstrumentCollection(
        VanillaGridSP grid):
    CInstrumentCollection(TYPE),
    grid(grid.clone())
{
    validatePop2Object();
}

VanillaGridConstSP VanillaGridInstrumentCollection::getGrid() const {
    return grid;
}

const vector<pair<int, int> >& VanillaGridInstrumentCollection::indices()
        const {
    if (_indices.empty()) {
        const DoubleMatrix* w = grid->getWeights().get();

        if (w) {
            ASSERT(w->numRows() == grid->strikes->size());
            ASSERT(w->numCols() == grid->maturities->size());

            for (int i = 0; i < w->numRows(); ++i) {
                for (int j = 0; j < w->numCols(); ++j) {
                    if (Maths::isPositive((*w)[j][i])) {
                        _indices.push_back(make_pair(i, j));
                    }
                }
            }
        }
        else {
            _indices.resize(grid->maturities->size() * grid->strikes->size());

            for (int i = 0; i < grid->strikes->size(); ++i) {
                for (int j = 0; j < grid->maturities->size(); ++j) {
                    _indices[j + grid->maturities->size() * i] =
                        make_pair(i, j);
                }
            }
        }
    }

    return _indices;
}

void VanillaGridInstrumentCollection::validatePop2Object() {
    try {
        ASSERT(!!grid);
        indices();
    }
    catch (exception& e) {
        throw ModelException(
            e, "VanillaGridInstrumentCollection::validatePop2Object()");
    }
}                

int VanillaGridInstrumentCollection::size() const {
    return indices().size();
}

DoubleArrayConstSP VanillaGridInstrumentCollection::weights() const {
    if (!_weights) {
        const DoubleMatrix* w = grid->getWeights().get();

        if (!w) {
            _weights.reset(new DoubleArray(size(), 1.));
        }
        else {
            _weights.reset(new DoubleArray(size()));
            DoubleArray& _w = const_cast<DoubleArray&>(*_weights);

            for (int p = 0; p < _w.size(); ++p) {
                const pair<int, int>& ij = _indices[p];
                _w[p] = (*w)[ij.second][ij.first];
            }
        }
    }

    return _weights;
}

InstrumentSP VanillaGridInstrumentCollection::operator [](int v) {
    ASSERT(0 <= v && v < size());

    ScheduleSP schedule(new Schedule(
        DateTimeArray(1, (*grid->maturities)[indices()[v].second]),
        DoubleArray(1, (*grid->instStrikesUsed)[indices()[v].second][indices()[v].first]),
        string("N")));

    return InstrumentSP(CVanilla::make(grid->valueDate,
                                       (*(*grid->instTypesUsed)[indices()[v].second])[indices()[v].first],
                                       false,
                                       grid->oneContract,
                                       grid->notional,
                                       grid->initialSpot,
                                       schedule.get(),
                                       grid->asset,
                                       grid->discount,
                                       grid->instSettle.get(),
                                       false));
}

void VanillaGridInstrumentCollection::Validate() {
    grid->Validate();
}

void VanillaGridInstrumentCollection::GetMarket(const IModel* model,
                                                const CMarketDataSP market) {
    grid->GetMarket(model, market);
}

DateTime VanillaGridInstrumentCollection::getValueDate() const {
    return grid->getValueDate();
}

string VanillaGridInstrumentCollection::discountYieldCurveName() const {
    return grid->discountYieldCurveName();
}

DateTime VanillaGridInstrumentCollection::endDate(
        const Sensitivity* sensControl) const {
    return grid->endDate(sensControl);
}

IMCProduct* VanillaGridInstrumentCollection::createProduct(
        const MonteCarlo* model) const {
    return grid->createProduct(model);
}

void VanillaGridInstrumentCollection::Price(IModel *model,
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
            int i = indices()[v].first;
            int j = indices()[v].second;

            ASSERT(!prices->isError(j, i));

            ResultsSP results = (*resultss)[v];

            results->storePrice(prices->getValue(j, i),
                                grid->discountYieldCurveName());

            if (!!stderrs && control->isPricing()) {
                ASSERT(!stderrs->isError(j, i));

                OutputRequest* request = 
                    control->requestsOutput(OutputRequest::VALUE_STE);
                if (request) {
                    results->storeRequestResult(request,
                                                stderrs->getValue(j, i));
                }
            }
        }            
    }
    catch (exception& e) {
        throw ModelException(e, "VanillaGridInstrumentCollection::Price()");
    }
}

void VanillaGridInstrumentCollection::scaleOutputs(
    CControlSP control, CResultsArraySP results)
{}

IObject* VanillaGridInstrumentCollection::
        defaultVanillaGridInstrumentCollection() {
    return new VanillaGridInstrumentCollection();
}

void VanillaGridInstrumentCollection::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VanillaGridInstrumentCollection, clazz);
    SUPERCLASS(CInstrumentCollection);
    EMPTY_SHELL_METHOD(defaultVanillaGridInstrumentCollection);
    IMPLEMENTS(IMCIntoProduct);
    FIELD(grid, "grid");
}

CClassConstSP VanillaGridInstrumentCollection::TYPE =
    CClass::registerClassLoadMethod(
        "VanillaGridInstrumentCollection", typeid(VanillaGridInstrumentCollection),
        VanillaGridInstrumentCollection::load);

bool VanillaGridInstrumentCollectionLoad() {
    return VanillaGridInstrumentCollection::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
