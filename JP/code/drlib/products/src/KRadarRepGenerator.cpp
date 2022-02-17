#include "edginc/config.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"

#include "edginc/TreeSlice.hpp"
#include "edginc/TreeSliceExpr.hpp"
#include "edginc/TreeSliceOper.hpp"

#include "edginc/RadarRepCreator.hpp"
//#include "edginc/IRadarRepCreator.hpp"
#include "edginc/IRadarRep.hpp"
#include "edginc/PolynomialBasis.hpp"
#include "edginc/Results.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/RadarRepWrapper.hpp"
#include "edginc/RadarRepDealWrapper.hpp"
#include "edginc/KFloatLeg.hpp"
#include "edginc/KSum.hpp"
#include "edginc/KOption.hpp"
#include "edginc/KRadarRepGenerator.hpp"
#include <fstream>
#include <vector>

DRLIB_BEGIN_NAMESPACE

using namespace std;
const char* OBS_SLICE_FILE = "C:\\temp\\obs.txt";
const char* FITTING_SLICE_FILE = "C:\\temp\\fitting.txt";
const char* PREDICTED_SLICE_FILE = "C:\\temp\\predicted.txt";

/*Kranthi Gade:
RadarRep generation: Given an array/vector of observation dates, an observable
FDProduct, a vector of fitting variable FDProducts, a collection of basis functions, this function
returns a set of regression coefficients each of which best represents the value of the observable product 
at a particular time step best in terms of the values of the fitting variable products at that
particular time step */

 class  RadarRepMarker : public SliceMarker< RadarRepMarker >  {
    public:
        // Gets called on each tp
        // TODO:  Need to collect state variables from product and take those 
        // into account.  These are used for path dependency

        RadarRepMarker(
            const TreeSlice* obsSlice,
            vector<const TreeSlice* >  fittingSliceVec,
            IRadarRepCreatorSP pc,
            vector<FittingArray>& fArrSlice,
            vector<double>& obsValSlice)
            :
            sliceCount(1 + fittingSliceVec.size()),
            obsSlice(obsSlice),
            fittingSliceVec(fittingSliceVec),
            pc(pc),
            fArrSlice(fArrSlice),
            obsValSlice(obsValSlice)
         {
            TreeSlice::loopOnSlices(*this);
         }

        // TreeSlice "expression template" primitives
        // TODO:  Take state variables into account for path dependency
        // Add all slices to list (S**)
        /*template< typename S >
            const S** listSlices(const S** list) const
        {
            int n = fittingSliceVec.size();

            for( int i = 0; i < n; ++i )
                list = fittingSliceVec[ i ]->listSlices( list );

            return obsSlice->listSlices( list );
        }*/
        template< typename S >
        const S** listInputSlices(const S** list) const
        {
            for( int i = 0; i < fittingSliceVec.size(); ++i )
                 list = fittingSliceVec[ i ]->listSlices( list );

            *list++ = obsSlice;

             return list;
        }

        template< typename S >
        S** listOutputSlices(S** list) const
        {
            //*list = obsSlice;
            return list;
        }

        //inline double calc() const
        void compute() const
        {
            int n = fittingSliceVec.size();

            ASSERT( n >= 1 );
            FittingArray fArr;
            for (int i=0; i<n; i++) {
                fArr.push_back(fittingSliceVec[i]->calc());
            }
            double obsVal = obsSlice->calc();
            pc->addSample(fArr, obsVal);

            /*For debugging purposes -- dump the slices to files just to see them */
            fArrSlice.push_back(fArr);
            obsValSlice.push_back(obsVal);

            ofstream fittingOut;
            fittingOut.open(FITTING_SLICE_FILE, ios_base::app);
            for( int i = 0; i < n; ++i ) {
                fittingOut << fArr[i] << " ";
            }
            fittingOut << endl;
            fittingOut.close();
            
            ofstream obsOut;
            obsOut.open(OBS_SLICE_FILE, ios_base::app);
            obsOut << obsVal;
            obsOut << endl;
            obsOut.close();
            
            /*End debugging */
        }
        void printDebug(char *s) const {
            strcat(s, "(RadarRepMarker)");
        }

        const int sliceCount;
    private:
        const TreeSlice* obsSlice;
        vector<const TreeSlice*>  fittingSliceVec;
        IRadarRepCreatorSP pc;

        mutable vector<FittingArray>& fArrSlice;
        mutable vector<double>& obsValSlice;
        //vector<vector<vector<double> > >& output;
};
//int RadarRepMarker::sliceCount = 1;

class  KRadarRepGeneratorTree : public FDProduct {
public:
    KRadarRepGeneratorConstSP  inst;
    FDProductArray  fittingProdVec;
    FDProductSP obsProd;

    /******************** methods ********************/
    KRadarRepGeneratorTree(const KRadarRepGeneratorConstSP &inst, FDModel* model);
    virtual DateTime getStartDate(void) const { return model->getDate(0);}

    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
        recordSliceToOutputName(ctrl, results, model, 
            inst->isTopComponent(), inst->outputName, "", getValue(0));
        inst->recordExtraOutput(ctrl, results);
        //results->storeResult("radar",OutputNameConstSP(new OutputName("RadarRep")), iob);
        //KRadarRepGeneratorSP kRadarRep;
        /*IRadarRepSP iRadarRep = pc->getRadarRep();
        //RadarRepSP myRadarRep = RadarRepSP(dynamic_cast<RadarRep*>(iRadarRep.get()));
        RadarRepSP myRadarRep = RadarRepSP::dynamicCast(iRadarRep);*/
        
        OutputRequest* request = ctrl->requestsOutput(OutputRequest::RADAR_REP);
        if (request == NULL) {
            throw ModelException(OutputRequest::RADAR_REP+" request has not been added to Control object");
        }
        /*RegressionCoeff coefVec = myRadarRep->getCoef();
        DoubleArraySP arr = DoubleArraySP(new DoubleArray(coefVec.begin(), coefVec.end()));
        results->storeRequestResult(request, RadarRepWrapperSP(new RadarRepWrapper(arr, inst->funcBasis, inst->fittingTransform)));*/
        RadarRepDealSP deal = RadarRepDealSP(new RadarRepDeal(radarMap));
        results->storeRequestResult(request, RadarRepDealWrapperSP(new RadarRepDealWrapper(deal)));
    }

    virtual void init(Control*) const {}
    virtual void initProd(void);

    virtual void update(int& step, UpdateType type) {
        DateTime stepDate = model->getDate(step);
        if ( stepDate.findUpper( inst->obsDates ) != stepDate.findLower(inst->obsDates) ) 
            return;   /*Return if date corresponding to 'step' is not one of the observation dates.*/

        const TreeSlice* obsSlice = &(obsProd->getValue(step, stepDate));
        vector<const TreeSlice*> fittingSliceVec;
        for (size_t i=0; i<inst->fittingKVec.size(); i++) {
            fittingSliceVec.push_back(&(fittingProdVec[i]->getValue(step, stepDate)));
        }
        
        ofstream obsOut;
        obsOut.open(OBS_SLICE_FILE, ios_base::app);
        obsOut << stepDate.toString() << endl;
        
        ofstream fittingOut;
        fittingOut.open(FITTING_SLICE_FILE, ios_base::app);
        fittingOut << stepDate.toString() << endl;

        pc->clearSamples();
        fArrSlice.clear();
        obsValSlice.clear();
        RadarRepMarker(obsSlice, fittingSliceVec, pc, fArrSlice, obsValSlice);

        IRadarRepSP myRadarRep = pc->getRadarRep();
        radarMap[stepDate] = myRadarRep;
        
        ofstream predicted; 
        predicted.open(PREDICTED_SLICE_FILE, ios_base::app);
        predicted << stepDate.toString() << endl;
        for( int i = 0; i < fArrSlice.size(); ++i ) {
            predicted << (*myRadarRep)(fArrSlice[i]) << endl;
        }
        predicted.close();
    } 
    virtual const TreeSlice & getValue(int step) const;
    virtual const TreeSlice & getValue(int step, DateTime stepDate) const;

private:
    TreeSliceSP slice;
    IRadarRepCreatorSP pc;
    map<DateTime, IRadarRepSP> radarMap;

    vector<FittingArray> fArrSlice;
    vector<double> obsValSlice;
};

/****************************** KRadarRepGeneratorTree ********************************/
KRadarRepGeneratorTree::KRadarRepGeneratorTree(const KRadarRepGeneratorConstSP &inst, FDModel* model) :
FDProduct(model), inst(inst) {    
    int i;
    try {
        int numFitting=inst->fittingKVec.size();
        //fittingProdVec.resize(s);
        for (i=0; i<numFitting; ++i) {
            fittingProdVec.push_back(model->createProduct(inst->fittingKVec[i]));
            fittingProdVec.back()->addModelResetDates(inst->obsDates);
        }
        obsProd = model->createProduct(inst->obsK);
        pc = IRadarRepCreatorSP(new RadarRepCreator(inst->funcBasis->getBasis(), inst->fittingTransform->getTransform()));        
    }
    catch (exception& e){
        string errMsg = "Error creating dependent product of KRadarRepGenerator (fittingKVec array index = " +
            Format::toString(i) + "), object outputName = " + inst->outputName;
        throw ModelException(e, "KRadarRepGeneratorTree::KRadarRepGeneratorTree", errMsg);
    }
}

/** initializing and setting product variables */
void KRadarRepGeneratorTree::initProd(void){
}
const TreeSlice & KRadarRepGeneratorTree::getValue(int step) const {
      /*Not important.  getValue() for a top-level component gets
        called at step 0 only.  The results from this method are not
        useful in any fashion. */
    try {        
        DateTime stepDate = model->getDate(step);
        return getValue(step, stepDate);
    }
    catch (exception& e) {
        throw ModelException(e, "KRadarRepGeneratorTree::getValue, timeStep = " + Format::toString(step) +
                             + " Date= " + model->getDate(step).toString() +
                             ", outputName = " + getOutputName());
    }  
}

const TreeSlice & KRadarRepGeneratorTree::getValue(int step, DateTime stepDate) const {
    /*Not important.  getValue() for a top-level component gets called at step 0 only.  The results
    from this method are not useful in any fashion. */
    try {        
        const TreeSlice& slice = obsProd->getValue(step, stepDate);
        return slice;
    }
    catch (exception& e) {
        throw ModelException(e, "KRadarRepGeneratorTree::getValue, timeStep = " + Format::toString(step) +
            ", outputName = " + getOutputName());
    }
}

/*********************************** KRadarRepGenerator ***********************************/

void KRadarRepGenerator::add(IProdCreatorSP el) {
    fittingKVec.push_back(el);
}

/*void KRadarRepGenerator::addValueDates(const DateTimeArray &dates) {
    valueDates.insert(valueDates.end(), dates.begin(), dates.end());
}*/

void KRadarRepGenerator::validatePop2Object(void) {
}

void KRadarRepGenerator::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        validatePop2Object();

        string disc = discountYieldCurveName();

        vector<DateTimeArray> valueDatesVec;
        for (int i=0; i<fittingKVec.size(); ++i) {

            fittingKVec[i]->setup(model, market);
            valueDatesVec.push_back(getResetDates(fittingKVec[i]));
        }
        obsK->setup(model,market);
        valueDatesVec.push_back(getResetDates(obsK));
        DateTimeArray valueDates = DateTime::merge(valueDatesVec);
        addDatesIndexSpec(fittingKVec, valueDates);
    }
    catch (exception& e) {
        throw ModelException(e, "KRadarRepGenerator::setup "+outputName);
    }
}

double KRadarRepGenerator::getValue(DateTime date, bool &canDo, CashflowInfo &cfi) const {
    return 0;
}


FDProductSP KRadarRepGenerator::createProduct(FDModel * model) const {
    if (!setupCalled) 
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KRadarRepGeneratorTree(KRadarRepGeneratorConstSP(this), model));
}

DateTimeArray KRadarRepGenerator::getResetDates(IProdCreatorSP prod) {
    KComponent* kComp = dynamic_cast<KComponent *>(prod.get());
    if (kComp == NULL) {
        return DateTimeArray();
    }
    else {
        KFloatLeg* kFloat = dynamic_cast<KFloatLeg*>(kComp);
        if (kFloat != NULL) {
            return kFloat->sched->resetEff.getDates();
        }
        else {
            KSum* kSum = dynamic_cast<KSum*>(kComp);
            if (kSum != NULL) {
                vector<DateTimeArray> valueDates;
                CModel::IProdCreatorArray listK = kSum->getListK();
                for (size_t i=0; i<listK.size(); i++) {
                    valueDates.push_back(getResetDates(listK[i]));
                }
                return DateTime::merge(valueDates);
            }
            else {
                KOption* kOption = dynamic_cast<KOption*>(kComp);
                if (kOption != NULL) {
                    return getResetDates(kOption->und);
                }
                else {
                    return DateTimeArray();
                }
            }
        }
    }
}

void KRadarRepGenerator::addDatesIndexSpec(CModel::IProdCreatorArray fittingKVec, DateTimeArray resetDates) {
    for (size_t i=0; i<fittingKVec.size(); i++) {
        if (dynamic_cast<IndexSpec *>(fittingKVec[i].get()) != NULL) {
            IndexSpecSP index = IndexSpecSP::dynamicCast(fittingKVec[i]);
            index->addResetDates(resetDates);
        }
    }
}

void KRadarRepGenerator::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("A radar generator for trees");
    REGISTER(KRadarRepGenerator, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(fittingKVec, "Fitting variable product vector");
    
    FIELD(obsDates, "Observation dates for the radar generator");
    //FIELD_MAKE_OPTIONAL(obsDates);
    //FIELD_INLINE(obsKVec, "Observable variable product vector");
    FIELD(obsK, "Observable variable product");
    FIELD(funcBasis, "Which basis function to be used");
    FIELD(fittingTransform, "How to transform fitting variables");
    Addin::registerConstructor(Addin::UTILITIES, KRadarRepGenerator::TYPE);
}

CClassConstSP const KRadarRepGenerator::TYPE = CClass::registerClassLoadMethod(
    "KRadarRepGenerator", typeid(KRadarRepGenerator), KRadarRepGenerator::load);

/******************************/
// for type linking
bool KRadarRepGeneratorLoad(void){
    return (KRadarRepGenerator::TYPE != 0);
}

DRLIB_END_NAMESPACE
