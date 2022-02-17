//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PutCliquet.cpp
//
//   Description : Put Cliquet instrument
//
//   Author      : Louis A Moussu
//
//   Date        : 13 Decembre 2005
//
//   
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/Black.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/XCB.hpp"
#include "edginc/RootFinder.hpp"
#include "edginc/FunctionWrapper.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/EndDateCollector.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/VolSVJ.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

// Model for closed form
class PutCliquetCFModel : public CModel {

    friend class PutCliquet;
    friend class PutCliquetCFPricer;

public:
        
    static CClassConstSP const TYPE;

private:

    bool allowRiskMapping;

public:

    virtual void Price(CInstrument* instrument,
                       CControl*	control,
                       CResults*	results){
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        auto_ptr<IProduct> product(intoProd.createProduct(this));
        product->price(this,control,results);
    };

    PutCliquetCFModel():
        CModel(TYPE),
        allowRiskMapping(false)
    {}

    virtual MarketDataFetcherSP createMDF() const {
        MDFAssetVolSP mdf(new  MDFAssetVol("VolSVJ"));
        return mdf;
    }
		
    IModel::WantsRiskMapping wantsRiskMapping() const {
        return allowRiskMapping ? riskMappingAllowed : riskMappingDisallowed;
    }

	// Product
    class IProduct{

    public:

        virtual void price(PutCliquetCFModel* model,
						   Control* control,
						   Results* results) const = 0;
        virtual ~IProduct() {};
    };

    // Product Interface
	class IIntoProduct : virtual public CModel::IModelIntoProduct{

		public:
			static CClassConstSP const TYPE;
			virtual IProduct* createProduct(const PutCliquetCFModel* model) const = 0;
    };

	static void IntoProduct_load(CClassSP& clazz){
		clazz->setPublic();
		REGISTER_INTERFACE(Model::IModelIntoProduct, clazz);
		EXTENDS(Model::IModelIntoProduct);
	}

	static void load(CClassSP& clazz){
		clazz->setPublic();
		REGISTER(PutCliquetCFModel, clazz);
		SUPERCLASS(CModel);
		EMPTY_SHELL_METHOD(defaultPutCliquetCFModel);
        FIELD(allowRiskMapping,
                     "Support B-S greeks if a suitable RiskMappingMatrix "
                         "is available");
	}

	static void loadIntoProduct(CClassSP& clazz){
		clazz->setPublic();
		REGISTER_INTERFACE(PutCliquetCFModel::IIntoProduct, clazz);
		EXTENDS(Model::IModelIntoProduct);
	}

	static IObject* defaultPutCliquetCFModel(){
		return new PutCliquetCFModel();
	}
};

typedef smartPtr<PutCliquetCFModel> PutCliquetCFModelSP;

CClassConstSP const PutCliquetCFModel::TYPE = CClass::registerClassLoadMethod(
    "PutCliquetCFModel", typeid(PutCliquetCFModel), PutCliquetCFModel::load);


// Put Cliquet
class PutCliquet: public Generic1Factor,
                  public IMCIntoProduct,
				  public PutCliquetCFModel::IIntoProduct {

    public:

        static CClassConstSP const TYPE;

        virtual void validatePop2Object();

        virtual void Validate();
        virtual DateTime getValueDate() const;
        virtual string discountYieldCurveName() const;
        
        virtual IMCProduct* createProduct(const MonteCarlo* model) const;
		virtual PutCliquetCFModel::IProduct* createProduct(const PutCliquetCFModel* model) const;

    private:

        friend class PutCliquetHelper;
        friend class PutCliquetMC;
		friend class PutCliquetCFPricer;

        PutCliquet(): Generic1Factor(TYPE) {} 

        double                      strike;
        DateTimeArray               maturities;
        DateTimeArray               qMaturities;
        double                      coupon;

    protected:

        IRefLevelSP             refLevel; // $unregistered
        IPastValuesSP           pastValues;   // $unregistered
};





typedef array<smartPtr<PutCliquet>, PutCliquet> PutCliquetArray;

DEFINE_TEMPLATE_TYPE(PutCliquetArray);

// Monte Carlo product
class PutCliquetMC : public MCProductClient {

    private:
    
        const PutCliquet*     inst;
        bool                  active;     

        // state vars
        SVGenSpotSP                  spotGen;      //!< Generator for spot
        SVGenSpot::IStateVarSP       spotSV;       //!< Spot state variable
        SVGenDiscFactorSP            dfGen;        //!< Generator for discount factors
        SVDiscFactorSP		  dfSV;         //!< Df state variable
        SVGenDiscFactorSP            dfGenQ;        //!< Generator for discount factors
        SVDiscFactorSP		  dfSVQ;         //!< Df state variable

        // Override default method on IMCProduct. This method is called every time
        // the path generator is changed (which is, at the moment, when the
        // past path generator is created, and then when the future path
        // generator is created
        virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
            static const string routine = "VanillaMCSV::pathGenUpdated";

            try {
                spotSV = spotGen->getSpotSV(newPathGen);
                dfSV = dfGen->getSVDiscFactor(newPathGen);
                dfSVQ = dfGenQ->getSVDiscFactor(newPathGen);
            } catch (exception& e) {
                throw ModelException(e, routine);
            }
        };

    public:    
        // Appends 'true' (ie non derived) state variable generators
        // required to the supplied collector.
        virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
            // ask for a reference level State Variable
            svCollector->append(spotGen.get());
            svCollector->append(dfGen.get());
            svCollector->append(dfGenQ.get());
        }

        // Use this opportunity to do any LogNormal driven initialisation
        // of the instrument before the main MC loop. e.g closed form barrier adjustment */
        void initialiseLN(const  IMCPathGenerator*  pathGen)const{
            // empty
        }

        PutCliquetMC(const PutCliquet*       inst,
                     const SimSeriesSP       simSeries,
                     const IRefLevelSP       refLevel,
                     const IPastValuesSP     pastValues):
        MCProductClient(inst->asset.get(),
                        inst->valueDate,
                        inst->discount.get(),
                        refLevel,
                        simSeries,
                        pastValues,
                        inst->instSettle.get(),
                        simSeries->getLastDate()),
        inst(inst),
        spotGen(new SVGenSpot(simSeries)),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->instSettle, inst->maturities)),
        dfGenQ(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->instSettle, inst->qMaturities)),
        active(true){}

        void payoff(const IPathGenerator*  pathGen,
                    IMCPrices&                prices) {
            
            double varLeg;
            double fixLeg;
            if (!active){
                varLeg = 0;
                fixLeg = 0;
            } else {
                int iAsset = 0;
                const SVPath& path = spotSV->path(iAsset);
                const SVPath& dfPath = dfSV->path();
                const SVPath& dfPathQ = dfSVQ->path();
            
                // Variable Leg
                bool stop = false;
                int iStart, iStop;
                if (path.begin()==0){
                    iStart = 1;
                } else {
                    iStart = path.begin();
                }
                int iStep = iStart;
                while ((!stop)&&(iStep<path.end())){
                    if (path[iStep]/path[iStep-1]<inst->strike){
                        stop = true;
                        iStop = iStep;
                    }
                    iStep += 1;
                }
                
                if (!stop) {
                    varLeg = 0;
                } else {
                    varLeg = inst->strike - path[iStop]/path[iStop-1];
                    varLeg *= dfPath[iStop];
                }

                // Fix Leg
                int iStartQ; 
                if (dfPathQ.begin()==0){
                    iStartQ = 1;
                } else {
                    iStartQ = dfPathQ.begin();
                }
                if (!stop){
                    fixLeg = 0;
                    for (iStep = iStartQ; iStep<dfPathQ.end(); iStep++){
                        fixLeg += inst->coupon*inst->qMaturities[iStep-1].yearFrac(inst->qMaturities[iStep])
                            *dfPathQ[iStep];
                    }
                } else {                    
                    fixLeg = 0;
                    iStep = iStartQ;
                    while ((inst->qMaturities[iStep].isLess(inst->maturities[iStop]))
                                &&(iStep<inst->qMaturities.size())){
                                fixLeg += inst->coupon*inst->qMaturities[iStep].yearFrac(inst->qMaturities[iStep+1])
                                *dfPathQ[iStep+1];
                                iStep++;
                            }
                    fixLeg += inst->coupon*inst->qMaturities[iStep-1].yearFrac(inst->maturities[iStop])
                        *dfPath[iStop];
          
                }

                if ((pathGen->doingPast())
                    &&(stop)){
                    active = false;
                }
            }
            prices.add(varLeg-fixLeg);
        }    
};




// Registration of PutCliquet
class PutCliquetHelper {

    public:
    // Invoked when Class is 'loaded'
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(PutCliquet, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultPutCliquet);
        FIELD(strike, "Strike");
        FIELD(maturities, "Maturities");
        FIELD(qMaturities, "Quarterly Maturities");
        FIELD(coupon, "Coupon");
    }

    static IObject* defaultPutCliquet(){
        return new PutCliquet();
    }
};

// Put Cliquet Product
class PutCliquetCFPricer : public PutCliquetCFModel::IProduct{

    private:

        const PutCliquet*     inst;
		
	public:

        PutCliquetCFPricer(const PutCliquet*     inst):
            IProduct(),
            inst(inst){};

		virtual void price(PutCliquetCFModel* model,
						   Control* control,
						   Results* results) const {

			// Get VolSVJ
			VolSVJSP vol = VolSVJSP::dynamicCast(VolRequestRaw::copyVolBase(*(inst->asset).get()));	

			double price = vol->calcPricePutCliquet(inst->strike,
								   inst->coupon,
								   inst->valueDate,
								   inst->maturities,
								   inst->qMaturities,
								   inst->discount,
								   inst->asset);
			results->storePrice(price,inst->discount->getCcy());
		};        
};

void PutCliquet::Validate()
{

}

void PutCliquet::validatePop2Object()
{

}


DateTime PutCliquet::getValueDate() const
{
  return valueDate;
}

string PutCliquet::discountYieldCurveName() const {
    return discount.getName();
}


// Implementation of MonteCarlo::IIntoProduct interface
IMCProduct* PutCliquet::createProduct(const MonteCarlo* model) const {

    SimSeriesSP simSeries(new SimSeries(1));
    simSeries->addDates(maturities);
    // to be changed !
    IRefLevelSP refLevel(IRefLevel::Util::makeTrivialAverage(startDate));
    IPastValuesSP pastValues(IPastValues::Util::makeTrivial(startDate, 1));


    return new PutCliquetMC(this, simSeries, refLevel, pastValues);
}

// Implementation of PutCliquetCFModel::IIntoProduct interface
PutCliquetCFModel::IProduct* PutCliquet::createProduct(const PutCliquetCFModel* model) const{
	return new PutCliquetCFPricer(this);
}


CClassConstSP const PutCliquet::TYPE = CClass::registerClassLoadMethod(
                                        "PutCliquet", typeid(PutCliquet), PutCliquetHelper::load);

class ObjFuncPutCliquet: public Calibrator::ObjFuncLeastSquare{

	public:

		friend class ObjFuncPutCliquetHelper;

		static CClassConstSP const TYPE;

		ObjFuncPutCliquet():Calibrator::ObjFuncLeastSquare(TYPE),
										model(PutCliquetCFModelSP(new PutCliquetCFModel)){};

		virtual void getMarket(MarketData* market){
			MarketDataSP marketSP(market);
			for (int i=0; i < putCliquets.size(); i++){
                putCliquets[i]->GetMarket(model.get(), marketSP);              
			}
		}

		virtual IObjectSP getAdjustableGroup(){
			return IObjectSP::attachToRef(this);
		}
		
		virtual int getNbFuncs() const{
			return putCliquets.size();
		}

		virtual void calcValue(CDoubleArray&   funcvals) const{
			CControlSP control(new Control());
			for (int i=0; i<putCliquets.size(); i++){
				PutCliquetCFPricer putCliquetCFPricer(putCliquets[i].get());
				ResultsSP results(new Results());
				putCliquetCFPricer.price(model.get(), control.get(), results.get());
				funcvals[i] = results->retrievePrice();
			}
		}

	private:

		PutCliquetArray			putCliquets;
		PutCliquetCFModelSP		model; // $unregistered
};

class ObjFuncPutCliquetHelper {

	public:

		static void load(CClassSP& clazz){
			clazz->setPublic(); 
			REGISTER(ObjFuncPutCliquet, clazz);
			SUPERCLASS(Calibrator::ObjFuncLeastSquare);
			EMPTY_SHELL_METHOD(defaultObjFuncPutCliquet);
			FIELD(putCliquets, "Put Cliquets");
		}

		static IObject* defaultObjFuncPutCliquet(){
			return new ObjFuncPutCliquet();
		}
};

CClassConstSP const ObjFuncPutCliquet::TYPE = CClass::registerClassLoadMethod(
    "ObjFuncPutCliquet", typeid(ObjFuncPutCliquet), ObjFuncPutCliquetHelper::load);

bool  PutCliquetLoad() {
    return (PutCliquet::TYPE != 0);
}

DRLIB_END_NAMESPACE
