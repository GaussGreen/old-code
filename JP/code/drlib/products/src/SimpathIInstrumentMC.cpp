//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SimPathIInstrumentMC.cpp
//
//   Description : MonteCarlo impl. of a Pseudo-instrument for RM applications
//
//   Author      : Afshin Bayrooti
//
//   Date        : October 2005
//
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"

#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/SVGenAggregatedSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedEnergyFuture.hpp"
#include "edginc/SVGenExpectedBasisForward.hpp"

#include "edginc/SimpathIInstrumentMC.hpp"

DRLIB_BEGIN_NAMESPACE

SimpathIInstrumentMC::SimpathIInstrumentMC(
        const IMultiMarketFactors* mFactors,
        const DateTime& today,
        const YieldCurve* discount,
        const SimPathICollectionMap& _collectionMap,
        const DateTimeArray& _timeLine,
        IMCWriterSP _mcWriter,
        const map<string, int>& _asset2rng,
	    bool _doPrices) :
    MCProductClient(mFactors, today, discount),
    collectionMap(_collectionMap),
    timeLine(_timeLine),
    mcWriter(_mcWriter),
    asset2rng(_asset2rng),
    doPrices(_doPrices)
{
    mcWriter->notifyStartHeader();
}

void SimpathIInstrumentMC::pathGenUpdated(
        IStateVariableGen::IStateGen* newPathGen)
{
    static const string routine = "SimpathIInstrument::pathGenUpdated";
    try {

        SVGenSpot::IStateVarSP fxSV;

        for (SimPathICollectionMapIterator it = collectionMap.begin();
            it != collectionMap.end();
            ++it)
        {
            if (!(fxSV.get()) && it->second->getFxEqGen().get())
            {
                // TODO: eventually introduce single-asset state vars for FX
                // to be initialized only once for now
                fxSV = it->second->getFxEqGen()->getSpotSV(newPathGen);
            }
            it->second->initSV(newPathGen);
            it->second->setFxEqSV(fxSV);
        }

    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

void SimpathIInstrumentMC::collectStateVars(
        IStateVariableCollectorSP svCollector) const
{
    for (SimPathICollectionMapIterator it = collectionMap.begin();
        it != collectionMap.end();
        ++it)
        it->second->appendSVGen(svCollector);
}

RM_Assets::IMCWriterSP SimpathIInstrumentMC::getWriter() const
{
	return mcWriter;
}

void SimpathIInstrumentMC::finalize()
{
	getWriter()->finalize();
}

//////////////////////////////////////////////////////////////////////////////////

class SimpathiPayoffEvent : public IPayoffEvent
{
    IDiffusionRecordSP rec;
	int pathIdx;
	int iStep;
public:
	SimpathiPayoffEvent(int _pathIdx, int _iStep, IDiffusionRecordSP r) : pathIdx(_pathIdx), iStep(_iStep), rec(r) {}
	const IDiffusionRecordSP getRecord() const {return rec;}
	int getPathIdx() const { return pathIdx;}
	int getIStep() const {return iStep;}
};

void SimpathIInstrumentMC::payoff(
        const MCPathGenerator* pathGen,
        IMCPrices& _prices)
{
    const int pathIdx = pathGen->getPathIndex();
	IGenericPrices* prices = dynamic_cast<IGenericPrices* >(& _prices);
	
	mcWriter->notifyStartPath(pathIdx); // notify that we generate a new path

    for(int iStep=0; iStep != timeLine.size(); ++iStep)
    {
		const DateTime & iStepDate = timeLine[iStep];
		
        for (SimPathICollectionMapIterator it = collectionMap.begin();
            it != collectionMap.end();
            ++it)
        {
			/* The logic is the following:
			   we loop over all assetInfos and update them (i.e. store diffusion results for the current path);
			   at the end of each date we flush data into the writer (i.e. go over all assets and call write_one_record(...))

			   The plan is to refactor it as follows:
			   Results of diffusion are captured in an IDiffusionRecord object (derived from CObject).
			   Writer updates state of assetInfo using this object only (i.e. updateDiffusedState accepts IDiffusionRecord).
			   The DiffusionResults live in Qlib's realm and are collected to be returned as "Prices".

			   Currently, to make modifications minimal, the IDiffusionRecord is created from the assetInfo.
			 */
			
            ISimPathICollectionSP collection = it->second;
            collection->updateDiffusedAssetInfo(iStep, iStepDate);

			// In the regression test framework we also want to report values back.
			// In SAMPRAS we definitely do not want to do this:
			// 1. the amount of data if huge (100M per path)
			// 2. we're in the inner loop, so it can be the performance problem;
			// How updatePrices is set ? If no writers is requested => we set doPrices.
			if (prices != NULL) {
				IDiffusionRecordSP rec = collection->getAssetInfo()->create_record();
				prices->add(SimpathiPayoffEvent(pathIdx, iStep, rec), 1.0);
			}
        }

        mcWriter->notifyEndDate(pathIdx, iStep, iStepDate);

    }
    mcWriter->notifyEndPath(pathIdx); // notify that we generate a new path
}

class SimpathiFakePricesMC : public IMCPrices {
public:
	virtual void add(double) {}
	virtual void reset() {}
	virtual int storagePerPath(IMCProduct* product) const { return -1;}
	virtual void configureCache(const IntArray& changedAssets) {}
	/** Returns a deep copy of this object */
	virtual IMCPrices* clone() const {
		QLIB_VERIFY(0, "Not implemented");
		return NULL;
	}
    virtual IMCPrices* emptyConstructor() const {
        return new SimpathiFakePricesMC;
    }
};
	
//////////////////////////////////////////////////////////////////////////////////////

class SimpathiPricesMC : public IMCPricesGeneric {
	IDiffusionRecordArrayArrayArraySP events; // indexed via [path][time point][asset]
	void process(const SimpathiPayoffEvent& rec, double weight) {
		int i = rec.getPathIdx();
		int j = rec.getIStep();
		// boring code to ensure elements[path][timepoint][asset] is valid
		// 1. ensure  elements[i] points to an array
		if (events->size() <= i)
			events->resize(i+1);
		if ((*events)[i].get() == NULL)
			(*events)[i] = IDiffusionRecordArrayArraySP(new IDiffusionRecordArrayArray);
		// 1. ensure  elements[i][j] points to an array
		if ((*events)[i]->size() <= j)
			(*events)[i]->resize(j+1);
		if ((*(*events)[i])[j].get() == NULL)
			(*(*events)[i])[j] = IDiffusionRecordArraySP(new IDiffusionRecordArray);
		
		(*(*events)[i])[j]->push_back(rec.getRecord());
	}

	int             nbIter;
	int             nbSubSamples;
	int             numMatDates;
public:
	SimpathiPricesMC(
		int             NbIter,
		int             NbSubSamples,
		int             NumMatDates // number of the maturity dates
		):
		nbIter(NbIter),
		nbSubSamples(NbSubSamples),
		numMatDates(NumMatDates),
		events(new IDiffusionRecordArrayArrayArray)
		{
			reset();
		}
	//////////////  IGenericPrices interface 
	virtual void add(const IPayoffEvent& ev, double weight )
		{
			const SimpathiPayoffEvent& mattr =
				dynamic_cast<const SimpathiPayoffEvent&> (ev);
			process(mattr, weight);
	}

	virtual IObjectSP getResult() const {
		return events;
	}

	////////////////// IGreeks interface ////////////
	virtual void reset() {
		events = IDiffusionRecordArrayArrayArraySP(new IDiffusionRecordArrayArrayArray);
	}

	virtual int storagePerPath(IMCProduct* product) const { return -1;}
	virtual void configureCache(const IntArray& changedAssets) {}

	///////////////// IMCPrices interface ////////
public:
	/** Returns a deep copy of this object */
	virtual IMCPrices* clone() const {
		QLIB_VERIFY(0, "Not implemented");
		return NULL;
	}

protected:
	/** Ease cloning */
	virtual IMCPrices* emptyConstructor() const { 
		return new SimpathiPricesMC(1, 1, numMatDates);
	}
public:
	// to be called from MCProduct::recordExtraOutput
	void recordExtraOutput(CControl *control, Results * result) const
	{
		// FIXME
	}
};

IMCPrices*  SimpathIInstrumentMC::createOrigPrices(int nbIter,
                                                   int  nbSubSamples,
                                                   int  mode)
{
	// doPrices is set when no writers is requested; in this case we will accumulate prices and return them as an IObject
	// for technical reasons, in the other case we still need to return something.
    return (doPrices) ? static_cast<IMCPrices*>(new SimpathiPricesMC(nbIter, nbSubSamples, timeLine.size())) : 
                        static_cast<IMCPrices*>(new SimpathiFakePricesMC());
}

IObjectSP SimpathIInstrumentMC::getResults()
{
    IObjectSP results(CBool::create(false));
    RM_Assets::IMCWriter * writer = mcWriter.get();
//    RM_Assets::MCXMLWriter * xmlWriter = dynamic_cast<RM_Assets::MCXMLWriter *> (writer);
    RM_Assets::MCRegTestWriter * xmlWriter = dynamic_cast<RM_Assets::MCRegTestWriter *> (writer);
    if (xmlWriter != NULL)
        results = xmlWriter->getResults();
    return results;
}

int SimpathIInstrumentMC::getRngIdx(int globalIdx, const string& assetName) const
{
	if (asset2rng.find(assetName) == asset2rng.end())
	{
		cerr << "ERROR: No RNG slot information available for the asset [" + assetName + "]" << endl;
		cerr << "Available slots are:\n";
		for(map<string, int>::const_iterator it= asset2rng.begin();
			it != asset2rng.end();
			++it)
			cerr << "Name =[" << it->first << "] rngSlot= " << it->second << endl;
	}
	
	QLIB_VERIFY(asset2rng.find(assetName) != asset2rng.end(), "ERROR: No RNG slot information available for the asset [" + assetName + "]");
	int rngSlot = asset2rng.find(assetName)->second;
	QLIB_VERIFY(rngSlot >= 0, "ERROR: asset [" + assetName +"] has negative slot index= " + Format::toString(rngSlot));
	
#if 0
    cerr << "SimpathIInstrumentMC::getRngIdx(" << assetName << ")= " << rngSlot << endl;
#endif
	
	return rngSlot;
}
DRLIB_END_NAMESPACE

