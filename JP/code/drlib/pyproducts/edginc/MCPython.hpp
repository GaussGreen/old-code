#ifndef QLIB_MCPYTHON_HPP
#define QLIB_MCPYTHON_HPP

#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/SVGenDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

class PYPRODUCTS_DLL MCPython: public GenericNFBase,
                               virtual public IMCIntoProduct
{
private:
    friend class MCPythonProdPy;
    /// exported fields ////////
    bool                isCall;
    bool                isBest;
    double              strike;
    DateTimeArray       averageOutDates;

    CStringArray        srcCode; // python source code

public:
    static CClassConstSP const TYPE;
    friend class MCPythonProd;
   
    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return averageOutDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    MCPython(): GenericNFBase(TYPE) {} // for reflection
    MCPython(const MCPython& rhs); // not implemented
    MCPython& operator=(const MCPython& rhs); // not implemented

    static IObject* defaultMCPython(){
        return new MCPython();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(MCPython, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultMCPython);
        FIELD(isCall,     "is it a call option");
        FIELD(isBest,     "True for Best, False for Worst.");
        FIELD(strike,     "strike level");
        FIELD(averageOutDates,  "averageOutDates");
        FIELD(srcCode,     "source code text");
        FIELD_MAKE_OPTIONAL(srcCode);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super MCPython */
class MCPythonProd: public MCProductClient,
                  virtual public IMCProductLN{
private:
    friend class MCPythonProdPy;
    const MCPython*         inst;    // reference to original instrument
    DoubleArray             sumOut;  // [iAsset] working area
    DoubleArray             sumSoFar; // preservation area (past)
    double                  CoP;     // +1 or -1 for call/put 
    int                     nbAssets; // for convenience

    bool                    useSrcCode; // for testing performance, to be removed

    SVGenSpot::IStateVarSP    spotSV;        // asset state variable
    SVGenSpotSP               spotGen;       // generator for asset
    IRefLevel::IStateVarGenSP    refLevelGen;  //!< Generator for ref level
    SVGenDiscFactorSP            dfGen;        //!< Generator for discount factors
    IRefLevel::IStateVarSP       refLevelSV;   //!< Ref level state variable
    SVDiscFactorSP               dfSV;         //!< Df state variable

    /** Returns the maximum delta scaling factor for each asset (it is the
        largest absolute derivative of the payoff wrt each simulated asset)
        Notional is excluded from the calculation */
    DoubleArray maxDeltaFactors(const IPathGenerator* futurePathGen) const{
        // here the largest weight, by definition, is 100%
        DoubleArray maxDeltaFactor(getNumAssets(), 1.0);
        for (int i = 0; i < maxDeltaFactor.size(); i++) {
            maxDeltaFactor[i] /=
                futurePathGen->refLevel(i, 0 /* path irrelevant*/);
        } 
        return maxDeltaFactor;
    }
public:
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor
        and then (apart from storing reference to MCPython) create the
        IMCPerf object which is the combination of the performance data as
        specified by the instrument together with the list of all simulation
        dates and their historic values */
    MCPythonProd(const MCPython*            inst,
                  const SimSeriesSP&        simSeries):
        MCProductClient(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
                  inst(inst),
                  nbAssets(getNumAssets()),
      // state var support
        spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), 
                                                      getToday())),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
            inst->instSettle, simSeries->getLastDate()))
    {
            sumOut.resize(nbAssets);
            sumSoFar.resize(nbAssets);
            CoP = inst->isCall? 1 : -1;
            //if (inst->srcCode.size()==0)
            //    throw ModelException("no python src code supplied");
            useSrcCode = inst->srcCode.size()>0;
    }

    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        svCollector->append(refLevelGen.get());
        svCollector->append(spotGen.get());
        svCollector->append(dfGen.get());
    }
    void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen);

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*, IMCPrices& prices);

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray volReq(1);
        const DateTime&    startDate = getRefLevel()->getAllDates().front();
        bool               fwdStarting = startDate.isGreater(getToday());
        double             interpLevel;
        const DateTime&    lastSimDate = getSimSeries()->getLastDate();

        if (fwdStarting){
            interpLevel = inst->strike;
        } else {
            /* not forward starting - some samples have fixed already
                (this includes averaging in) */
            int numDates = inst->averageOutDates.size();
            int numRemaining = getToday().numFutureDates(inst->averageOutDates);
            
            interpLevel = (numDates * inst->strike * pathGen->refLevel(iAsset, 0)
                            - sumSoFar[iAsset])/ numRemaining;
        }
        volReq[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interpLevel,
                                                                    startDate,
                                                                    lastSimDate,
                                                                    fwdStarting));

        return volReq;
    }
};

DRLIB_END_NAMESPACE

#endif