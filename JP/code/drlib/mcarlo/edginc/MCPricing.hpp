#ifndef EDR_MCPRICING_HPP
#define EDR_MCPRICING_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MCPathGenerator.hpp"  // Why pass in MCPathGeneratorSP to createFuturePrices
#include "edginc/LRGenerator.hpp"
#include "edginc/MCPrices.hpp"

DRLIB_BEGIN_NAMESPACE

//class MCPricesSP;
class MonteCarlo;
class IMCProduct;
class Control;
class Results;
class IMCQuickXGamma;


/** Class to allow finer control of the main pricing loop eg
    quuick greeks, quick X gamma etc. Default version is used unless
    the model overrides it */
class MCARLO_DLL MCPricing{
public:
    virtual ~MCPricing();

    /** Creates IMCPrices object for 'past' simulation run */
    virtual IMCPrices* createPastPrices(
        const MonteCarlo*          mcarlo,  // the model
        IMCProduct*                 product);

    /** whether we should ask the path generator to save data such
        that random access to paths can be provided */
    virtual bool needRandomAccessToPaths(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                 product,
        Control*                   control);

    /** Creates IMCPrices object for future simulation run. */
    virtual IMCPricesSP createFuturePrices(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                 product,
        Control*                   control,
        const MCPathGeneratorSP&   pathGen); // path generator

    /** Run the simulation */
    virtual void runSim(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                 product,
        Control*                   control,
        MCPathGenerator*           pathGen, // path generator
        IMCPrices&         prices);

    /** Store the results of the simulation in the results. This
        must include the price */
    virtual void postResults(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                 product,
        Control*                   control,
        IMCPrices&         prices,
        const string&              discountISOCode,
        Results*                   results);

    /** MCPricing structure maintains the list of simulation dates which is
        the master timeline for the diffusion.  Default implementation
        is provided here for backwards compatibility (FIX ME).  */
    virtual DateTimeArray& getTimeLine()
    {
        return simDates;
    }

protected:
    DateTimeArray simDates; // FIX ME
};

typedef refCountPtr<MCPricing> MCPricingSP;

//////////////////////////////////////////////////////////////////////
// Pricers based on stateless payoffs
//////////////////////////////////////////////////////////////////////
class MCARLO_DLL MCSlicePricing: public MCPricing {
public:
    void runSim( MonteCarlo* mcarlo,  // the model
                 IMCProduct* product,
                 Control* control,
                 MCPathGenerator* pathGen, // path generator
                 IMCPrices& prices);

    virtual DateTimeArray& getTimeLine();

private:
     DateTimeArray simDates;
};


class MCARLO_DLL MCPathPricing: public MCPricing {
public:
    void runSim( MonteCarlo* mcarlo,  // the model
                 IMCProduct* product,
                 Control* control,
                 MCPathGenerator* pathGen, // path generator
                 IMCPrices& prices);
    virtual DateTimeArray& getTimeLine();
private:
    DateTimeArray simDates;
};


class MCARLO_DLL MCPastPricing: public MCPricing {
public:
    void runSim( MonteCarlo* mcarlo,  // the model
                 IMCProduct* product,
                 Control* control,
                 MCPathGenerator* pathGen, // path generator
                 IMCPrices& prices) ;
};


//////////////////////////////////////////////////////////////////////
// Specialised QGPricing class for quick greeks
//////////////////////////////////////////////////////////////////////
class MCARLO_DLL QGPricing: public MCPricing{
    /** whether we should ask the path generator to save data such
        that random access to paths can be provided */
    bool needRandomAccessToPaths(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                 product,
        Control*                   control);

    /** Creates IMCPrices object for future simulation run. */
    IMCPricesSP createFuturePrices(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                 product,
        Control*                   control,
        const MCPathGeneratorSP&   pathGen);

    /** Run the simulation */
    void runSim(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                 product,
        Control*                   control,
        MCPathGenerator* pathGen, // path generator
        IMCPrices&         prices);
};
//////////////////////////////////////////////////////////////////////
// end of Specialised QGPricing class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Specialised QXGammaPricing class for quick X gamma
//////////////////////////////////////////////////////////////////////
class MCARLO_DLL QXGamma: public QGPricing{
private:
    IMCPricesSP       qxgPrices;
public:
    QXGamma(const IMCPricesSP& qxgPrices);

    /** Creates IMCPrices object for future simulation run. */
    IMCPricesSP createFuturePrices(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                 product,
        Control*                   control,
        const MCPathGeneratorSP&   pathGen);
};
//////////////////////////////////////////////////////////////////////
// end of Specialised QGPricing class
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// Specialised SubGamma class - for SubGamma (see quick X Gamma)
//////////////////////////////////////////////////////////////////////
/** Class is for calculating what we call SubGamma. SubGamma for asset
    A is the same as Gamma for asset A but it is only done on a subset
    of paths for which the cross gamma wrt assets A and B is non
    zero. Hence for each B there is a different SubGamma for A */
class MCARLO_DLL SubGamma: public MCPricing{
private:
    MCPricingSP  origPricing; // what delta would be calculated with
public:
    bool                        doRealDelta; /* true: calculate real delta
                                                 too */
    vector<MCPricesGeneralSP>   gammaPrices; /* num assets - but
                                                 gammaVals[assetIdx] not
                                                 used */
    DoubleMatrix                 gammaVals;   /* num assets * num assets */
    int                          assetIdx;
    bool                         firstPricingCall;
    vector<bool>                 completed; // per asset
#ifdef _DEBUG
    int                          debugCount; // 2 * number of pricings
#endif

    SubGamma(const MCPricingSP&         origPricing,
             bool                       doRealDelta,
             const MCPricesGeneralSP&   origPrices,
             int                        numAssets);

    /** Creates IMCPrices object for future simulation run. */
    IMCPricesSP createFuturePrices(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                product,
        Control*                   control,
        const MCPathGeneratorSP&   pathGen);

    /** Run the simulation */
    void runSim(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                product,
        Control*                   control,
        MCPathGenerator* pathGen, // path generator
        IMCPrices&         pricesRef);

    /** Store the results of the simulation in the results. This
        must include the price */
    void postResults(
        MonteCarlo*                mcarlo,  // the model
        IMCProduct*                 product,
        Control*                   control,
        IMCPrices&         prices,
        const string&              discountISOCode,
        Results*                   results);

    /** Prepare for calculating another set of subGammas wrt specified
        asset */
    void resetForAsset(
        int                                 assetIdx,
        IMCQuickXGamma*            qckXGamma,
        const MCPathGeneratorSP&  futurePathGen,
        const ScalarShiftArray&             theShifts);

    /** Calcuates all the sub gammas that we need */
    void calculateSubGamma(
        MonteCarlo*                         mcarlo,
        IMCQuickXGamma*            qckXGamma,
        const MCPathGeneratorSP&  futurePathGen,
        Control*                            control,
        CInstrument*                        instrument,
        CResults*                           results,
        const ScalarShiftArray&             theShifts);

};
//////////////////////////////////////////////////////////////////////
// end of SubGamma class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// LRPricing class (do price + greeks using LR method)
//////////////////////////////////////////////////////////////////////
class MCARLO_DLL LRPricing: public MCPricing{
private:
    MCPricingSP              mainPricing;//cooperate with this one (except on sim)
    LRGenerator::GreekCalculatorArray greekCalcs; // holds LR greek calculators

public:
    LRPricing(const MCPricingSP&         mainPricing);

    /** whether we should ask the path generator to save data such
        that random access to paths can be provided */
    bool needRandomAccessToPaths(
        MonteCarlo*                 mcarlo,  // the model
        IMCProduct*                 product,
        Control*                    control);

    /** Creates IMCPrices object for future simulation run. */
    IMCPricesSP createFuturePrices(
        MonteCarlo*                 mcarlo,  // the model
        IMCProduct*                 product,
        Control*                    control,
        const MCPathGeneratorSP&    pathGen);

    /** Run the simulation */
    void runSim(
        MonteCarlo*                 mcarlo,  // the model
        IMCProduct*                 product,
        Control*                    control,
        MCPathGenerator*            pathGen, // path generator
        IMCPrices&                  prices);

    /** Store the results of the simulation in the results. This
        must include the price */
    void postResults(
        MonteCarlo*                 mcarlo,  // the model
        IMCProduct*                 product,
        Control*                    control,
        IMCPrices&                  prices,
        const string&               discountISOCode,
        Results*                    results);
};

//////////////////////////////////////////////////////////////////////
// end of Specialised LRPricing class
//////////////////////////////////////////////////////////////////////

DRLIB_END_NAMESPACE
#endif // EDR_MCPRICING_HPP
