
#ifndef EDG_RISK_MGR_H
#define EDG_RISK_MGR_H
#include "edginc/Control.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

/** Directs calculation of price + sensitivities etc */
class RISKMGR_DLL CRiskMgr{
public:  
    /** Calculates price and greeks etc for possibly two ways of viewing the
        world. One is totally theoretical and the other is using mark to market
        for instruments which could be priced theortically but are also traded
        and therefore have a mark to market (mtm) value. */
    static CResultsArraySP calculateTheoAndMtm(
        IModelSP                mdl,
        IInstrumentCollectionSP instruments,
        CControlSP              control);

    /** Constructs a MultiRiskMgrInterface object and writes it to file. Such an
        object can be used directly with the regression tester */
    static void writeInputs(IModelSP                model,
                            IInstrumentCollectionSP instruments,
                            CControlSP              control,
                            MarketDataSP            market);

    /** Constructs a RiskMgrInterface object and writes it to file. Such an
        object can be used directly with the regression tester */
    static void writeInputs(IModelSP                model,
                            CInstrumentSP           instrument,
                            CControlSP              control,
                            MarketDataSP            market);

    /** Writes an output file containing supplied object. Name is
        based on fileName in control */
    static void writeOutputs(Control*      control,
                             IObjectSP     output);

private:
    CRiskMgr(); // don't allow to be instantiated

    /** Only does one set of calculations (for theo or mtm). Inst must 
        be set up correctly already */
    static CResultsArraySP calculateTheoOrMtm(bool useTheoAssetPrice,
                                              bool derivAssetsExist,
                                              IModelSP model,
                                              IInstrumentCollectionSP instruments,
                                              CControlSP control);

    static void validate(const string& method,
                         IModelSP model,
                         IInstrumentCollectionSP instruments,
                         CControlSP control,
                         MarketDataSP market);

    static string getOutputFileName(Control* control);
};

typedef CRiskMgr RiskMgr;

DRLIB_END_NAMESPACE

#endif
