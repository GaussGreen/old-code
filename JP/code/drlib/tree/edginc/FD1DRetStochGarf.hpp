//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRetStochGarf.hpp
//
//   Description : mixing of one factor trinomial tree for local vol process.
//
//----------------------------------------------------------------------------

#ifndef EDG_FD1D_RET_STOCHGARF_HPP
#define EDG_FD1D_RET_STOCHGARF_HPP

#include "edginc/FD1DRet.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/FXVolBase.hpp"
#include "edginc/FD1DRetLV.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1DRetStochGarf : public FD1DRetLV{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    static IObject* defaultFD1DRetStochGarf(){
        return new FD1DRetStochGarf();
    }
    
    friend class FD1DRetStochGarfSolver;

    FD1DRetStochGarf();
    virtual ~FD1DRetStochGarf();

    //virtual TimeMetricConstSP getTimeMetric() const;

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    
    /** this is the main entry point for a model */
    virtual void Price(CInstrument* instrument, 
                                CControl*    control, 
                                CResults*    results);

    // create FD1DRetLVs at here.
    virtual void initModel();

    // call FD1DRetLVs->prod at here.
    virtual void finaliseModel(CControl*    control);

protected:
    /** create a solver*/
    virtual IFDSolverSP createSolver();

private:

    int numOfVols; // $unregistered
    FD1DRetLVArray FD1DRetLVs; // $unregistered
//    string volType;
    vector<FD1DRetLVSP > modelsToUse; // $unregistered

//  get at Price, and use it at each functions...
    CControl* myControl; // $unregistered
    CInstrument* myInst; // $unregistered
};

DRLIB_END_NAMESPACE
#endif
