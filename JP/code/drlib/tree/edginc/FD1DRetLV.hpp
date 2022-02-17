//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRetLV.cpp
//
//   Description : one factor trinomial tree for local vol process.
//
//   Author      : Xiaolan ZHNANG
//
//----------------------------------------------------------------------------

#ifndef FD1D_RET_LV_HPP
#define FD1D_RET_LV_HPP

#include "edginc/FD1DRet.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/FXVolBase.hpp"


DRLIB_BEGIN_NAMESPACE


class TREE_DLL FD1DRetLV : public FD1DRet{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    static IObject* defaultFD1DRetLV(){
        return new FD1DRetLV();
    }
    
    friend class FD1DRetSolver;
    friend class FD1DRetStochGarf;

    FD1DRetLV();
    FD1DRetLV(CClassConstSP type);
    virtual ~FD1DRetLV();

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;
    
    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant: this is a plain old log-normal model.
     * See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    /** get the initial conditions and corresponding indexes in the slices
        for forward induction */
    virtual void getInitialConditions(IntArray& initialIndex, 
                                        DoubleArray& initialValue) const;

protected:

    /** retrieving market data from MDF */    
    virtual void retrieveFactor();

    /**-----------------------------------------------------
        model initialisation  (second part)
        most memory init should be here
        finalise model initialisation once after initModel() and product init() 
        called by validate 
    -----------------------------------------------------*/

    virtual void finaliseModel(CControl*    control);

    /**-----------------------------------------------------
        each derived model can set diff. boundaries 
        based on the dynamics of underlying,
        alpha is the input truncation1D (nb of std) 
        outLowB, outUpB are low and up boundaries for FD 
    -----------------------------------------------------*/

    virtual     void setFdBounds(double& volForBound,  
                                double alpha, 
                                double& outLowB, 
                                double& outUpB);

    /**----------------------------------------------------------------
        PDE 1 factor is
        U_t + a*U_x + c*U_xx + f*U + Jump(x) = g;
    ----------------------------------------------------------------*/

    virtual void pdeCoeff(int step, double** coeff, 
                int bot1, int top1) ;

    /**----------------------------------------------------------------
        calculate vol for current step for a set of spot levels
        returns number of vol calculated - one (flat for all node) or num 
    ----------------------------------------------------------------*/

    int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end);

    void setVolLV(CVolProcessedDVFSP newVolProcessed);

protected:
    CVolProcessedDVFSP VolLV; // $unregistered
    string volType;

    //to see if we need them for Fd
    bool            useTweakingForTimeDerivs; // $unregistered
    double          tweakStrikeUnscaled; // $unregistered
    double          tweakTimeUnscaled; // $unregistered
    double          probDensRatioMin; // $unregistered
    bool            useMidPoint; // $unregistered

    FXAssetConstSP fxAsset; // $unregistered
 
    const Correlation* eqFXCorr; // $unregistered

    CVolProcessedBSSP volFXBS; // $unregistered

    refCountPtr<CVolProcessedDVF::IVolCalculator> volCalculator; // $unregistered

    /**------------------------------------------------------------------------------
        set up variance array 
        For FD, also do initialize the caching of forward values for the local volatility 
    ------------------------------------------------------------------------------*/

    void prepare(double& volForBound);
};

typedef smartPtr<FD1DRetLV> FD1DRetLVSP;
typedef array<FD1DRetLVSP, FD1DRetLV> FD1DRetLVArray;

DRLIB_END_NAMESPACE
#endif
