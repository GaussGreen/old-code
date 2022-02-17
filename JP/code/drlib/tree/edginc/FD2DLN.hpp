//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DLN.hpp
//
//   Description : two factor finite difference engine for two LN
//                    asset factors   dX_i = drift_i*dt + sigma_i * dW_i with i=1,2
//
//   Author      : Xiaolan Zhang
//
//   Date        : May 27, 2005
//
//----------------------------------------------------------------------------

#ifndef FD2DLN_HPP
#define FD2DLN_HPP

#include "edginc/FD2D.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/VolProcessedBS.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD2DLN : public FD2D{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    static IObject* defaultFD2DLN(){
        return new FD2DLN();
    }
    
    friend class FD2DSolver;

    FD2DLN();
    virtual ~FD2DLN();

    /** validate some inputs */
    virtual void validatePop2Object();

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

    CVolProcessedBSConstSP getProcessedVol();
    
    
    /** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
    virtual double getPrice0( const TreeSlice & price ) const;
    
    
    /** get the initial conditions and corresponding indexes in the slices
        for forward induction */
    virtual void getInitialConditions(IntArray& initialIndex, DoubleArray& initialValue) const{};
    
    
protected:
    /** Invoked after instrument has got its market data. */
    virtual void getMarket(const MarketData* market, IInstrumentCollectionSP instruments);

    /** retrieving market data from MDF */    
    virtual void retrieveFactor();
    
    virtual void initModel();

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
    
    virtual  void setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1,
                              double& volForBound2, double alpha2, double& outLowB2, double& outUpB2);


    /**-----------------------------*/
    /** update payoffIndex          */
    /** ----------------------------*/
   
    virtual void payoffIndexUpdate (int& step, FDProduct::UpdateType);

    
    /**----------------------------------------------------------------
       setup the coefficients of the PIDE needed to be solved
       Assume the most general PDE is 
       U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
       this ft is called at each time step, so they're time depd.
       ----------------------------------------------------------------*/
    
    /** assuming that DimX (vertical) is stock and DimY is Variance*/
    
    void pdeCoeff(int step, double*** coeff,  
                  int bot1, int top1, int bot2, int top2 );
    
    /** override FDModel::makeProduct to retrieve spot index product */
    virtual FDProductSP makeProduct( const IProdCreatorSP & creator );

    // ------------------ data ------------------------------------------------------

    CAssetConstSP   underlying2; // $unregistered

    FDProductSP     payoffIndex2; // $unregistered

    // log-normal vol
    string volType; // $unregistered

    CVolProcessedBSSP volLN1; // $unregistered
    CVolProcessedBSSP volLN2; // $unregistered

    //correlation
    double rho; // $unregistered

private:
    void generateDriftAndVols(const DateTimeArray& futurePathDates);

    vector<DoubleArray>       vols; // $unregistered
    vector<CVolRequestLNSP> volRequests;    // [numAssets]  $unregistered
    vector<DoubleArray>       drifts;     /* x[numAssets] each with // $unregistered
                                              NbVolReq(iAsset) cols x
                                              NbSteps rows */
};

typedef smartPtr<FD2DLN> FD2DLNSP;

DRLIB_END_NAMESPACE
#endif
