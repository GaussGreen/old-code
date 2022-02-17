//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : TranchePricerForCalibrator.hpp
//
//   Description : A tranche pricer for SCID calibration
//
//   Author      : Adrian Bozdog
//
//----------------------------------------------------------------------------

#ifndef QR_TRANCHEPRICERFORCALIBRATOR_HPP
#define QR_TRANCHEPRICERFORCALIBRATOR_HPP

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SCID.hpp"
#include "edginc/SCIDparameters.hpp"
#include "edginc/SCIDCalibparameters.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Array.hpp"


DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(TranchePricerForCalibrator)
FORWARD_DECLARE(ICDS)
FORWARD_DECLARE_WRAPPER(ICDSParSpreads)
FORWARD_DECLARE_WRAPPER(YieldCurve)

/** Computes the  MTMs introduced in the Hessian matrix used in the minimization procedure 
 *  Computes MTMs at both tranche level anf portfolio level**/

/** TODO: write a function that updates the matrix **/
    
    
class PRODUCTS_DLL TranchePricerForCalibrator :    public CInstrument,
                                                   public virtual SCID::IIntoProduct {
public:
    static CClassConstSP const TYPE;

    virtual ~TranchePricerForCalibrator();
    /*=========================================================================
     * I/CInstrument Interface
     *=======================================================================*/
    virtual void GetMarket(const IModel* model, const CMarketDataSP);
    virtual void Validate();
    virtual DateTime getValueDate() const;
    virtual string discountYieldCurveName() const;

    /*=========================================================================
     * Pricing model interface implementations
     *=======================================================================*/
    virtual SCID::IProduct* createProduct(SCID* model) const;

    
protected:
    TranchePricerForCalibrator();

private:
    TranchePricerForCalibrator(const TranchePricerForCalibrator& rhs);
    TranchePricerForCalibrator& operator=(const TranchePricerForCalibrator& rhs);

    DoubleArraySP getCalibratedPEL(DoubleArraySP& initialWeights);

    int ComputeMarketConstraints(DoubleMatrix & EqIneq, DoubleArray &EqIneqVal, DoubleMatrix &Min, vector<DoubleArrayArray> &riskyAnnuity, vector<DoubleArrayArray> &defaultLeg); 
	int ComputeWeights(DoubleArraySP &imslWeights, vector<DoubleArrayArray> &riskyAnnuity, vector<DoubleArrayArray> &defaultLeg);
	void getCalibratedLegs(DoubleArraySP& imslWeights, DoubleArrayArraySP& resultRA, DoubleArrayArraySP& resultDL, vector<DoubleArrayArray> &riskyAnnuity, vector<DoubleArrayArray> &defaultLeg);   
    
    /*=========================================================================
     * DATA FIELDS
     *=======================================================================*/
    
    /**Trade value date*/
    DateTime        valueDate;
    
    mutable YieldCurveSP            rfCrv; // $unregistered
    YieldCurveWrapper               discount;        
	int nbLoops;


	// data obtain from the model, and stored for simplicity
	MaturityPeriodSP		lossCalculationInterval;
	int						seed;
	double					timeSteps;
	int						nbPathsNoJump;
	int						nbPathsAtLeastOneJump;
	int						convolutionNoJump;
	int						convolutionAtLeastOneJump;

    DoubleMatrix pel;

    
	SCIDparametersWrapper				sCIDparam;
	SCIDCalibParametersWrapper			sCIDCalibParam;
    /*=========================================================================
     * OTHER METHODS
     *=======================================================================*/
    /**Returns the closed form price of the instrument*/
    void priceClosedForm(
        CResults* results, 
        Control* control, 
        SCID* model);

    /*=========================================================================
     * FRIENDS
     *=======================================================================*/
    friend class TranchePricerForCalibratorHelper;
    friend class TranchePricerForCalibratorSCID;
};

DRLIB_END_NAMESPACE
#endif




