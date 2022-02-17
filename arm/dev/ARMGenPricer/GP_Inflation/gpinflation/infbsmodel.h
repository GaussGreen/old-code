/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: infbsmodel.h,v $
 * Revision 1.10  2003/12/15 18:20:36  ebenhamou
 * change for new compiler
 *
 * Revision 1.8  2003/11/20 15:27:03  ebenhamou
 * version with dates.
 *
 * Revision 1.7  2003/09/23 17:35:26  ebenhamou
 * convexity adjustment
 *
 * Revision 1.6  2003/09/22 18:11:54  ebenhamou
 *  added functionality for adj
 *
 * Revision 1.5  2003/09/11 11:01:32  ebenhamou
 * factorisation of code and past reset
 *
 * Revision 1.4  2003/09/10 17:09:06  ebenhamou
 * the CanPriceInflation has to be in the class.. remove multiple inheritance ambiguity
 *
 * Revision 1.2  2003/09/09 11:45:28  ebenhamou
 * added interface for option pricing
 *
 * Revision 1.1  2003/09/02 17:35:37  ebenhamou
 * Initial revision
 *
 *
 */

/*----------------------------------------------------------------------------*/
 
/*! \file infbsmodel.h
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \brief Very simple model, based on Back Scholes that knows a volcube a curve
 *	of forward CPI and a discount curve
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFLATION_INFBSMODEL_H
#define _INGPINFLATION_INFBSMODEL_H

#include <mod/bsmodel.h>			/// because of inheritance
#include "infoptionmodel.h"		/// because of inheritance
#include "infswopvol.h"		/// because of inheritance
#include "gpinflation/infcapfloor.h"
#include "gpbase/gpvector.h"

/// forward declaration of ARM Kernel objects
class ARM_ZeroCurve;
class ARM_VolCurve;
class ARM_IRIndex;
class ARM_CorrelManager;
class ARM_VolCube;
struct ARM_CorrelMatrix;

/// ARM namespace
CC_BEGIN_NAMESPACE( ARM )

/// forward declaration of ARM_GP objects
class ARM_InfCurv;
class ARM_InfCapFloorContext;

////////////////////////////////////////////////////
///	Class  : StoreVolInfo
///	Action : Helper function to store 
///				some pricing information
////////////////////////////////////////////////////

class StoreVolInfo : public StoreInfoObj
{
public:
	StoreVolInfo( double strike, double timeToStart, double tenor, 
		int k, ARM_InfCapFloor* infCapFloor )
	: itsStrike( strike ), itsTimeToStart( timeToStart ), 
	itsTenor( tenor ), itsK( k ),	itsInfCapFloor( infCapFloor )
	{}
	
	virtual void Store( double* data )
	{
		double vol				= data[0];
		double volLookupStrike	= data[1];
		itsInfCapFloor->StoreVol( vol, volLookupStrike, itsStrike, itsTimeToStart, itsTenor, itsK );
	}
	ARM_InfCapFloor* GetInfCapFloor() const { return itsInfCapFloor; }
private:
	double	itsStrike;
	double	itsTimeToStart;
	double	itsTenor;
	int		itsK;
	ARM_InfCapFloor* itsInfCapFloor;
};

//////////////////////////////////////////////////////////
///  \class ARM_InfBSModel
///  \author  Eric Benhamou
///  \version 1.0
///  \date August 2003
/// 
///  \brief an inflation bs model is composed of
///  a yield curve, a cpi forward curve
///  and a vol cube, which is derived from an ARM_ZeroCurve
///  It is a very simple model used to price only vanilla inflation
///  products
/// 
///  the assumptions are the one of a standard BS model, i.e.
///  \f[
/// 		\frac{dCPI_{t}}{CPI_{t}}=\sigma \left\( t,K\right\) dW_{t}
///  \f]
///  with \f$CPI_{t}\f$ the forward CPI
/// 
//////////////////////////////////////////////////////////

class ARM_InfBSModel : public ARM_BSModel, public InfOptionModel
{        
protected :
	/// the cap vol curve is saved into the Black Scholes model
	ARM_CorrelManager*	itsCorrelManager;
	ARM_BSModel*		itsIRModel;
	ARM_VolCurve*		itsInfSwoptVolCurve;
	ARM_VolCurve*		itsIRSwoptVolCurve;
	ARM_VolCube*		itsIRSwoptVolCube;

	bool				switchToIRSwoptVolCurve;
	bool				infVolHasBeenComputed;
	bool				irVolHasBeenComputed;
	const ARM_Date		itsAsOfDate;
	bool				isCorrelMatrixValidated;
	bool				isGoodCorrelMatrix;

	/// function for the cap vol and strike
	void ComputeVolAndStrikeForCap(
		double CPIForward,
		double strike,
		int callPut,
		ARM_InfIdx* infIdx,			
		StoreInfoObj& storeInfo,
		const ARM_Date& numDate,	
		const ARM_Date& denomDate,	
		int optionType,				
		double renormalisationFactor,
		double& totalVol,
		double& pricingStrike,
		double& tenor );

	void CleanUp();	

public:
	ARM_InfBSModel(){}
	//// correlManager, IRModel and SwoptCurve have null arguments to
	/// allow a simple version of the InfBSModel
	ARM_InfBSModel(	const ARM_Date&		asOfDate,
					ARM_ZeroCurve*		discountCurv, 
					ARM_InfCurv*		infFwdCurv,
					ARM_VolCurve*		infCapVolCurv,
					ARM_CorrelManager*	CorrelManager	= NULL,
					ARM_BSModel*		IRModel			= NULL,
					ARM_VolCurve*		infSwoptCurve	= NULL,
					ARM_VolCurve*		IRSwoptCurve	= NULL );

	//// another constructor added to include the smile to the IR part
	ARM_InfBSModel(	const ARM_Date&		asOfDate,
					ARM_ZeroCurve*		discountCurv, 
					ARM_InfCurv*		infFwdCurv,
					ARM_VolCurve*		infCapVolCurv,
					ARM_CorrelManager*	CorrelManager	= NULL,
					ARM_BSModel*		IRModel			= NULL,
					ARM_VolCurve*		infSwoptCurve	= NULL,
					ARM_VolCube*		IRSwoptCube		= NULL );

	ARM_InfBSModel(const ARM_InfBSModel& rhs);
	ARM_InfBSModel& operator = (const ARM_InfBSModel &rhs );
	virtual ~ARM_InfBSModel();

	//// standard ARM Object support
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	virtual ARM_GP_Vector CptFwdCPIRatio(	const ARM_Date& numDate,
											const ARM_Date& denomDate,
											const ARM_Date& paymentDate,
											long			dailyInterpType,
											double			denomFixing	= GETDEFAULTVALUE,
											ARM_InfIdx*		infIdx		= NULL){ 
		return FwdCPIRatio( numDate,denomDate,	paymentDate, 1.0,-1.0,	dailyInterpType, denomFixing,	infIdx); }
	
	/// standard pricing of forward CPI
	virtual ARM_GP_Vector FwdCPIRatio( 
		const ARM_Date& numDate,
		const ARM_Date& denomDate,
		const ARM_Date& paymentDate,
		double multiple,
		double spread,
		long dailyInterpType,
		double denomFixing		= GETDEFAULTVALUE,
		ARM_InfIdx* infIdx		= NULL );

	virtual double DiscountedCPI(
		const ARM_Date& resetDate, 
		const ARM_Date& paymentDate, 
		long dailyInterpType,
		const string& CPILag,
		const string& DCFLag,
		ARM_InfIdx* infIdx		= NULL );


	virtual ARM_INF_PRICING_INFO CanPriceInflation() const { return PRICE_FWDNOPTION; }

	virtual double GetModelTimeWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx );

	virtual double SingleAssetOptionPrice(
		double CPIForward,
		double strike,
		int callPut,
		double discounting, 
		ARM_InfIdx* infIdx,
		ARM_Object* optionContext,
		StoreInfoObj& storeInfo	);
	
	virtual double TwoAssetsOptionPrice(
		double CPIForward,
		double secondAssetFwd, 
		double strike,
		int callPut, 
		double discounting,
		ARM_InfIdx* infIdx,
		ARM_IRIndex* secondIndex,
		ARM_Object* optionContext, 
		StoreInfoObj& storeInfo	);

	virtual ARM_CorrelMatrix* GetInfIRCorrelMatrix ( ARM_InfIdx* infIdx, ARM_IRIndex* otherIndex, const string& mktTag  ) const;
	virtual ARM_CorrelMatrix* GetInfSwpCorrelMatrix( ARM_InfIdx* infIdx,  const string& mktTag  ) const;
	virtual ARM_CorrelMatrix* GetIRSwpCorrelMatrix( ARM_IRIndex*  irIdx,  const string& mktTag  ) const;

	virtual double GetInfIRCorrel( double TjFromAsOf, double TiFromAsOf, 
		ARM_InfIdx* infIdx, ARM_IRIndex* otherIndex, const string& mktTag  ) const;
	virtual double GetInfSwpCorrel( double tenor1, double tenor2, ARM_InfIdx*	infIdx,	const string& mktTag  ) const;
	virtual double GetIRSwpCorrel( double tenor1, double tenor2,  ARM_IRIndex*	 irIdx, const string& mktTag  ) const;

	virtual ARM_CorrelMatrix* GetInfInfCorrelMatrix( ARM_InfIdx* infIdx1, ARM_InfIdx* infIdx2, const string& mktTag ) const;
	virtual double GetInfInfCorrel( double TjFromAsOf, double TiFromAsOf, 
		ARM_InfIdx* infIdx1, ARM_InfIdx* infIdx2, const string& mktTag  ) const;

	static double ComputeVolForSwaption(
		double strike,
		double optionMaturity,
		double swaptionMaturity,
		ARM_VolCurve* volCurve,
		const string& volName );

	static double ComputeVolForSwaptionWithSmile(
		double strike,
		double optionMaturity,
		double swaptionMaturity,
		ARM_VolCube* volCube,
		const string& volName );

	static double ComputeVolForOATSwaption(
		double strike,
		double coupon,
		double optionMaturity,
		double swaptionMaturity,
		ARM_VolCurve* volCurve,
		const string& volName );

	/// specific to the variable notional case
	// inf part
	double ComputeInfVolWithVarNotional(
		const ARM_GP_Vector& communNotional,
		double floatStartTime,		
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatPayTimes,
		const ARM_GP_Vector* dfTerms,	
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdRates,
		double strike,
		ARM_InfIdx* infIdx,
		ARM_Object* optionContext,
		ARM_GP_Vector& coeffs,
		ARM_GP_Vector& tenors,
		ARM_GP_Vector& vols);


	static void ComputeInfFlowsFromSwapRates(
		const ARM_GP_Vector&  dfTerms,
		const ARM_GP_Vector&  fixPayPeriods,
		const ARM_GP_Vector&  swpRates,
		ARM_GP_Vector& infFlows);

	static void ComputeInfSensiCoefs(
		const ARM_GP_Vector&  dfTerms,
		const ARM_GP_Vector&  fixPayPeriods,
		const ARM_GP_Vector&  infCoefs,
		ARM_GP_Vector&  swpRates,
		const double& shift,
		const double& initVlue,
		ARM_GP_Vector& sensiCoefs);

	// ir part
	double ComputeIRVolWithVarNotional(
		const ARM_GP_Vector& irNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& irFixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		double strike,
		ARM_IRIndex* irIdx,
		ARM_Object* optionContext,
		ARM_GP_Vector& coeffs,
		ARM_GP_Vector& tenors,
		ARM_GP_Vector& vols);

	double ComputeInfIrCovarWithVarNotional(
			ARM_GP_Vector& infcoeffs,
			ARM_GP_Vector& inftenors,
			ARM_GP_Vector& infvols,
			ARM_GP_Vector& ircoeffs,
			ARM_GP_Vector& irtenors,
			ARM_GP_Vector& irvols,
			ARM_InfIdx* infIdx,
			ARM_IRIndex* otherIndex);

	static void ComputeDiscountFactorsFromSwapRates (
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& swapRates,
		double startDf, 
		ARM_GP_Vector& dfs);

	static void ComputeIRSensiCoefs(
		const double dfCoefStart,
		const double dfStart,
		const ARM_GP_Vector& fixCoefs,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector&  dfCoefs,
		ARM_GP_Vector&  swpRates,
		const double& shift,
		const double& initVlue,
		ARM_GP_Vector& sensiCoefs);

	/// Commun
	ARM_VolCurve* GetInfSwoptVolCurve()   const { return itsInfSwoptVolCurve; }
	ARM_VolCurve* GetIRSwoptVolCurve()	  const { return itsIRSwoptVolCurve; }
	ARM_VolCube* GetIRVolCube()			  const { return itsIRSwoptVolCube; }
	ARM_VolCube* GetIRGVolCube()		  const { return itsIRSwoptVolCube; }
	ARM_CorrelManager* GetCorrelManager() const { return itsCorrelManager; }
	void SetIRSwoptVolCurve( ARM_VolCurve* vol)	{ itsIRSwoptVolCurve = vol; }

	ARM_BSModel* GetIRModel ( ) const				{ return itsIRModel; }
	void		 SetIRModel ( ARM_BSModel* mod )	{ itsIRModel=mod ;  }
	
    virtual ARM_VolCurve* GetVolatility(int mode = -1 ) const;
	void ValidateCorrelManagerForVarNotionalSwaption(ARM_InfIdx* infIdx, ARM_IRIndex* otherIndex );
	bool GetCorrelMatrixValidated() {return isCorrelMatrixValidated;}
	void SetCorrelMatrixValidated(bool YesOrNo) {isCorrelMatrixValidated =YesOrNo;}

	bool GetisGoodCorrelMatrix() {return isGoodCorrelMatrix;}
	void SetisGoodCorrelMatrix(bool YesOrNo) {isGoodCorrelMatrix =YesOrNo;}

	ARM_Date GetModelDateWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx );


	inline void SetAsOfDate( const ARM_Date& asOfDate ){ SetStartDate( const_cast<ARM_Date&>(asOfDate) );	}
		
};


CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
