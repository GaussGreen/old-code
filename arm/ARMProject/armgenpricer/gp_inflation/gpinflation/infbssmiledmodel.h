/*----------------------------------------------------------------------------*/
 
/*! \file infbssmiledmodel.h
 * Copyright (c) NATIXIS CM April 2007 Paris
 *
 *  \brief Sabr model for the inflation volatility 
 *  + vol Atm for the interest rate
 *  + corelation between 2 atm volatilities ( Ir and Inf) 
 *
 *	\author  Mathieu Bernardo
 *	\version 1.0
 *	\date April 2007
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFLATION_INFBSSMILEDMODEL_H
#define _INGPINFLATION_INFBSSMILEDMODEL_H

#include "mod/bssmiled.h"			/// because of inheritance

#include "mod/bsmodel.h"			/// because of inheritance
#include "infoptionmodel.h"			/// because of inheritance
#include "gpinflation/infidx.h"
#include "gpinflation/infcapfloor_.h"

/// forward declaration of ARM Kernel objects
class ARM_ZeroCurve;
class ARM_VolCurve;

CC_BEGIN_NAMESPACE( ARM )



//////////////////////////////////////////////////////////
///  \class ARM_InfBSSmiledModel
///  \author  Mathieu Bernardo
///  \version 1.0
///  \date April 2007
/// 
///  compute cap and swaption
/// 
//////////////////////////////////////////////////////////

class ARM_InfBSSmiledModel : public ARM_BSSmiledModel, public InfOptionModel{     
	
private :

	ARM_VolCurve*		itsCorrelVolCurve;
	ARM_VolCurve*		itsIrCapVolCurve;
	ARM_VolCurve*		itsCorrelAdjVolCurve;

public:
	ARM_InfBSSmiledModel(	ARM_Date&			asOfDate,
							ARM_ZeroCurve*		discountCurve, 
							ARM_InfCurv*		infFwdCurve,
							ARM_VolCurve*		irCapVolCurve,
							ARM_VolCurve*		correlVolCurve,
							ARM_VolCurve*		correlAdjVolCurve,
							ARM_VolCurve*		infCapVolCurve,
							ARM_VolCurve*		nuCurve,
							ARM_VolCurve*		rhoCurve,
							ARM_VolCurve*		betaCurve,
							int					sabrMode			= 3,
							int					alphaOrSigmaInput	= 1);

	ARM_InfBSSmiledModel(const ARM_InfBSSmiledModel& rhs);
	ASSIGN_OPERATOR		(	ARM_InfBSSmiledModel		);
	virtual ~ARM_InfBSSmiledModel();
	virtual ARM_Object* Clone(){ return new ARM_InfBSSmiledModel(*this);}
	virtual void View(char* id = NULL, FILE* ficOut = NULL);
	
	ARM_VolCurve*	GetAdjConvCorrel() const { return itsCorrelAdjVolCurve; }

	double SingleAssetOptionPrice(			double			CPIForward,
											double			strike,
											int				callPut,
											double			discounting, 
											ARM_InfIdx* infIdx,
											InfCapFloorContext* optionContext,
											StoreVolInfoWithNewCF * storeInfo	);


	virtual double SingleAssetOptionPrice(	double CPIForward,
											double strike,
											int callPut,
											double discounting, 
											ARM_InfIdx* infIdx,
											ARM_Object* optionContext,
											StoreInfoObj& storeInfo	);

	virtual double TwoAssetsOptionPrice(	double			CPIForward, 
											double			secondAssetFwd, 
											double			strike,
											int				callPut, 
											double			discounting, 
											ARM_InfIdx*		infIdx, 
											ARM_IRIndex*	secondIndex,
											ARM_Object*		optionContext, 
											StoreInfoObj&	storeInfo	){
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSSmiledModel::TwoAssetsOptionPrice : not implemented");	}

	virtual double			GetModelTimeWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx);
	virtual ARM_Date		GetModelDateWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx);

	virtual ARM_GP_Vector CptFwdCPIRatio(	const double& numDate,
											const double& demDate,
											const double& paymentDate,
											long		  dailyInterpType,
											double		  denomFixing	= GETDEFAULTVALUE,
											ARM_InfIdx*	  infIdx		=NULL);

	virtual ARM_GP_Vector FwdCPIRatio(		const ARM_Date& numDate,
											const ARM_Date& denomDate,
											const ARM_Date& paymentDate,
											double			multiple,
											double			spread,
											long			dailyInterpType,
											double			denomFixing	= GETDEFAULTVALUE,
											ARM_InfIdx*		infIdx		=NULL);

	void ComputeVolAndStrikeForCap(			double			CPIForward,
											double			strike,
											int				callput,
											ARM_InfIdx* infIdx,			
											StoreVolInfoWithNewCF* storeInfo,
											const ARM_Date& numDate,	
											const ARM_Date& denomDate,
											double&			totalVol,
											double&			pricingStrike,
											double&			tenor) ;

	void ComputeVolAndStrikeForCap(			double			CPIForward,
											double			strike,
											int				callput,
											ARM_InfIdx*		infIdx,			
											StoreInfoObj&	storeInfo,
											const ARM_Date& numDate,	
											const ARM_Date& denomDate,	
											int				optionType,				
											double&			totalVol,
											double&			pricingStrike,
											double&			tenor ) ;
	

	virtual double ComputeVol		(		double	maturity, 
											double	volTenor, 
											double	fwd, 
											double	strike, 
											int		mode = -1);

	virtual double DiscountedCPI(			const ARM_Date& resetDate, 
											const ARM_Date& paymentDate, 
											long			dailyInterpType,
											const string&	CPILag,
											const string&	DCFLag,
											ARM_InfIdx*		infIdx		 = NULL ) ;
	
	virtual ARM_INF_PRICING_INFO CanPriceInflation() const { return PRICE_FWDNOPTION; }

private:

	void	fillTimesForCalculation(	const ARM_Date& denomDate,  
										const ARM_Date& numDate, 
										const ARM_Date& paymentDate, 
										ARM_InfIdx* infIdx,
										double& demLag,
										double& numLag,
										double& payLag, 	
										double& tenLag);

	
	ARM_VolCurve*	GenerateAdjCorrelMatrix( );

};



CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
