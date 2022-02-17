/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: infcurvmodel.h,v $
 * Revision 1.10  2003/09/23 17:35:44  ebenhamou
 * interface for convexity adjustment
 *
 * Revision 1.9  2003/09/23 13:52:38  ebenhamou
 * dos2unix
 *
 * Revision 1.8  2003/09/22 18:11:25  ebenhamou
 * added default constructor..
 *
 * Revision 1.7  2003/09/11 11:01:54  ebenhamou
 * factorisation of code to the base class
 *
 * Revision 1.6  2003/09/10 17:08:05  ebenhamou
 * the CanPriceInflation has to be in the class.. remove multiple inheritance ambiguity
 *
 * Revision 1.4  2003/09/02 17:37:47  ebenhamou
 * support for CanPriceInflation
 *
 * Revision 1.3  2003/08/20 08:36:46  ebenhamou
 * added virtual and view method
 *
 * Revision 1.2  2003/08/06 09:28:06  ebenhamou
 * added fast fwfCPI func
 *
 * Revision 1.1  2003/08/05 08:25:10  ebenhamou
 * Initial revision
 *
 */


/*----------------------------------------------------------------------------*/
 
/*! \file infcurvmodel.h
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \brief Class to compute Discounted CPI and FWDCPI
 * combination of infCurve and ZeroCurve
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFLATION_INFCURVMODEL_H
#define _INGPINFLATION_INFCURVMODEL_H

#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"

#include <mod/model.h>
#include "infcurv.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_InfCurv;
class ARM_InfIdx;

/*!
 * \class	InfFwdModel
 * \brief dummy class for pricing of inflation forwards
 *		pure virtual to make it an abstract class
 *		with some default functionalities!
 *		and factorise code between infycmodel and bsmodel
 *
 * \author  Eric Benhamou
 * \version 1.0
 * \date August 2003
 *
 */

class InfFwdModel
{
private:
	/// for safety cloned!
	ARM_InfCurv* itsInfFwdCurve;
	void CleanUp();

public:
	InfFwdModel( ARM_InfCurv* infFwdCurv = NULL );
	InfFwdModel( const InfFwdModel& rhs );
	InfFwdModel& operator=( const InfFwdModel& rhs );
	virtual ~InfFwdModel();

	void SetInfFwdCurv( ARM_InfCurv* infFwdCurv);
	ARM_InfCurv* GetInfFwdCurv( ) const;

	virtual double DiscountedCPI(
		const ARM_Date& resetDate, 
		const ARM_Date& paymentDate, 
		long dailyInterpType,
		const string& CPILag,
		const string& DCFLag,
		ARM_InfIdx* infIdx = NULL ) =0;

	virtual double FwdCPI(
		const ARM_Date& resetDate, 
		long dailyInterpType,
		const string& CPILag,
		const string& DCFLag,
		ARM_InfIdx* infIdx = NULL );

	virtual double FwdCPI(
		const ARM_Date& resetDate, 
		long dailyInterpType,
		ARM_InfIdx* infIdx = NULL );

	/// a CPI Index ratio is the ratio of two CPI Fixing
	
	virtual ARM_GP_Vector FwdCPIRatio( 
		const ARM_Date& numDate,
		const ARM_Date& denomDate,
		const ARM_Date& paymentDate,
		double multiple,
		double spread,
	 	long dailyInterpType,
		double denomFixing		= GETDEFAULTVALUE,
		ARM_InfIdx* infIdx		= NULL )=0 ;

	virtual ARM_GP_Vector CptFwdCPIRatio( 
		const ARM_Date& numDate,
		const ARM_Date& denomDate,
		const ARM_Date& paymentDate,
		double multiple,
		double spread,
	 	long dailyInterpType,
		double denomFixing		= GETDEFAULTVALUE,
		ARM_InfIdx* infIdx		= NULL ){
			return FwdCPIRatio( numDate,denomDate,	paymentDate, 1.0,-1.0,	dailyInterpType, denomFixing,	infIdx); }

	virtual double		GetCPIIndexValue() const;
	virtual ARM_Date	GetCPIIndexDate()  const{	return	itsInfFwdCurve->GetCPIIndexDate();} 

};




/*!
 * \class	ARM_InfBSModel
 * \author  Eric Benhamou
 * \version 1.0
 * \date August 2003
 *
 * \brief an inflation curve model is composed of
 * a yield curve and a cpi forward curve
 *
 */
class ARM_InfCurvModel : public ARM_Model, public InfFwdModel
{
public:
	ARM_InfCurvModel( ARM_ZeroCurve*	discountCurve, ARM_InfCurv*	infFwdCurv );
	ARM_InfCurvModel( const ARM_InfCurvModel& rhs);
	ARM_InfCurvModel& operator = (const ARM_InfCurvModel &rhs );
	virtual ~ARM_InfCurvModel();

	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	/// pricing of the interest rates part
	virtual double ForwardYield(
		double calcDate, 
		double resetDate,
        double maturityDate, 
		double yieldMaturity, 
        int compMeth = 0, 
		int dayCount = KACTUAL_365, 
        int DomOrFrgRate = 1);

	virtual double ExpectedFwdYield(
		ARM_Date& resetDate, 
        ARM_Date& maturity, 
        ARM_Date& payDate, 
        int compMeth = 0, 
        int dayCount = KACTUAL_365, 
        int DomOrFrgRate = 1,
        int discYC = 1,
        int YieldDecomp = K_COMP_PROP,
        double Margin = 0.0,
        StoreFwdRateInfo* StoreInfo = NULL,
		int IndexFreq = -1);

	/// Pricing functions.. the pure inflation functions
	/// these ones have to be redefined since it is pure virtual!
	
	//TODO : à virer 
	virtual double DiscountedCPI(
		const ARM_Date& resetDate, 
		const ARM_Date& paymentDate, 
		long dailyInterpType,
		const string& CPILag,
		const string& DCFLag,
		ARM_InfIdx* infIdx = NULL); //TODO : à virer 

	virtual ARM_GP_Vector FwdCPIRatio( 
		const ARM_Date& numDate,
		const ARM_Date& denomDate,
		const ARM_Date& paymentDate,
		double multiple,
		double spread,
	 	long dailyInterpType,
		double denomFixing		= GETDEFAULTVALUE,
		ARM_InfIdx* infIdx		= NULL);//TODO : à virer  
	//TODO : à virer 
	virtual ARM_GP_Vector CptFwdCPIRatio( 
		const ARM_Date& numDate,
		const ARM_Date& denomDate,
		const ARM_Date& paymentDate,
	 	long dailyInterpType,
		double denomFixing		= GETDEFAULTVALUE,
		ARM_InfIdx* infIdx		= NULL) { //TODO : à virer 
			return FwdCPIRatio( numDate,denomDate,	paymentDate, 1.0,-1.0,	dailyInterpType, denomFixing,	infIdx); }

	virtual ARM_INF_PRICING_INFO CanPriceInflation() const { return PRICE_FWD_CPI; }
};

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
