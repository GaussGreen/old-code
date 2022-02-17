/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file CallableQuantoCalculator.h
 *
 *  \brief
 *
 *	\author  J. Messines
 *	\version 1.0
 *	\date apr 2007
 */

#ifndef _INGPCALCULATORS_CALLABLEQUANTOCALCULATOR_H
#define _INGPCALCULATORS_CALLABLEQUANTOCALCULATOR_H

#include "gencalculator.h"


CC_BEGIN_NAMESPACE( ARM )

const string LOCAL_MODEL_NAME	= "LOCAL";

class ARM_CallableQuantoCalculator : public ARM_GenCalculator
{
protected:

	static const string CallableQuantoColNamesTable [];
	static const int CallableQuantoProductToPriceColumns [];

public:
	
	enum CallableQuantoColAlias
	{
		ResetDate = 0,
		StartDate,
		EndDate,
		Fees,
		Bermuda,
		Option,

		NbCols,
	};

	enum CallableQuantoProductToPriceAlias
	{
		OptionPrice,
		
		NbProductsToPrice,
	};

	enum mdmKeysAlias
    {
        YcDomKey=0,
        YcForKey,
		ForexKey,
		YcBasisDomKey,
        YcBasisForKey,
		OswDomModelKey,
        OswForModelKey,
		CfDomModelKey,
        CfForModelKey,
        CorrelMatrixKey,
        MrsDomKey,
        MrsForKey,
		FXVolKey,
        
        NbKeys
    };

	enum ScheduleCount
	{
		ExerSched = 0,
		
		NbSched,
	};
    

public:
	ARM_CallableQuantoCalculator(const ARM_CallableQuantoCalculator&);
	
	ARM_CallableQuantoCalculator& operator=(const ARM_CallableQuantoCalculator&);

	ARM_CallableQuantoCalculator(	const ARM_Currency& ccyDom,
									const ARM_Currency& ccyFor,
									ARM_Date&	startDate,
									ARM_Date&	endDate,
									int payRec,
									const ARM_ReferenceValue&  notional,
									int callFreq,
									int callNotice,
									const string& callCal,
									const ARM_ReferenceValue&  callFees,
									int fundFreq,
									int fundDayCount,
									const ARM_ReferenceValue& fundSpread,
									int cpnPayFreq,
									int cpnDayCount,
									const string& cpnResetCal,
									const string& cpnPayCal,
									const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice
								);

	~ARM_CallableQuantoCalculator();

	void CleanUp();
    void CopyNoCleanUp(const ARM_CallableQuantoCalculator& rhs);

	/// ----------------------------------------------------------
	/// ------------- pure virtual functions

	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_RowInfo MiddleRows( size_t i, const ARM_DateStripCombiner& datesStructure ) const;
	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;

	virtual void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;

	/// pricing function
	virtual void ComputePricingData() const;
	virtual void CreateAndSetModel();
	virtual void CreateAndSetCalibration() {};							
    virtual void UpdateModel();
    virtual void UpdateCalibration(bool isUpdateStrike=true) {};		
	virtual void Calibrate() {};										
	virtual double Price();
    virtual void CheckData();
	virtual void CheckMktData();

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  { return ARM_DateStripPtr(NULL);};
	inline virtual ARM_DateStripPtr GetRefFundDateStrip() const { return ARM_DateStripPtr(NULL);};
	inline virtual ARM_DateStripPtr GetRefDateStrip() const { return ARM_DateStripPtr(NULL);};

	/// Vector
	inline virtual ARM_GP_Vector GetvCpnNominal() const  { return ARM_GP_Vector(0);};
	inline virtual ARM_GP_Vector GetvFundNominal() const { return ARM_GP_Vector(0);};
	inline virtual ARM_GP_Vector GetvFundSpread() const { return ARM_GP_Vector(0);};

	/// the discount curve
	inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const  { return NULL;};
	inline virtual ARM_ZeroCurve* GetForeignZeroCurve() const  { return NULL;};
	inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return NULL;};
	inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve() const  { return NULL;};
	inline virtual ARM_Forex* GetForex() const  { return NULL;};


	/// -------------- end of pure virtual functions

	///	to do 
	virtual size_t	ExtraDescSize() const = 0;
	virtual void	UpdateExtraDesc(size_t i, const ARM_DateStripCombiner& datesStructure, vector<string>&, vector<ARM_GP_VALUE_TYPE>&) const = 0;
	virtual void	UpdateExtraColNames(vector<string>&) const = 0;
	virtual size_t	ExtraNbProductsToPrice() const = 0;
	virtual string	ColNames(size_t i) const = 0;

protected:

	// GENERAL DATAS
	ARM_Date					itsStartDate;           
	ARM_Date					itsEndDate;
	int							itsPayRec;
	ARM_ReferenceValue			itsNotional;

	// CALL DATAS
	int							itsCallFreq;
	int							itsCallNotice; // Positive value
	string						itsCallCal;
	ARM_ReferenceValue			itsCallFees;   // NoticeDates, Associated Fees 

	// FUNDING LEG / Swap Var leg
	int							itsFundFreq;
	ARM_ReferenceValue			itsFundSpread;
	int							itsFundDayCount; //ACT360,ACT365
	
	// COUPON LEG  / Swap fix leg
	int							itsCpnPayFreq;
	int							itsCpnDayCount;
	string						itsCpnResetCal; //Used in case of variable payment: Fixing calendar
	string                      itsCpnPayCal;	//Payment calendar

	// Pricing Data
    std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/				itsProductsToPrice;
	vector<double>				itsProductPrices;
	bool						itsHasBeenPriced;
	
	ARM_GP_Vector       itsSchedulerDatas;
    ARM_GP_Vector       itsTruncatorDatas;

	inline int	GetFundFreq()		const {return itsFundFreq;}
	inline int	GetFundDayCount()	const {return itsFundDayCount;}
	inline int	GetPayRec()			const {return itsPayRec;}
	inline int	GetCpnPayFreq()		const {return itsCpnPayFreq;}
	inline int	GetCpnDayCount()	const {return itsCpnDayCount;}

	inline const ARM_ReferenceValue& GetFundSpread()	const {return itsFundSpread;}
	inline const ARM_ReferenceValue& GetNotional()		const {return itsNotional;}

	inline void						SetProductsToPrice(const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice){itsProductsToPrice = productsToPrice;};
	inline const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/&		GetProductsToPrice() const	{return itsProductsToPrice;}

	inline const ARM_Date&			GetStartDate() const	{return itsStartDate;};
	inline const ARM_Date&			GetEndDate() const	{return itsEndDate;};
	
};

CC_END_NAMESPACE()

#endif