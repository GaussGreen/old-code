/*Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file CraQuantoCalculator.h
 *
 *  \brief
 *
 *	\author  J. Messines
 *	\version 1.0
 *	\date apr 2007
 */

#ifndef _INGPCALCULATORS_CRAQUANTOCALCULATOR_H
#define _INGPCALCULATORS_CRAQUANTOCALCULATOR_H

#include "callablequantocalculator.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_CraQuantoCalculator : public ARM_CallableQuantoCalculator
{
protected:

	static const string CraQuantoColNamesTable [];
	static const int CraQuantoProductToPriceColumns [];

public:
	
	enum CraQuantoColAlias
	{
		Funding = 0,
		Corridor,
		Swap,
		
		NbCols,
	};

	enum CraQuantoProductToPriceAlias
	{
		FundingPrice,
		CorridorPrice,
		SwapPrice,
		
		NbProductsToPrice,
	};

public:
	ARM_CraQuantoCalculator(const ARM_CraQuantoCalculator&);
	
	ARM_CraQuantoCalculator& operator=(const ARM_CraQuantoCalculator&);

	ARM_CraQuantoCalculator(	const ARM_Currency& ccyDom,
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
								int payIndex,
								int payResetTiming,
								const ARM_ReferenceValue& fixRate,
								const ARM_ReferenceValue& barrierDown,
								const ARM_ReferenceValue& barrierUp,
								int refResetFreq,
								int refResetTiming,
								int refIndex,
								const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice);

	~ARM_CraQuantoCalculator();

	virtual ARM_Object*			Clone() const;
	
	void CleanUp();
    void CopyNoCleanUp(const ARM_CraQuantoCalculator& rhs);

	virtual ARM_CstManagerPtr	CreateCstManager();
	virtual ARM_StringVector	PricedColumnNames() const;
	
	// to be redefined for each callable quanto
	virtual size_t				ExtraDescSize() const;
	virtual void				UpdateExtraDesc(size_t i, const ARM_DateStripCombiner& datesStructure, vector<string>&, vector<ARM_GP_VALUE_TYPE>&) const;
	virtual void				UpdateExtraColNames(vector<string>&) const;
	virtual size_t				ExtraNbProductsToPrice() const;
	virtual string				ColNames(size_t i) const;

	virtual void				CreateAndSetCalibration();
	virtual void				Calibrate();
	void						CreateEmptyCalibration();
	ARM_StdPortfolioPtr			CreateSwaptionPortfolio(ARM_Currency& ccy);

	void						UpdateCalibration(bool isUpdateStrike=true);
	void						UpdateSwaptionPrices();
	void						ComputeSwaptionPrices(ARM_StdPortfolioPtr& pf, ARM_BSModel* bsModel);

private:

	int					itsPayIndex;
	int					itsPayResetTiming;   
	ARM_ReferenceValue	itsFixRate;
	ARM_ReferenceValue	itsPayIndexMult;

	int					itsRefIndex;
	int					itsRefResetFreq;   
	int					itsRefResetTiming;
	ARM_ReferenceValue	itsBarrierDown;
	ARM_ReferenceValue	itsBarrierUp;

	ARM_ReferenceValue	itsFundSpread;
	ARM_ReferenceValue	itsNotional;

	inline int GetPayIndex()							const {return itsPayIndex;}
	inline int GetPayResetTiming()						const {return itsPayResetTiming;}
	inline int GetRefIndex()							const {return itsRefIndex;}
	inline int GetRefResetTiming()						const {return itsRefResetTiming;}
	inline int GetRefResetFreq()						const {return itsRefResetFreq;}

	inline const ARM_ReferenceValue& GetBarrierDown()	const {return itsBarrierDown;}
	inline const ARM_ReferenceValue& GetBarrierUp()		const {return itsBarrierUp;}
	inline const ARM_ReferenceValue& GetFixRate()		const {return itsFixRate;}
	inline const ARM_ReferenceValue& GetPayIndexMult()	const {return itsPayIndexMult;}

};

CC_END_NAMESPACE()

#endif