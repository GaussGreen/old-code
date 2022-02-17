//////////////////////////////////////////////////////////////////////
// ICM_CPPI.h: interface for the ICM_CPPI class.
//
//////////////////////////////////////////////////////////////////////

#pragma warning (disable : 4786 )

#if !defined(_ICM_CPPI_H_)
#define _ICM_CPPI_H_

/********************************************************************************/
/*! \file icm_cppi.h
 * 
 *  \brief Describes a Constant Proportion Portfolio Insurance
 *  \author 
 *	\version 1.0
 *	\date   April 2005 */
/*
 *********************************************************************************/


#include "ICMKernel\inst\icm_security.h"
#include "ICMKernel\util\icm_RangeFactor.h"

class ICM_Cppi : public ARM_Security
{
private:
	ARM_Date itsStartDate;
	ARM_Date itsEndDate;
	ARM_Security* itsSecurity;
	ARM_Currency* itsCurrency;
	ICM_RangeFactor* itsRangeFactor;
	double itsNotional;						// en Euros
	double itsProtectedAmount;				// en Euros
	double itsManagementCost;				// en % (ex : 0.5%)
	double itsAdditionalLeverage;			// en % du notional (ex : 1%)
	double itsDesactivateCushion;			// en % du notional (ex : 1.5%)
	
	double itsDefTime;
	double itsAvgLossValue;
	double itsNbLoss;
	double itsAvgFinalValue;

	string itsCorrelName;

public:
	void Init();

	// -----------------------------------
	// Constructor / Destructor
	// -----------------------------------
	ICM_Cppi()	{ Init();}

	ICM_Cppi(const ARM_Date& startDate,
			 const ARM_Date& endDate,
			 ARM_Security* sec,
			 const double& notional,
			 const double& ProtectedNotional,
			 const double& AdditionalLeverage,
			 ICM_RangeFactor* riskyFactor,
			 const double& DesactivateCushion /** = 0.0 **/ ,
			 const double& ManagementCost /** = 0.0**/ ,
			 ARM_Currency* ccy /** = ARM_DEFAULT_CURRENCY**/ ,
			 string CorrelName /** = "TRAXX" **/ )
	{
		Init();
		
		Set(startDate,
			endDate,
			sec,
			notional,
			ProtectedNotional,
			AdditionalLeverage,
			riskyFactor,
			DesactivateCushion,
			ManagementCost,
			ccy,
			CorrelName);
		
	}

	void Set(const ARM_Date& startDate,
			 const ARM_Date& endDate,
			 ARM_Security* sec,
			 const double& notional,
			 const double& ProtectedNotional,
			 const double& AdditionalLeverage,
			 ICM_RangeFactor* riskyFactor,
			 const double& DesactivateCushion/**  = 0.0**/ ,
			 const double& ManagementCost /** = 0.0**/ ,
			 ARM_Currency* ccy /** = ARM_DEFAULT_CURRENCY **/ ,
			 string CorrelName /** = "TRAXX" **/ );
			 

	virtual ~ICM_Cppi() 
	{	
		if (itsSecurity)
			delete itsSecurity;
		itsSecurity = NULL;

		if (itsCurrency)
			delete itsCurrency;
		itsCurrency = NULL;

		if (itsRangeFactor)
			delete itsRangeFactor;
		itsRangeFactor = NULL;
	}

	// -----------------------------------
	// Copy and Clone Methode
	// -----------------------------------
	void BitwiseCopy(const ARM_Object* src);

	void Copy(const ARM_Object* src);

	ARM_Object* Clone(void);

	// -----------------------------------
	// Get ET Set
	// -----------------------------------
	ARM_Date GetStartDate()				{ return itsStartDate;}
	ARM_Date GetEndDate()				{ return itsEndDate;}
	ARM_Security* GetSecurity()			{ return itsSecurity;}
	ARM_Currency* GetCurency()			{ return itsCurrency;}
	ICM_RangeFactor* GetRangeFactor()	{ return itsRangeFactor;}
	double GetNotional()				{ return itsNotional;}
	double GetProtectedAmount()	 		{ return itsProtectedAmount;}
	double GetManagementCost()			{ return itsManagementCost;}
	double GetAdditionalLeverage()		{ return itsAdditionalLeverage;}
	double GetDesactivateCushion()		{ return itsDesactivateCushion;}
	double GetDefTime()					{ return itsDefTime;}
	double GetAvgLossValue()			{ return itsAvgLossValue;}
	double GetNbLoss()					{ return itsNbLoss;}
	double GetAvgFinalValue()			{ return itsAvgFinalValue;}
	string GetCorrelName()				{ return itsCorrelName;}
	


	void SetStartDate(const ARM_Date& startdate)	{itsStartDate = startdate;}
	void SetEndDate(const ARM_Date& enddate)		{itsEndDate = enddate;}
	void SetCurrency(ARM_Currency*  ccy)			{itsCurrency = ccy;}
	void SetRangeFactor(ICM_RangeFactor* Rfactor)	{itsRangeFactor = Rfactor;}
	void SetNotional(const double& value)			{itsNotional = value;}
	void SetProtectedAmount(const double& value)	{itsProtectedAmount = value;}
	void SetManagementCost(const double& value)		{itsManagementCost = value;}
	void SetAdditionalLeverage(const double& value)	{itsAdditionalLeverage = value;}
	void SetDesactivateCushion(const double& value)	{itsDesactivateCushion = value;}
	void SetDefTime(const double& value)			{itsDefTime = value;}
	void SetAvgLossValue(const double& value)		{itsAvgLossValue = value;}
	void SetNbLoss(const double& value)				{itsNbLoss = value;}
	void SetAvgFinalValue(const double& value)		{itsAvgFinalValue = value;}
	void SetCorrelName(const string& value)			{itsCorrelName = value;}

	// --------------
	//	View Method
	// --------------
	void View(char* id, FILE* ficOut);
	
};
#endif