/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file binaryprdccalculator.h
 *
 *  \brief file for the Power Reverse Dual Currencies Calculator
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date May 2007
 */


#ifndef _INGPCALCULATORS_BINARYPRCSCALCULATOR_H
#define _INGPCALCULATORS_BINARYPRCSCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "prdcalculator.h"

#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

///-----------------------------------------------------------------------------
/// \class ARM_BinaryPRDCCalculator
/// \brief
///  Class that implements a Power Reverse Dual Currencies calculator
///-----------------------------------------------------------------------------
class ARM_BinaryPRDCCalculator : public ARM_PRDCalculator 
{
private: 
    static const string BiPRDCColNamesTable [];
	static const string BiPRDCProfileNamesTable [];

	double itsEpsilon;
	ARM_DigitType  itsDigitType;
	ARM_Curve itsBarrier;
	ARM_Curve itsDigitFixCpn;

	ARM_GP_Vector itsvBarrier;
	ARM_GP_Vector itsvDigitLvge;
	ARM_BoolVector	itsvCpnIsDigit;

public:

	enum BiPRDCProfileType
    {
        barrierUpProfile = 0,
		barrierDownProfile,
		digitLvgeProfile,
	};

    enum BiPRCSColAlias
    {
		FxDownCallStrip =0,
		FxUpCallStrip,
		DigitalStrip,
		PRDCBermuda,
        PRDCOption
    };

	/// constructor for deal from term sheet
	ARM_BinaryPRDCCalculator(const ARM_Date& asofDate,
		const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& DomCcy,
		const ARM_Currency& ForCcy,
		const ARM_Currency& FundCcy,
		int cpnDayCount,
		int cpnFreq,
		int FxResetGap,
		int stubRule,
		int resetTiming,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& domesticCpn,
		const ARM_Curve& foreignCpn,
		const ARM_Curve& minCpn,
		const ARM_Curve& barrier,
		const ARM_Curve& initialFx,
		const ARM_Curve& fixCpn,
		int fundFreq,
		int fundDayCount,
		int compFreq,
		int compType,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,		
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNCall,
		const ARM_Curve& fees,
		int redemptionGap ,
		double redemptionStrike,
		const ARM_RedemptionType& redemptionType    = ARM_PRCSRedemptionType::standard,		
		const ARM_StringVector& productsToPrice		= ARM_StringVector(1,BiPRDCColNamesTable[PRDCOption]),
		bool fxLocalModelFlag						= true,
		bool basisIRCalibFlag					    = true,
		const ARM_Curve& fundlevrage                = ARM_FlatCurve(1),
		ARM_DigitType digitType                     = ARM_FXDigitType::centred,
		double epsilon                              = 1.0e-6);

	///copy constructor, assignment constructor, destructor
	ARM_BinaryPRDCCalculator( const ARM_BinaryPRDCCalculator& rhs );
	ASSIGN_OPERATOR(ARM_BinaryPRDCCalculator)
	~ARM_BinaryPRDCCalculator(){};

    /// PRDC deal description creation functions
	void CheckData();

	ARM_CstManagerPtr ComputeCstManager();

	/// Convert ARM_Curves into relevant vectors
	void ComputeProductVectorsFromCurves();
	
	void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;
	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;


    /// Standard ARM support
    virtual ARM_Object* Clone() const { return new ARM_BinaryPRDCCalculator(*this); }

	virtual string ExportShortName() const { return "LPRDC";}
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

