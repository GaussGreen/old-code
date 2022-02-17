/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file prcscalculator.h
 *
 *  \brief file for the Power Reverse Dual Currencies Calculator
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date February 2006
 */


#ifndef _INGPCALCULATORS_PRCSCALCULATOR_H
#define _INGPCALCULATORS_PRCSCALCULATOR_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "prdcalculator.h"

#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

///-----------------------------------------------------------------------------
/// \class ARM_PRDCCalculator
/// \brief
///  Class that implements a Power Reverse Dual Currencies calculator
///-----------------------------------------------------------------------------
class ARM_PRDCCalculator : public ARM_PRDCalculator 
{
private: 
    static const string PRDCColNamesTable [];

public:

    enum PRCSColAlias
    {
        PRDCBermuda=0,
        PRDCOption
    };

	/// constructor for deal from term sheet
	ARM_PRDCCalculator(const ARM_Date& asofDate,
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
		const ARM_Curve& MinCpn,
		const ARM_Curve& MaxCpn,
		const ARM_Curve& initialFx,
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
		const ARM_StringVector& productsToPrice		= ARM_StringVector(1,PRDCColNamesTable[PRDCOption]),
		bool fxLocalModelFlag						= true,
		bool basisIRCalibFlag					    = true,
		const ARM_Curve& fundlevrage                = ARM_FlatCurve(1));

	/// Constructor for trade id from Summit
	ARM_PRDCCalculator(ARM_PowerReverse* powRev,
		ARM_DFBSModel* model,
		const ARM_ObjectVector& otherMktDatas,
		const ARM_GP_Vector& schedulerDatas = ARM_GP_Vector(0),
		const ARM_GP_Vector& truncatorDatas = ARM_GP_Vector(0),
		const ARM_StringVector& columnsToPrice=ARM_StringVector(1,PRDCColNamesTable[PRDCOption]),
		bool markovianDriftSamplerFlag = true,
		bool fxLocalModelFlag = false,
		ARM_PRDCCalibType calibType = ARM_PRCSCalibTypes::ATMCalib, 
		const ARM_GP_Vector& fxATSCalibDatas = ARM_GP_Vector(0),
		bool basisIRCalibFlag = true,
		ARM_BasisType basisType = ARM_PRCSBasisType::flowByflow);

	///copy constructor, assignment constructor, destructor
	ARM_PRDCCalculator( const ARM_PRDCCalculator& rhs );
	ASSIGN_OPERATOR(ARM_PRDCCalculator)
	~ARM_PRDCCalculator(){};


    /// PRDC deal description creation functions
	void CheckData();
	
	void InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const;
	virtual ARM_RowInfo ColumnNames() const;
	virtual ARM_RowInfo MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;


    /// Standard ARM support
    virtual ARM_Object* Clone() const { return new ARM_PRDCCalculator(*this); }

	virtual string ExportShortName() const { return "LPRDC";}
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

