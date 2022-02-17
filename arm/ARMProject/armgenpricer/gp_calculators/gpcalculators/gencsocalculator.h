/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *      \file gencsocalculator.h
 *
 *  \brief
 *
 *	\author  JP Riaudel
 *	\version 1.0
 *	\date July 2005
 */

#ifndef _INGPCALCULATORS_GENCSOCALCULATOE_H
#define _INGPCALCULATORS_GENCSOCALCULATOE_H

#include "gpbase/port.h"
#include "gencalculator.h"

class ARM_SwapLeg;
class ARM_SpreadOption;

CC_BEGIN_NAMESPACE( ARM )

/// Default MDM key names
const string YC_KEY_NAME            = "YC_";
const string OSWMODEL_KEY_NAME      = "OSWMOD_";
const string OSWSMIMODEL_KEY_NAME   = "OSWSMIMOD_";
const string CFMODEL_KEY_NAME       = "CFMOD_";
const string SOMODEL_KEY_NAME       = "SOMOD_";
const string MRS_KEY_NAME           = "MRS_";
const string CORREL_KEY_NAME        = "CORREL_";
const string VOLRATIO_KEY_NAME      = "VOLRATIO_";
const string MRSSPREAD_KEY_NAME     = "MRSSPREAD_";
const string YC_BASIS_KEY_NAME      = "YC_BASIS_";
const string FOREX_KEY_NAME         = "SPOT_FX_";
const string UNKNOWN_KEY_NAME       = "UNKNOWN";

const double NON_CALL_FEE			= 1.0e15;


class ARM_GenCSOCalculator : public ARM_GenCalculator
{

public:

	enum mdmKeysAlias
    {
        YcKey=0,
        MrsKey,
		CorrelKey,
        SoModelKey,
        CfModelKey,
        OswModelKey,
		VolRatioKey,
		MrsSpreadKey,
        FundingKey,
        BasisKey,
		fundBasisKey,
        ForexKey,

        NbKeys
    };

	/// constructor, copy constructor, assignment constructor, destructor
	ARM_GenCSOCalculator(const ARM_Date& startDate,
						 const ARM_Date& endDate,
						 int CMSLong,
						 int CMSShort,
						 int cpnDayCount,
						 int cpnFreq,
						 int cpnResetTiming,
						 const ARM_Curve& cpnnominal,
						 const ARM_Curve& fixCoupon,
						 const ARM_Curve& leverageLong,
						 const ARM_Curve& leverageShort,
						 const ARM_Curve& cpnMin,
						 const ARM_Curve& cpnMax,
						 const ARM_Curve& strike,
						 int fundFreq,
						 int fundDayCount,
						 const ARM_Curve& fundSpread,
						 const ARM_Curve& fundLeverage,
						 int exerciseFreq,
						 int noticeGap,
						 int payRec,
						 const ARM_Curve& fees,
						 std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
						 vector<double> ModelDatas,
						 const ARM_MarketData_ManagerRep& mktDataManager,
						 const ARM_StringVector& mdmKeys);

    // Contructeur used With InitFrom...
	ARM_GenCSOCalculator(const ARM_Date& asOfDate,
						 const ARM_Date& startDate,
                         const ARM_Date& fixEndDate,
						 const ARM_Date& endDate,
						 int CMSLong,
						 int CMSShort,
						 int cpnDayCount,
						 int cpnFreq,
						 int cpnResetTiming,
						 const ARM_Curve& cpnnominal,
						 const ARM_Curve& fixCoupon,
						 const ARM_Curve& leverageLong,
						 const ARM_Curve& leverageShort,
						 const ARM_Curve& cpnMin,
						 const ARM_Curve& cpnMax,
						 const ARM_Curve& strike,
						 int fundFreq,
						 int fundDayCount,
						 const ARM_Curve& fundSpread,
						 const ARM_Curve& fundLeverage,
						 int exerciseFreq,
						 int noticeGap,
						 int payRec,
						 const ARM_Curve& fees);

	/// For CSO bi currencies, it is ugly.....
	ARM_GenCSOCalculator(const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& CpnCcy,
		const ARM_Currency& FundCcy,
		int CMSLong,
		int CMSShort,
		int cpnDayCount,
		int cpnFreq,
		int cpnResetTiming,
		int stubRule,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& leverageLong,
		const ARM_Curve& leverageShort,
		const ARM_Curve& cpnMin,
		const ARM_Curve& cpnMax,
		const ARM_Curve& strikes,
		int fundFreq,
		int fundDayCount,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,
		const ARM_Curve& fundLeverage,
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNCall,
		const ARM_Curve& fees,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
		vector<double> ModelDatas,
		const ARM_MarketData_ManagerRep& mktDataManager);

	ARM_GenCSOCalculator(const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& CpnCcy,
		const ARM_Currency& FundCcy,
		int CMSLong,
		int CMSShort,
		int cpnDayCount,
		int cpnFreq,
		int cpnResetTiming,
		int stubRule,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& leverageLong,
		const ARM_Curve& leverageShort,
		const ARM_Curve& cpnMin,
		const ARM_Curve& cpnMax,
		const ARM_Curve& strikes,
		int fundFreq,
		int fundDayCount,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,
		const ARM_Curve& fundLeverage,
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNCall,
		const ARM_Curve& fees);

	ARM_GenCSOCalculator(const ARM_GenCSOCalculator& rhs);
	ARM_GenCSOCalculator& operator=(const ARM_GenCSOCalculator&);
	virtual ~ARM_GenCSOCalculator(); 

	/// Create the underlying and MktModel
	void CreateUnderlying();

    /// Set the keys or name alias of model (usefull in multi-currency context)
    void SetModelKeys();
	virtual ARM_Vector* ComputeAll(); /// Underlying computing

	virtual ARM_DateStripCombiner DatesStructure() const;
	virtual ARM_DateStripCombiner CustomDatesStructure(const ARM_DateStripVector& dateStrips = ARM_DateStripVector(0)) const;
	inline virtual ARM_DateStripPtr GetFundDateStrip() const { return itsFundDateStrip;};
	inline virtual ARM_DateStripPtr GetStructDateStrip() const { return itsStructDateStrip;};

	/// Specialised version for datas consistency
    //virtual void CheckData() = 0;
	//virtual void CheckMktData() = 0;

	/// Dates Strip
	inline virtual ARM_DateStripPtr GetOriginalFundDateStrip() const  { return itsFundDateStrip;};
	inline virtual ARM_DateStripPtr GetRefFundDateStrip() const { return itsFundDateStrip;};
	inline virtual ARM_DateStripPtr GetRefDateStrip() const { return itsFundDateStrip;};

	/// Vector
	inline virtual std::vector<double> GetvCpnNominal() const  { return itsCpnNominal.GetOrdinates().GetValues();};
	inline virtual std::vector<double> GetvFundNominal() const { return itsFundNominal.GetOrdinates().GetValues();};
	inline virtual std::vector<double> GetvFundSpread() const { return itsFundSpread.GetOrdinates().GetValues();};

	/// the discount curve
	inline virtual ARM_ZeroCurve* GetDomesticZeroCurve() const			{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));			};
	inline virtual ARM_ZeroCurve* GetForeignZeroCurve()	const			{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));     };
	inline virtual ARM_ZeroCurve* GetDomesticDiscountZeroCurve() const  { return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));		};
	inline virtual ARM_ZeroCurve* GetForeignDiscountZeroCurve()	const	{ return NULL;}//dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[fundBasisKey]));   };
	inline virtual ARM_Forex* GetForex() const							{ return dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));			};

private:
	/// CSO financial datas
	ARM_Date                    itsStartDate;           // start date
	ARM_Date                    itsFixEndDate;
	ARM_Date                    itsEndDate;             // end date

	int							itsCMSLong;
	int							itsCMSShort;

    ARM_Curve					itsCpnNominal;          // cpn nominal curve
    ARM_Curve					itsFundNominal;         // funding nominal curve
    ARM_Curve					itsFixCpn;              // fix coupon rate or reverse coupon strike cst...

	int                         itsCpnDaycount;         // reverse coupon day count
    int                         itsCpnFreq;             // reset & pay frequency (f
	int							itsCpnResetTiming;
    string                      itsCpnResetCal;         // reverse reset calendar
    string                      itsCpnPayCal;           // reverse payment calendar
	int                         itsStubRule;            // ability to have shortstart .... 

    ARM_Curve					itsCpnFloor;
    ARM_Curve					itsCpnCap;
    ARM_Curve					itsLeverageShort;
    ARM_Curve					itsLeverageLong;
	ARM_Curve					itsStrike;

	int							itsFundFreq;
	int							itsFundDaycount;
    ARM_Curve					itsFundSpread;
	ARM_Curve					itsFundLeverage;

	int							itsExerciseFreq;
	int							itsNoticeGap;
	int							itsPayRec;
	int                         itsNbNoCall;

    ARM_Curve					itsFees;

	int							itsNbFixFlows;
    
    /// Names of model used in basis case (chosen among MdM keys)
    mdmKeysAlias itsCpnModelKey;
    mdmKeysAlias itsFundingModelKey;
    mdmKeysAlias itsBasisRefModelKey;

	ARM_SwapLeg*  itsFundingLeg;
	ARM_SwapLeg*  itsFixedLeg;
	ARM_SwapLeg*  itsInitFixedLeg;
	ARM_SpreadOption*   itsFloor;
	ARM_SpreadOption*   itsCap;

protected:

	ARM_DateStripPtr itsFundDateStrip;
	ARM_DateStripPtr itsStructDateStrip;
	ARM_DateStripPtr itsCalibDateStrip;
	ARM_DateStripPtr itsExerciseDateStrip;
	ARM_DateStripPtr itsExerciseDateStripUnadj;
	
    vector<double> itsModelDatas;

	/// Flag to specify product to price
    std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ itsProductsToPrice;
	bool itsHasBeenPriced;

	// Pricing Data
	double itsCSOPrice;
	double itsStructPrice;

	/// To memorise the line in the deal description corresponding to the first call
    /// Mutable because may be updated in const functions
    CC_IS_MUTABLE size_t itsFirstEventIdx;

    /// Accessors
	inline const ARM_Curve&		GetFixCpn()			const { return itsFixCpn;};
	inline const ARM_Curve&		GetCpnCap()			const { return itsCpnCap;};
	inline const ARM_Curve&		GetCpnFloor()		const { return itsCpnFloor;};
	inline const ARM_Curve&		GetCpnNominal()		const { return itsCpnNominal;};
    inline const ARM_Curve&		GetFundNominal()	const { return itsFundNominal;};
	inline const ARM_Curve&		GetLeverageLong()	const { return itsLeverageLong;};
	inline const ARM_Curve&		GetLeverageShort()	const { return itsLeverageShort;};
	inline const ARM_Curve&		GetFundSpread()		const { return itsFundSpread;};
	inline const ARM_Curve&		GetFundLeverage()	const { return itsFundLeverage;};
	inline const ARM_Curve&		GetFees()			const { return itsFees;};
	inline const ARM_Curve&		GetStrike()			const { return itsStrike;};

	inline void	 SetFixCpn(const ARM_Curve& fixCoupon) { itsFixCpn = fixCoupon;};
	inline void	 SetFundNominal(const ARM_Curve& fundNominal) { itsFundNominal = fundNominal;};

    ///Mkt keys alias
    inline mdmKeysAlias GetCpnModelKey() const      {return  itsCpnModelKey;}
    inline mdmKeysAlias GetFundingModelKey() const  {return  itsFundingModelKey;}
    inline mdmKeysAlias GetBasisRefModelKey() const {return  itsBasisRefModelKey;}

    /// Dates
	inline ARM_DateStripPtr GetCalibDateStrip()		const  {return itsCalibDateStrip;	}
	inline ARM_DateStripPtr GetExerciseDateStrip()	const  {return itsExerciseDateStrip;}
	inline ARM_DateStripPtr GetExerciseDateStripUnadj()	const  {return itsExerciseDateStripUnadj;}

	inline int GetCMSLong()	const       { return itsCMSLong;		}
	inline int GetCMSShort() const      { return itsCMSShort;		}
	inline int GetFundFreq() const      { return itsFundFreq;		}
	inline int GetCpnFreq()	const       { return itsCpnFreq;		}
	inline int GetCpnResetTiming() const{ return itsCpnResetTiming;	}
	inline int GetFundDaycount() const  { return itsFundDaycount;	}
	inline int GetCpnDaycount()	const   { return itsCpnDaycount;	}
	inline int GetExerciseFreq() const  { return itsExerciseFreq;	}
	inline int GetNoticeGap() const     { return itsNoticeGap;		}
	inline int GetPayRec() const		{ return itsPayRec;			}
	

	inline ARM_Date GetStartDate() const	{ return itsStartDate;	}
	inline ARM_Date GetEndDate() const		{ return itsEndDate;	}
    inline ARM_Date GetFixEndDate() const	{ return itsFixEndDate;	}

	inline int GetNbFixFlows() const { return itsNbFixFlows;};
	inline void SetNbFixFlows (int i) { itsNbFixFlows = i;};

	///
	inline int GetNbNoCall() const { return itsNbNoCall;}

    const ARM_INDEX_TYPE const GetIndexType()  { return GetCurrencyUnit()->GetVanillaIndexType();}

	inline string GetPayCal()  const {return itsCpnPayCal;}
	inline string GetResetCal()  const {return itsCpnResetCal;}

};

CC_END_NAMESPACE()

#endif //_INGPCALCULATORS_GENCSOCALCULATOE_H