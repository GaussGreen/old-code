/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: enumgencalculator.h,v $
 * Revision 1.1  2005/07/23 14:51:19  emezzine
 * Initial revision
 *
 *
 */

/*! \file enumgencalculator.h
 *
 *  \brief 
 *
 *	\author  E.M Ezzine 
 *	\version 1.0
 *	\date July 2005
 */


#ifndef _INGPCALCULATORS_ENUMGENCALCULATOR_H
#define _INGPCALCULATORS_ENUMGENCALCULATOR_H

#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_GenCalculatorCcyType
{
	enum CcyType
		{
            DomCurrency = 0,
		    ForCurrency,
		    FundCurrency,
		    CpnCurrency,
            Unknown,
		};
};

struct ARM_MRSCalibrationType
{
	enum MRSCalibType
		{
            diagstartfwd = 0,
		    stmfirstcolumn,
            Unknown,
		};
};

struct ARM_MRSStrikeCalibrationType
{
	enum MRSStrikeCalibType
		{
            strikeEquivalent = 0,
		    strikeATM,
            strikeKeepStdMoyenessCst,
            Unknown,
		};
};


struct ARM_SigmaCalibrationType
{
	enum SigmaCalibType
		{
            strikeEquivalent = 0,
		    strikeATM,
			StrikeDeltaApproxi,
            Unknown,
		};
};


struct ARM_PRCSRedemptionType
{
	enum PRCSRedemptionType
	{
		standard =0,
		mandatoryRedemption,
		dualOptionRedemption,
	};
};

struct ARM_PRCSBasisType
{
	enum PRCSBasisType
	{
		flowByflow =0,
		average,
	};
};

struct ARM_PRCSCalibTypes
{
	enum PRDCCalibType
    {
        ATMCalib=0,
        ATSFxCalib,
        ATSFxMixedCalib,
        ATSFxProfileCalib,
        ATSFxMoneynessCalib,
        ATSFxShiftedCalib,
        ATSFxMinVolCalib,
        ATMDoubleCalib,
        ATSFxEquivCalib,
		HybridBasketCalib,
		ATSFxBarrierMoneyness,
    };
};


struct ARM_FXVanillaType
{
	enum FXVanillaType
	{
		vanilla=0,
		spread,
		basket,
		digit,
		perf,
		digitspread,
		quotient,
		FXBall,
		FXBallPerf,
	};
};

struct ARM_FXBasketType
{
	enum BasketType
	{
		min = -1,
		max = 1,
	};
};

struct ARM_MixPRDNoticeType
{
	enum NoticeType
	{
		Call = 0,
		KO ,
	};
};


struct ARM_TARNFXPayoffType
{
	enum TARNFXPayoffType
	{
		TARNFX = 0,
		CHOOSERFX,
		INDIANFX,
		PRDKO,
	};
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/