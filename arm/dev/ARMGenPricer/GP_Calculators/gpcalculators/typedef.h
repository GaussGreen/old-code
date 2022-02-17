/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: typedef.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */

/*! \file typdef.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPCALCULATORS_TYPEDEF_H
#define _INGPCALCULATORS_TYPEDEF_H

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"
#include "enumgencalculator.h"


/// Start ARM namespace
CC_BEGIN_NAMESPACE( ARM )

///struct
struct ARM_GenCalculatorCcyType;
typedef ARM_GenCalculatorCcyType::CcyType ARM_CcyType;

struct ARM_MRSCalibrationType;
typedef ARM_MRSCalibrationType::MRSCalibType ARM_MRSCalibType;

struct ARM_MRSStrikeCalibrationType;
typedef ARM_MRSStrikeCalibrationType::MRSStrikeCalibType  ARM_MRSStrikeCalibType;

struct ARM_SigmaCalibrationType;
typedef ARM_SigmaCalibrationType::SigmaCalibType  ARM_SigmaCalibType;

struct ARM_PRCSRedemptionType;
typedef ARM_PRCSRedemptionType::PRCSRedemptionType  ARM_RedemptionType;

struct ARM_PRCSRedemptionType;
typedef ARM_PRCSRedemptionType::PRCSRedemptionType  ARM_RedemptionType;

struct ARM_PRCSBasisType;
typedef ARM_PRCSBasisType::PRCSBasisType  ARM_BasisType;

struct ARM_PRCSCalibTypes;
typedef ARM_PRCSCalibTypes::PRDCCalibType  ARM_PRDCCalibType;

struct ARM_FXVanillaType;
typedef ARM_FXVanillaType::FXVanillaType  ARM_VanillaType;

struct ARM_FXBasketType;
typedef ARM_FXBasketType::BasketType  ARM_BasketType;

struct ARM_MixPRDNoticeType;
typedef ARM_MixPRDNoticeType::NoticeType  ARM_PRDNoticeType;

struct ARM_TARNFXPayoffType;
typedef ARM_TARNFXPayoffType::TARNFXPayoffType  ARM_FXTARNPayoffType;



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
