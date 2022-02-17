/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: mktcst.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file mkcst.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPCALCULATORS_MKTDATACST_H
#define _INGPCALCULATORS_MKTDATACST_H

#include "gpbase/port.h"
#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

/// constant for Summit data
const string SUMMIT_DATA_FOLDER			= "P:\\Samba\\arm-import";
const string SUMMIT_FILE_LOCATION_HISTO = "histo\\";
const string SUMMIT_ZC_FILE_LOCATION	= "\\ZC\\";
const string SUMMIT_VOL_FILE_LOCATION	= "\\VOL\\";
const string SUMMIT_SMILE_FILE_LOCATION	= "\\SMILE\\";
const string SUMMIT_FXVOL_FILE_LOCATION	= "\\FXVOL\\";
const string SUMMIT_YLD_FILE_LOCATION	= "\\YLD\\";
const string SUMMIT_FX_FILE_LOCATION	= "\\FX\\";


/// constant for vol type
struct ARM_MarketData
{
	enum VolType
	{
		MKT_ATM_VOL,
		MKT_ATM_AND_SMILE_VOL
	};

	enum VolMktType
	{
		MKT_CAPORCAPLET_VOL,
		MKT_SWAPTION_VOL
	};
};

const string CAPLET_TENORS_TABLE[]		= { "1M", "3M", "6M", "12M" };
const size_t CAPLET_TENORS_TABLE_SIZE	= sizeof(CAPLET_TENORS_TABLE)/sizeof(CAPLET_TENORS_TABLE[0]);
const string CAP_TENORS_TABLE[]			= { "2Y","5Y", "10Y", "30Y" };
const size_t CAP_TENORS_TABLE_SIZE		= sizeof(CAP_TENORS_TABLE)/sizeof(CAP_TENORS_TABLE[0]);
const string SWO_TENORS_TABLE[]			= { "1Y", "2Y", "10Y", "30Y" };
const size_t SWO_TENORS_TABLE_SIZE		= sizeof(SWO_TENORS_TABLE)/sizeof(SWO_TENORS_TABLE[0]);

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
