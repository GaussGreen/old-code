/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: assetinfo.cpp,v $
 * Revision 1.1  2004/06/08 16:44:43  ebenhamou
 * Initial revision
 *
 */

/*! \file assetinfo.cpp
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */

#include "gpinflation/assetinfo.h"
#include "gpbase/ostringstream.h"
#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_SingleAssetInfo
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_SingleAssetInfo::toString( const string& indent ) const
{
	CC_Ostringstream os;
	os << indent << " Type		: ";
	switch ( itsType )
	{
		case K_FIXED_LEG:
			{
			os << "Fixed Leg\n";
			break;
			}
		case K_GENERICINF_LEG :
			os << "Inflation Leg\n";
			break;
		case K_FLOATING_LEG :
			os << "Floating Libor Leg\n";
			break;
		default:
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": unknown type, allowed is inflation, fixed, floating" );
	};

	os << indent << " Name		: " << itsName << "\n";
	os << indent << " Sign		: ";
	
	if( itsRcvOrPay == K_PAY )
		os << "PAY\n";
	else
		os << "REC\n";

	if ( (itsType == K_GENERICINF_LEG || itsType == K_FLOATING_LEG) )
	{
		os << indent    << " Forward Price	: " << itsFwd    << "\n";
		os << indent    << " BS Vol Used	: " << itsVol    << "\n";
		os << indent    << " Tenor       	: " << itsTenor  << "\n";
		os << indent    << " Expiry       	: " << itsExpiry << "\n";
		os << indent    << " StrikeForVol  	: " << itsStrikeForVol << "\n";
	}
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_SingleAssetInfo
///	Routine: toString
///	Returns: 
///	Action : says whether the asset is fixed or not!
////////////////////////////////////////////////////
bool ARM_SingleAssetInfo::isFixedAsset( ) const
{ 
	return itsType == K_FIXED_LEG;
}



////////////////////////////////////////////////////
///	Class  : ARM_TwoAssetsInfo
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_TwoAssetsInfo::toString( const string& indent ) const
{
	CC_Ostringstream os;
	os << indent << " 1) First Leg\n";
	os << indent << itsFirstAsset.toString( indent );
	os << "\n";

	os << indent << " 2) Second Leg\n";
	os << indent << itsSecondAsset.toString( indent );
	os << "\n";


	if( itsFirstAsset.isFixedAsset() || itsSecondAsset.isFixedAsset() )
	{
		os << indent << " 3) One Asset Pricing Info\n";
	}
	else
	{
		os << indent << " 3) Two Assets Pricing Info\n";
	}

	if( !itsFirstAsset.isFixedAsset() && !itsSecondAsset.isFixedAsset() )
		os << indent << " Correlation	: " << itsCorrelation << "\n";
	os << indent <<  " Option Maturity: " << itsOptionMaturity << "\n";
	os << indent <<  " Annuity	: " << itsAnnuity << "\n";

	os << indent << " Strike		: " << itsStrike << "\n";
	os << indent << " Pricing Strike	: " << itsPricingStrike << "\n";

	os << "\n";
	
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_TwoAssetsInfo
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : clone the object
////////////////////////////////////////////////////
ARM_Object* ARM_TwoAssetsInfo::Clone()
{
	return new ARM_TwoAssetsInfo(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_TwoAssetsInfo
///	Routine: View
///	Returns: void
///	Action : prints characteristic of the object 
///				in a file
////////////////////////////////////////////////////
void ARM_TwoAssetsInfo::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
	
	/// first determine that the file is not already opened
    if ( NULL == ficOut )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

	/// use the method to string
	/// and just says the type and what is in it
    fprintf(fOut, "%s", toString().c_str() );

	/// to allow to have nested view
    if ( NULL == ficOut )
       fclose(fOut);
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/