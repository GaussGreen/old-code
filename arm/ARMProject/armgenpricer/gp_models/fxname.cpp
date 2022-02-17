/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file forexvanilla.h
 *  \brief to create an object for Fx Name
 * 
 *	\author  K. Belkheir
 *  \ reviwer E.Ezzine
 *	\version 1.0
 *	\date January 2007
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/fxname.h"

CC_BEGIN_NAMESPACE( ARM )

///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
/// --- Class ARM_FXName ---
///$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


///  Tables of the FX

const string ARM_FXName::MktNamesTable [] =
{
	"EURJPY",
	"EURUSD", 	
	"EURCHF",
	"EURAUD",
	"EURTRY",
	"EURGBP",
	"EURMXN",
	"EURSEK",
	"USDJPY",
	"USDAUD", 
	"USDCHF",
	"USDBRL",
	"USDMXM",
	"USDSEK",
	"GBPUSD",
	"AUDJPY",
};

//--------------------------------------------------------------------------------------------------------------
/// standard constructor

ARM_FXName::ARM_FXName( const ARM_FXName& rhs )
	:
ARM_RootObject(rhs),
	itsIsInvMkt(rhs.itsIsInvMkt),
	itsMktName(rhs.itsMktName),
	itsMktNamesVect(rhs.itsMktNamesVect),
	itsInvMktNamesVect(rhs.itsInvMktNamesVect)	
{
}

ARM_FXName::ARM_FXName(const string& fxName)
:
	itsIsInvMkt(false),
	itsMktName(fxName),
	itsMktNamesVect(),
	itsInvMktNamesVect()
{
	/// creation of the itsMktNamesVect and itsInvMktNamesVect
	size_t MktNamesTableSize = sizeof(MktNamesTable)/sizeof(MktNamesTable[0]);
	for(int i=0;i<MktNamesTableSize; i++)
	{
		string name(MktNamesTable[i]);
		itsMktNamesVect.push_back(name);
		itsInvMktNamesVect.push_back(name.substr(3,3) + name.substr(0,3));
	}
	ARM_GP_StrVector:: const_iterator Iter = itsMktNamesVect.find(fxName);
	if(!itsMktNamesVect.contain(fxName))
	{
		if(itsInvMktNamesVect.contain(fxName)){
			itsMktName = fxName.substr(3,3)+ fxName.substr(0,3);
			itsIsInvMkt=true;
		}
		else 
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+": The FX " +  fxName + " does not exist in the market");
	}
}

//--------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------
//return true if FX isinverted false otherwise
bool ARM_FXName::IsInvMkt(const string& fxName) 
{
	return itsInvMktNamesVect.contain(fxName);
}

CC_END_NAMESPACE()

//-----------------------------------------------------------------------------
/*---- End of file ----*/



