/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: advanceduserstable.cpp,v $
 * Revision 1.1  2003/10/16 07:52:11  ebenhamou
 * Initial revision
 *
 */


/*! \file advanceduserstable.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpinfra/advanceduserstable.h"

CC_BEGIN_NAMESPACE( ARM )

extern const string GP_User_Level2_Table[]=
{
	/// exotic desk
	"pamiel",
	"fbourcier",
	"tcayrousse",
	"sgasquet",
	"delandier",
	"ypilchen",

	/// London quant
	"nbiala",
	"dkalafatis",
	"spannetier",


	/// Paris Quant
	"nbelgrade",
	"ebenhamou",
	"ocroissant"
	"emezzine",
	"rguillemot",
	"ykhlif",
	"jmprie",
	"arajona",

	/// Paris Lib
	"mabdelmoumni",
	"mcampet",
	"jpriaudel",

	/// other
	"apelletier"

};

extern const size_t GP_User_Level2_Table_Size= sizeof(GP_User_Level2_Table)/sizeof(GP_User_Level2_Table[0]);

bool UserControl_IsInTheList( const string list[], size_t size, const string& elem )
{
	bool isInTheList = false;
	
	for(size_t i=0;i<size; ++i)
		if( elem == list[i] )
		{
			isInTheList =true;
			break;
		}
	return isInTheList;
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

