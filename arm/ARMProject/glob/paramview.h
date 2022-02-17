/*
 * Copyright CDC IXIS CM : Paris July 2003
 *
 * $Log: paramview.h,v $
 * Revision 1.6  2003/08/06 13:19:41  ekoehler
 * use of dos2unix
 *
 * Revision 1.4  2003/08/01 14:34:22  ekoehler
 * typo correction: a * had gone
 *
 * Revision 1.1  2003/07/30 12:12:50  ekoehler
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*
    paramview.h
 
    ParamView framework

*----------------------------------------------------------------------------*/
#ifndef _PARAMVIEW_H
#define _PARAMVIEW_H

#include "armglob.h"

/*----------------------------------------------------------------------------*/
#define ARM_PARAMVIEW_ENDOFLINE			-1111111111
#define ARM_PARAMVIEW_ENDOFLINE_CHAR	"ENDOFLINE"

#define ARM_PARAMTABLESIZE 200

/*----------------------------------------------------------------------------*/

/* for a given category of constants K_SOMETHING, table of mapping K_SOMETHING vs "Something"*/
typedef struct
{
	long	itsMethodFlag;
	char*	itsMethodName;
} ARM_TagArrayRow;


/* Table of all the mapping tables*/
typedef struct
{
	char*			itsGenericTag;
	ARM_TagArrayRow*	itsTagArray;
} ARM_ParamMappingRow;


class ARM_ParamView : public ARM_Object
{
		static ARM_ParamMappingRow itsParamMappingTable[ARM_PARAMTABLESIZE];
	public:
        ARM_ParamView() {};
		static char* GetMappingName( char* GenericTag, long MethodFlag );
        virtual ~ARM_ParamView() {};
};

#endif
/*---------------------------------------------------------------*/
/*--- End Of File ---*/

