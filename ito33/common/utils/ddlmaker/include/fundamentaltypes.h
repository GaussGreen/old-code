/*************************************************************************** 
*  Name:        ddlmaker/include/fundamentaltypes.h                        
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  $Id: fundamentaltypes.h,v 1.2 2005/07/06 09:28:30 cosmin Exp $
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __FUNDAMENTALTYPES_H__
#define __FUNDAMENTALTYPES_H__

/* this data types are used in the input script that
 * is sgbd agnostic */
enum BASE_TYPE_ENUM
{
	BASE_TYPE_INTEGER,
	BASE_TYPE_CHARACTER,
	BASE_TYPE_DECIMAL,
	BASE_TYPE_DATE,
	BASE_TYPE_FLOAT,
	BASE_TYPE_BLOB
};
#endif
