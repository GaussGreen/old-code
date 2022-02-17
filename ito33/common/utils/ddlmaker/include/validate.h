/*************************************************************************** 
*  Name:        ddlmaker/include/validate.h                                 
*  Purpose:	Declarations of functions that validate the input.                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-09-05                                                   
*  RCS-ID:      $Id: validate.h,v 1.1 2005/09/12 08:51:50 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __VALIDATE_H__
#define __VALIDATE_H__

#include <list>
#include <string>

typedef std::list<std::string> TKeywordList;

void
LoadKeywordList(const std::string& fileName);

bool
ReservedKeyword(const std::string& name);

#endif // __VALIDATE_H__
