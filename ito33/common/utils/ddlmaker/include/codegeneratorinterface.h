/*************************************************************************** 
*  Name:        ddlmaker/include/codegeneratorinterface.h                       
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  $Id: codegeneratorinterface.h,v 1.2 2005/07/06 09:28:30 cosmin Exp $
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __CODEGENERATORINTERFACE_H__
#define __CODEGENERATORINTERFACE_H__

#include "backend.h"

/*
 * Every class in our intermediate data hierarchy that holds
 * key information about the database structure should inherit
 * this class, implement the CodeGenerator method and call the apropriate
 * method of Backend.
 * 
 * This is just a classic application of the Visitor Design Pattern
 */

class CodeGeneratorInterface
{
public:
	virtual void CodeGenerator(Backend& backend) = 0;
};

#endif
