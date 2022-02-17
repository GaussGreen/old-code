/*************************************************************************** 
*  Name:        ddlmaker/include/backend.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: backend.h,v 1.3 2005/07/07 15:54:42 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __BACKEND_H__
#define __BACKEND_H__

#include <boost/shared_ptr.hpp>
#include <string>

class Column;
class Table;
class ConstraintBase;
class ConstraintPrimaryKey;
class ConstraintForeignKey;
class ConstraintUnique;
class ConstraintNull;
class TypeBase;

/*
 * This will be the class implemented in the plugins
 */

class Backend
{
public:
	/* prolog */
	virtual void Begin() = 0;
	
	virtual void HandleTable(Table& table) = 0;

	virtual void HandleColumn(Column& column) = 0;

	virtual void HandlePrimaryKey(ConstraintPrimaryKey& primarKey) = 0;

	virtual void HandleForeignKey(ConstraintForeignKey& foreignKey) = 0;

	virtual void HandleUnique(ConstraintUnique& unique) = 0;

	virtual void HandleNotNull(ConstraintNull& notNull) = 0;

	virtual void HandleType(TypeBase& type) = 0;

	virtual std::string GetOutput() = 0;
	
	/* epilog */
	virtual void End() = 0;
	
	virtual ~Backend()
	{
	}
};

#endif
