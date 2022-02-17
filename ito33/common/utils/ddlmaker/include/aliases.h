/*************************************************************************** 
*  Name:        ddlmaker/include/aliases.h                                 
*  Purpose:	This is a special kind of backend. It is used to generate
*  		an inverse table:
*  		sql_table_name.=SCHEMA_H_TABLE_NAME
*  		sql_table_name.sql_column_name=SCHEMA_H_COLUMN_NAME
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: aliases.h,v 1.1 2005/09/16 15:35:55 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __ALIASES_H__
#define __ALIASES_H__

#include "backend.h"

#include "table.h"
#include "column.h"
#include "constraint.h"

#include <string>
#include <map>
#include <sstream>

/*
 * This class will generate output like this:
 * table_name.=SQL_IDENT(table_name)
 * table_name.column_name=SQL_IDENT(column_name)
 *
 * Of course the output needs to be passed the by 
 * cpp which should define __SQL__ and take 
 * sqlschema.h as macro-definitions file
 */
class AliasesBackend
:public Backend
{
public:
	AliasesBackend()
	{
	}
	
	void Begin() {}
	
	void HandleTable(Table& table);

	void HandleType(TypeBase& type) {}

	void HandleColumn(Column& column);

	void HandlePrimaryKey(ConstraintPrimaryKey& primaryKey) {}
	void HandleForeignKey(ConstraintForeignKey& foreignKey) {}
	void HandleUnique(ConstraintUnique& unique) {}
	void HandleNotNull(ConstraintNull& notNull) {}

	std::string GetOutput();

	/* nothing to do */
	void End()
	{
	}

private:
	/* assemble all output here */
	std::ostringstream m_output;	
	
	/* make it non-copyable */
	AliasesBackend(const AliasesBackend&);
	AliasesBackend& operator=(const AliasesBackend&);
};

#endif
