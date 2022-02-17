/*************************************************************************** 
*  Name:        ddlmaker/include/sqlite.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: sqlite.h,v 1.5 2005/08/24 08:21:52 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __SQLITE_H__
#define __SQLITE_H__

#include "backend.h"

#include "table.h"
#include "column.h"
#include "constraint.h"

#include <string>
#include <map>
#include <sstream>

class SqliteBackend
:public Backend
{
public:
	SqliteBackend()
	{
	}
	
	void Begin();
	
	/* called for each table to generate sql statements */
	void HandleTable(Table& table);

	/* 
	 * called for all types in the symbol table, it 
	 * prepares the types for usage.
	 * For example, it will output the code to create
	 * domains.
	 */
	void HandleType(TypeBase& type);

	/*
	 * generates code for a column
	 */
	void HandleColumn(Column& column);

	/*
	 * we need to use the visitor patter to 
	 * generate code for the constraints as we
	 * use them in a polymorphical manner
	 */
	void HandlePrimaryKey(ConstraintPrimaryKey& primaryKey);
	void HandleForeignKey(ConstraintForeignKey& foreignKey);
	void HandleUnique(ConstraintUnique& unique);
	void HandleNotNull(ConstraintNull& notNull);

	std::string GetOutput();

	/* nothing to do */
	void End()
	{
	}

private:
	/* call the appropriate method for a certain type */
	std::string TypeSelector(TypeBase& type);
	
	/* 
	 * these functions output base types recognized 
	 * by the sqlite sgbd engine.
	 * they handle well-known attributes in order to 
	 * accomplish the later.
	 * For example a base type "char" that has the attribute
	 * "varying" set to true will generate the code "VARCHAR"
	 */
	std::string TypeInt(TypeBase& type);
	std::string TypeCharacter(TypeBase& type);
	std::string TypeDecimal(TypeBase& type);
	std::string TypeDate(TypeBase& type);

	/*
	 * float and blob are not recognized by sqlite */
	std::string TypeFloat(TypeBase& type);
	std::string TypeBlob(TypeBase& type);
	
	std::string ConstraintName(ConstraintBase& constraint);

	/* 
	 * for each column from a specific table output
	 * declaration code in this map.
	 * The key is the column name and the value is the column
	 * code.
	 */
	std::map<std::string, std::string> m_columnCode;

	/* 
	 * for those constraints that have more
	 * than one column in theirs primary list (thus table constraints)
	 * put the generation code in this map.
	 * The key is the constraint name and the value is the
	 * constraint code.
	 */
	std::map<std::string, std::string> m_tableConstraintCode;
	
	std::string ForeignKeyOnDelete(ConstraintForeignKey& foreignKey);
	
	/* assemble all output here */
	std::ostringstream m_output;	
	
	/* make it non-copyable */
	SqliteBackend(const SqliteBackend&);
	SqliteBackend& operator=(const SqliteBackend&);
};

#endif
