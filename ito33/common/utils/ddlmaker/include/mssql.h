/*************************************************************************** 
*  Name:        ddlmaker/include/mssql.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  $Id: mssql.h,v 1.6 2005/08/24 08:21:52 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __MSSQL_H__
#define __MSSQL_H__

#include "backend.h"

#include "table.h"
#include "column.h"
#include "constraint.h"

#include <string>
#include <map>
#include <sstream>

class MssqlBackend
:public Backend
{
public:
	MssqlBackend()
	{
	}
	
	/* nothing to do */
	void Begin()
	{
	}
	
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

	void End();
private:
	/* handle domain creation for fundamental types */
	std::string DomainDeclGen(TypeBase& type);

	/* call the appropriate method for a certain type */
	std::string TypeSelector(TypeBase& type);
	
	/* 
	 * these functions output base types recognized 
	 * by the mssql sgbd engine.
	 * they handle well-known attributes in order to 
	 * accomplish the later.
	 * For example a base type "char" that has the attribute
	 * "varying" set to true will generate the code "VARCHAR"
	 */
	std::string TypeInt(TypeBase& type);
	std::string TypeCharacter(TypeBase& type);
	std::string TypeDecimal(TypeBase& type);
	std::string TypeDate(TypeBase& type);
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
	
	/* 
	 * output type definitions (domains in pgsql case)
	 * only for the types that are actually used in tables.
	 * In code output flow we handle types first so we need
	 * a way to store the code for domains.
	 * The key is the type alias and the value is the
	 * type code.
	 */
	std::map<std::string, std::string> m_typeCode;
	
	/* 
	 * this holds all information that could have been
	 * dumped to output, however dependencies should
	 * be dumped first to output.
	 * For example: HandleColumn will output COMMENT on different
	 * columns. But at that moment the generated table was not
	 * dumped to m_output so we cannot write the COMMENTS 
	 * directly to m_output.
	 */
	std::ostringstream m_backpatch;

	/* output type code for the types used in a table */
	std::string OutputUsedTypes(Table& table);
	
	std::string ForeignKeyOnDelete(ConstraintForeignKey& foreignKey);

	std::string GenerateOnDeleteTriggerWhereClause(
                TColumnPtrList& sonColumns, TColumnPtrList& fatherColumns);

	std::string ConstructOnDeleteTriggers();
	
	/* assemble all output here */
	std::ostringstream m_output;	
	
	/* 
	 * for a table name remember its primary key.
	 * We use this when we try to emulate ON DELETE CASCADE
	 * for sql server. When outputing the trigger for a certain table
	 * we need to know which was the primary key for that table.
	 */
	std::map<std::string, ConstraintPrimaryKey> m_pkeys;
	
	/*
	 * for every table name remember the foreign keys that
	 * have this table as parent (that point to this table)
	 */
	std::map<std::string, std::list<ConstraintForeignKey> > m_fkeys;
		
	/* make it non-copyable */
	MssqlBackend(const MssqlBackend&);
	MssqlBackend& operator=(const MssqlBackend&);
};

#endif
