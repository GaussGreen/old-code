/*************************************************************************** 
*  Name:        ddlmaker/include/oracle.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: oracle.h,v 1.1 2005/09/05 09:04:34 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __ORACLE_H__
#define __ORACLE_H__

#include "backend.h"

#include "table.h"
#include "column.h"
#include "constraint.h"

#include <string>
#include <map>
#include <sstream>

class OracleBackend
:public Backend
{
public:
	OracleBackend()
	{
	}
	
	/* does nothing */
	void Begin()
	{
	}
	
	/* called for each table to generate sql statements */
	void HandleTable(Table& table);

	/* 
	 * called for all types in the symbol table, it 
	 * prepares the types for usage.
	 */
	void HandleType(TypeBase& type);

	/*
	 * generates code for a column
	 */
	void HandleColumn(Column& column);

	/*
	 * we need to use the visitor pattern to 
	 * generate code for the constraints as we
	 * use them in a polymorphical manner
	 */
	void HandlePrimaryKey(ConstraintPrimaryKey& primaryKey);
	void HandleForeignKey(ConstraintForeignKey& foreignKey);
	void HandleUnique(ConstraintUnique& unique);
	void HandleNotNull(ConstraintNull& notNull);

	std::string GetOutput();

	/* does nothing */
	void End()
	{
	}

private:
	/* 
	 * domain constraint creation for types.
	 * we know that oracle has no support whatsoever for user
	 * defined types (domains).
	 * we circumvent this by translating the domain definition
	 * into a column "check" constraint.
	 */
	std::string DomainDeclGen(TypeBase& type, 
			const std::string& tableName,
			const std::string& columnName);

	/* call the appropriate method for a certain type */
	std::string TypeSelector(TypeBase& type);
	
	/* 
	 * these functions output base types recognized 
	 * by the pgsql sgbd engine.
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
	 * Use this map to store check constraint info for a
	 * certain type.
	 */
	std::map<std::string, std::string> m_typeCode;
	
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
	 * this holds all information that could have been
	 * dumped to output, however dependencies should
	 * be dumped first to output.
	 * For example: HandleColumn will output COMMENT on different
	 * columns. But at that moment the generated table was not
	 * dumped to m_output so we cannot write the COMMENTS 
	 * directly to m_output.
	 */
	std::ostringstream m_backpatch;

	/* 
	 * when dealing with tables that have a PK we can
	 * generate the unique id by using a sequence.
	 * but that sequence needs to be created before we
	 * retrieve the id from it.
	 */
	std::string OutputSequenceCreation(Table& table);
	
	std::string ForeignKeyOnDelete(ConstraintForeignKey& foreignKey);
	
	/* assemble all output here */
	std::ostringstream m_output;	

	/* all identifiers in a oracle database must have less than 30 chars */
	void AdjustLength(std::string& str);

	/* make it non-copyable */
	OracleBackend(const OracleBackend&);
	OracleBackend& operator=(const OracleBackend&);
};
#endif
