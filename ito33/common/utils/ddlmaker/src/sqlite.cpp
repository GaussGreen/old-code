/*************************************************************************** 
*  Name:        ddlmaker/src/sqlite.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: sqlite.cpp,v 1.3 2005/07/07 15:54:42 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include <algorithm> 
#include <iterator>
#include <iostream>

#include "knownattributes.h"
#include "fundamentaltypes.h"
#include "sqlite.h"
#include "variable.h"
#include "tools.h"
#include "column.h"

void SqliteBackend::Begin()
{
	m_output << "pragma synchronous = off;" << std::endl;
}

void SqliteBackend::HandleType(TypeBase& type)
{
	/* 
	 * do nothing
	 * there is no preambul for types in sqlite,
	 * one cannot create domains or user
	 * defined types
	 */
}

std::string SqliteBackend::TypeInt(TypeBase& type)
{
	return "INTEGER";
}

std::string SqliteBackend::TypeCharacter(TypeBase& type)
{
	std::ostringstream buffer;
	
	TVarBasePtr varyingVal = GET_ATTR(type, VARYING_ATTRIBUTE);
	if(varyingVal && IsTrue(*varyingVal))
		buffer << "VARCHAR";
	else
		buffer << "CHAR";
	
	TVarBasePtr lengthVal = GET_ATTR(type, LENGTH_ATTRIBUTE);
	if(lengthVal)
		buffer << "(" << lengthVal->ToSql() << ")";
	
	return buffer.str();
}

std::string SqliteBackend::TypeDecimal(TypeBase& type)
{
	std::ostringstream buffer;

	buffer << "DECIMAL";
	
	TVarBasePtr precisionVal = GET_ATTR(type, PRECISION_ATTRIBUTE);
	if(precisionVal)
		buffer << "(" << precisionVal->ToSql();

	TVarBasePtr scaleVal = GET_ATTR(type, SCALE_ATTRIBUTE);
	if(scaleVal)
		buffer << ", " << scaleVal->ToSql();
	
	buffer << ")";

	return buffer.str();
}

std::string SqliteBackend::TypeDate(TypeBase& type)
{
	TVarBasePtr containsTimeVal = GET_ATTR(type, CONTAINS_TIME_ATTRIBUTE);
	if(containsTimeVal && IsTrue(*containsTimeVal))
		return "TIMESTAMP";
	else
		return "DATE";
}

std::string SqliteBackend::TypeFloat(TypeBase& type)
{
	return "REAL";
}

std::string SqliteBackend::TypeBlob(TypeBase& type)
{
	return "BLOB";
}

void SqliteBackend::HandlePrimaryKey(ConstraintPrimaryKey& primaryKey)
{
	std::ostringstream buffer;


	/* 
	 * if this has only one column in primary list and 
	 * then generate a column constraint (we know
	 * that the higher instance of HandleTable is the one
	 * calling us, and that it has already outputed code for
	 * declaring the column in the current code map).
	 * If this has more than one column in primary list
	 * then output a table constraint.
	 */
	if(primaryKey.GetPrimaryColumns().size() == 1)
	{
		buffer << ConstraintName(primaryKey);
		
		/* 
		 * check if the type assigned to this column has
		 * the the sequence attribute set to true.
		 * If the attribute is true add the "AUTOINCREMENT"
		 * clause to the PRIMARY KEY
		 */
		buffer << " PRIMARY KEY";

		/* 
		 * ignore serial columns for the moment, sqlite 2.*
		 * has no support for 
		 * them
		 */
		/*
		TVarBasePtr sequenceTypeVal = 
			GET_ATTR(*(primaryKey.GetPrimaryColumns().front()->GetType()),
					SEQUENCE_ATTRIBUTE);
		if(sequenceTypeVal && IsTrue(*sequenceTypeVal))
			buffer << " AUTOINCREMENT";
		*/

		m_columnCode[primaryKey.GetPrimaryColumns().front()->GetName()] +=
			buffer.str();
	}
	else
	{
		buffer << " PRIMARY KEY";
		buffer << "(";
		buffer << Column::ConcatenateNames(primaryKey.GetPrimaryColumns(), ", ");
		buffer << ")";
		m_tableConstraintCode[primaryKey.GetName()] = buffer.str();
	}
}

void SqliteBackend::HandleUnique(ConstraintUnique& unique)
{
	std::ostringstream buffer;

	if(unique.GetPrimaryColumns().size() == 1)
	{
		buffer << ConstraintName(unique);
		buffer << " UNIQUE";
		m_columnCode[unique.GetPrimaryColumns().front()->GetName()] +=
			buffer.str();
	}
	else
	{
		buffer << " UNIQUE";
		buffer << "(";
		buffer << Column::ConcatenateNames(unique.GetPrimaryColumns(), ", ");
		buffer << ")";
		m_tableConstraintCode[unique.GetName()] = buffer.str();
	}
}

void SqliteBackend::HandleForeignKey(ConstraintForeignKey& foreignKey)
{
	std::ostringstream buffer;

	if(foreignKey.GetPrimaryColumns().size() == 1)
	{
		buffer << ConstraintName(foreignKey);
		buffer << " REFERENCES ";
		buffer << foreignKey.GetFatherTable()->GetName();
		
		/* there can be only one father column for a single son column */
		buffer << "(" << foreignKey.GetSecondaryColumns().front()->GetName() << ")";

		buffer << ForeignKeyOnDelete(foreignKey);
		
		m_columnCode[foreignKey.GetPrimaryColumns().front()->GetName()] +=
			buffer.str();
	}
	else
	{
		buffer << " FOREIGN KEY ";
		
		buffer << "(";
		buffer << Column::ConcatenateNames(foreignKey.GetPrimaryColumns(), ", ");
		buffer << ")";
		
		buffer << " REFERENCES " << foreignKey.GetFatherTable()->GetName();
		
		buffer << "(" ;
		buffer << Column::ConcatenateNames(
				foreignKey.GetSecondaryColumns(), ", ");
		buffer << ")";

		m_tableConstraintCode[foreignKey.GetName()] = buffer.str();
	}
}

void SqliteBackend::HandleNotNull(ConstraintNull& notNull)
{
	std::ostringstream buffer;
	
	if(notNull.NotNull())
	{
		buffer << ConstraintName(notNull);
		buffer << " NOT NULL";
		m_columnCode[notNull.GetPrimaryColumns().front()->GetName()] +=
			buffer.str();
	}
}

std::string SqliteBackend::TypeSelector(TypeBase& type)
{
	std::string toReturn;

	switch(type.GetTypeId())
	{
		case BASE_TYPE_INTEGER:
			toReturn = TypeInt(type);
			break;
		case BASE_TYPE_CHARACTER:
			toReturn = TypeCharacter(type);
			break;
		case BASE_TYPE_DECIMAL:
			toReturn = TypeDecimal(type);
			break;
		case BASE_TYPE_DATE:
			toReturn = TypeDate(type);
			break;
		case BASE_TYPE_FLOAT:
			toReturn = TypeFloat(type);
			break;
		case BASE_TYPE_BLOB:
			toReturn = TypeBlob(type);
			break;
		default:
			Error("Unrecognized type");
	}

	return toReturn;
}

void SqliteBackend::HandleColumn(Column& column)
{
	std::string columnName = column.GetName();
	std::ostringstream buffer;
	
	/* we need to handle the constraints */
	buffer << columnName;
	buffer << " ";
	
	buffer << TypeSelector(*(column.GetType()));

	TVarBasePtr defaultVal = GET_ATTR(column, DEFAULT_ATTRIBUTE);
	if(defaultVal)
		buffer << " DEFAULT " << defaultVal->ToSql();
	
	m_columnCode[column.GetName()] = buffer.str();
	
	/* 
	 * constraint handling will add its own code to columnCode[column.GetName()] 
	 * (of course, only if the constraint is a column constraint
	 * and not a table constraint)
	 */
}

void SqliteBackend::HandleTable(Table& table)
{
	std::ostringstream buffer;

	buffer << "CREATE TABLE " << table.GetName() << std::endl;
	buffer << "(" << std::endl;

	TColumnPtrList::const_iterator tempIterator;
	for(tempIterator = table.GetColumns().begin();
			tempIterator != table.GetColumns().end();
			++tempIterator)
	{
		(*tempIterator)->CodeGenerator(*this);
	}

	TConstraintBasePtrList::const_iterator tempIterator2;
	for(tempIterator2 = table.GetConstraints().begin();
			tempIterator2 != table.GetConstraints().end();
			++tempIterator2)
	{
		(*tempIterator2)->CodeGenerator(*this);
	}

	std::map<std::string, std::string>::iterator tempIterator3;
	for(tempIterator3 = m_columnCode.begin();
			tempIterator3 != m_columnCode.end();
			++tempIterator3)
	{
		buffer << (*tempIterator3).second;
		if(std::distance(tempIterator3, m_columnCode.end()) > 1)
		{
			buffer << ", ";
			buffer << std::endl;
		}
	}
	
	if(m_tableConstraintCode.size())
		buffer << ", " << std::endl;
	
	std::map<std::string, std::string>::iterator tempIterator4;
	for(tempIterator4 = m_tableConstraintCode.begin();
			tempIterator4 != m_tableConstraintCode.end();
			++tempIterator4)
	{
		buffer << (*tempIterator4).second;
		if(std::distance(tempIterator4, m_tableConstraintCode.end()) > 1)
		{
			buffer << ", ";
			buffer << std::endl;
		}
	}

	buffer << std::endl << ");" << std::endl;
	
	/* dump everything output */
	m_output << buffer.str();

	m_columnCode.clear();
	m_tableConstraintCode.clear();
}

std::string SqliteBackend::ConstraintName(ConstraintBase& constraint)
{
	if(constraint.GetName().size() > 0)
	{
		return " CONSTRAINT " + constraint.GetName();
	}
	else
		return "";
}

std::string SqliteBackend::GetOutput()
{
	return m_output.str();
}

std::string SqliteBackend::ForeignKeyOnDelete(
		ConstraintForeignKey& foreignKey)
{
	std::string toReturn;
	
	TVarBasePtr onDeleteVal = GET_ATTR(foreignKey, ON_DELETE_ATTRIBUTE);
	if(onDeleteVal && onDeleteVal->ToString() == "cascade")
		toReturn = " ON DELETE CASCADE";
	
	return toReturn;
}
