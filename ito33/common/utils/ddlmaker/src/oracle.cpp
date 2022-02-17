/*************************************************************************** 
*  Name:        ddlmaker/src/oracle.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: oracle.cpp,v 1.1 2005/09/05 09:04:12 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include <algorithm> 
#include <iterator>
#include <iostream>

#include "knownattributes.h"
#include "fundamentaltypes.h"
#include "oracle.h"
#include "variable.h"
#include "tools.h"
#include "column.h"

/*
 * There is a major difference between Oracle's backend 
 * and the backends for the other sgbds. Oracle has no
 * support for user defined (domain) types.
 */
void OracleBackend::HandleType(TypeBase& type)
{
}

std::string OracleBackend::DomainDeclGen(TypeBase& type, 
		const std::string& tableName,
		const std::string& columnName)
{
	std::ostringstream buffer;

	buffer << TypeSelector(type);

	TVarBasePtr defaultVal = GET_ATTR(type, DEFAULT_ATTRIBUTE);
	if(defaultVal)
		buffer << " DEFAULT " << defaultVal->ToSql();

	TVarBasePtr notNullVal = GET_ATTR(type, NOT_NULL_ATTRIBUTE);
	if(notNullVal && IsTrue(*notNullVal))
		buffer << " NOT NULL";
	
	TVarBasePtr domainVal = GET_ATTR(type, DOMAIN_ATTRIBUTE);
	if(domainVal)
	{
		buffer << " CONSTRAINT ";
		std::string ckName("ck_");
		ckName.append(tableName + "_" + columnName);
		AdjustLength(ckName);
		buffer << ckName << " CHECK (" << columnName << 
			" IN (" << domainVal->ToSql() << "))";
	}

	return buffer.str();
}

std::string OracleBackend::TypeInt(TypeBase& type)
{
	return "INT";
}

std::string OracleBackend::TypeCharacter(TypeBase& type)
{
	std::ostringstream buffer;
	
	TVarBasePtr varyingVal = GET_ATTR(type, VARYING_ATTRIBUTE);
	if(varyingVal && IsTrue(*varyingVal))
		buffer << "VARCHAR2";
	else
		buffer << "CHAR";
	
	TVarBasePtr lengthVal = GET_ATTR(type, LENGTH_ATTRIBUTE);
	if(lengthVal)
		buffer << "(" << lengthVal->ToSql() << ")";
	
	return buffer.str();
}

std::string OracleBackend::TypeDecimal(TypeBase& type)
{
	std::ostringstream buffer;

	buffer << "NUMERIC";
	
	TVarBasePtr precisionVal = GET_ATTR(type, PRECISION_ATTRIBUTE);
	if(precisionVal)
		buffer << "(" << precisionVal->ToSql();

	TVarBasePtr scaleVal = GET_ATTR(type, SCALE_ATTRIBUTE);
	if(scaleVal)
		buffer << ", " << scaleVal->ToSql();
	
	buffer << ")";

	return buffer.str();
}

std::string OracleBackend::TypeDate(TypeBase& type)
{
	TVarBasePtr containsTimeVal = GET_ATTR(type, CONTAINS_TIME_ATTRIBUTE);
	if(containsTimeVal && IsTrue(*containsTimeVal))
		return "TIMESTAMP";
	else
		return "DATE";
}

std::string OracleBackend::TypeFloat(TypeBase& type)
{
	TVarBasePtr precisionVal = GET_ATTR(type, PRECISION_ATTRIBUTE);

	if(precisionVal)
		return "FLOAT("+precisionVal->ToSql()+")";
	else
		return "FLOAT";
}

std::string OracleBackend::TypeBlob(TypeBase& type)
{
	std::string toReturn;
	TVarBasePtr containsTextVal = GET_ATTR(type, CONTAINS_TEXT_ATTRIBUTE);

	if(containsTextVal && IsTrue(*containsTextVal))
		toReturn = "TEXT";
	else
	{
		Error("The Postgresql backend doesn't have support for blobs (except for TEXT)");
	}

	return toReturn;
}

void OracleBackend::HandlePrimaryKey(ConstraintPrimaryKey& primaryKey)
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
	buffer << ConstraintName(primaryKey);
	buffer << " PRIMARY KEY";
	if(primaryKey.GetPrimaryColumns().size() == 1)
	{
		m_columnCode[primaryKey.GetPrimaryColumns().front()->GetName()] +=
			buffer.str();
	}
	else
	{
		buffer << "(";
		buffer << Column::ConcatenateNames(primaryKey.GetPrimaryColumns(), ", ");
		buffer << ")";
		m_tableConstraintCode[primaryKey.GetName()] = buffer.str();
	}
}

void OracleBackend::HandleUnique(ConstraintUnique& unique)
{
	std::ostringstream buffer;

	buffer << ConstraintName(unique);
	buffer << " UNIQUE";
	if(unique.GetPrimaryColumns().size() == 1)
	{
		m_columnCode[unique.GetPrimaryColumns().front()->GetName()] +=
			buffer.str();
	}
	else
	{
		buffer << "(";
		buffer << Column::ConcatenateNames(unique.GetPrimaryColumns(), ", ");
		buffer << ")";
		m_tableConstraintCode[unique.GetName()] = buffer.str();
	}
}

void OracleBackend::HandleForeignKey(ConstraintForeignKey& foreignKey)
{
	std::ostringstream buffer;

	buffer << ConstraintName(foreignKey);
	if(foreignKey.GetPrimaryColumns().size() == 1)
	{
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

void OracleBackend::HandleNotNull(ConstraintNull& notNull)
{
	std::ostringstream buffer;
	std::string columnName = notNull.GetPrimaryColumns().front()->GetName();

	if(notNull.NotNull())
	{
		buffer << ConstraintName(notNull);
		buffer << " NOT NULL";
		
		// do not add a NOT NULL constraint if it already exists
		std::map<std::string, std::string>::const_iterator pos;
		pos = m_columnCode.find(columnName);
		if ( pos == m_columnCode.end() ||
			pos->second.find("NOT NULL") == std::string::npos )
		{
			m_columnCode[columnName] +=
				buffer.str();
		}
	}
}

std::string OracleBackend::TypeSelector(TypeBase& type)
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

void OracleBackend::HandleColumn(Column& column)
{
	std::string columnName = column.GetName();
	std::ostringstream buffer;
	
	/* we need to handle the constraints */
	buffer << columnName;
	buffer << " ";
	
	/* 
	 * either output the alias for a type which had its domain
	 * instantiated or simply try to inline the type
	 */
	TVarBasePtr domainInstanceVal = 
		GET_ATTR(*(column.GetType()), DOMAIN_INSTANCE_ATTRIBUTE);
	if(domainInstanceVal && IsTrue(*domainInstanceVal))
		// if this is a domain output check constraint
		buffer << DomainDeclGen(*(column.GetType()), 
				column.GetOwnerTable()->GetName(),
				column.GetName());
	else
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
	
	TVarBasePtr commentVal = GET_ATTR(column, COMMENT_ATTRIBUTE);
	if(commentVal)
	{
		m_backpatch << "COMMENT ON COLUMN ";
		m_backpatch << column.GetOwnerTable()->GetName() << "." 
			<< columnName;
		m_backpatch << " IS " << commentVal->ToSql() << ";" << std::endl;
	}
}

void OracleBackend::HandleTable(Table& table)
{
	std::ostringstream buffer;

	buffer << OutputSequenceCreation(table);
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
	m_output << m_backpatch.str();

	m_backpatch.str("");
	m_columnCode.clear();
	m_tableConstraintCode.clear();
}

std::string OracleBackend::OutputSequenceCreation(Table& table)
{
	std::ostringstream buffer;
	TColumnPtrList::const_iterator iterator;

	for(iterator = table.GetColumns().begin(); 
			iterator != table.GetColumns().end(); ++iterator)
	{
		std::map<std::string, std::string>::iterator typeCodeIterator =
			m_typeCode.find((*iterator)->GetType()->GetAlias());
		
		TVarBasePtr sequenceVal = GET_ATTR(*((*iterator)->GetType()), 
				SEQUENCE_ATTRIBUTE);
		if ( sequenceVal && IsTrue(*sequenceVal) )
		{
			buffer << "CREATE SEQUENCE " << table.GetName() <<
				"_" << (*iterator)->GetName() << "_seq;" <<
				std::endl;
		}
	}

	return buffer.str();
}

std::string OracleBackend::ConstraintName(ConstraintBase& constraint)
{
	if(constraint.GetName().size() > 0)
	{
		std::string constraintName(constraint.GetName());
		AdjustLength(constraintName);
		return " CONSTRAINT " + constraintName;
	}
	else
		return "";
}

std::string OracleBackend::GetOutput()
{
	return m_output.str();
}

std::string OracleBackend::ForeignKeyOnDelete(
		ConstraintForeignKey& foreignKey)
{
	std::string toReturn;
	
	TVarBasePtr onDeleteVal = GET_ATTR(foreignKey, ON_DELETE_ATTRIBUTE);
	if(onDeleteVal && onDeleteVal->ToString() == "cascade")
		toReturn = " ON DELETE CASCADE";
	
	return toReturn;
}

void OracleBackend::AdjustLength(std::string& str)
{
	str.assign(str, 0, 30);
}

