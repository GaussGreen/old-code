/*************************************************************************** 
*  Name:        ddlmaker/src/pgsql.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: pgsql.cpp,v 1.2 2005/07/06 09:03:02 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include <algorithm> 
#include <iterator>
#include <iostream>

#include "knownattributes.h"
#include "fundamentaltypes.h"
#include "pgsql.h"
#include "variable.h"
#include "tools.h"
#include "column.h"

void PostgresqlBackend::HandleType(TypeBase& type)
{
	/*
	 * Postgres has a funny way (the simplest) of handling domains.
	 * The types are not hierarchical.
	 * When creating a domain we cannot create it using
	 * an existing domain as starting point. We NEED to create it
	 * using a base recognizable type.
	 * Example: you cannot do:
	 * CREATE DOMAIN CHAR_1 AS CHAR(1);
	 * CREATE DOMAIN SOME_DOMAIN AS CHAR_1;
	 *
	 * But instead you need to CREATE DOMAIN SOME_DOMAIN AS CHAR(1).
	 * This is why we use a form of backpatching for type code. We put
	 * all the code in a map and before each table we see which
	 * types we use in it. We output exactly those types.
	 * A direct effect is that we only output the types used.
	 */

	/* 
	 * for certain type (for example SERIAL) you cannot create
	 * a domain. This imposes a new attribute in the frontend
	 * for types "domainInstance". If "domainInstance" is true
	 * then a new domain will be instantiated for the aliased
	 * type. If not we try to inline the base type at the 
	 * column declaration point.
	 */
	TVarBasePtr domainInstanceVal = GET_ATTR(type, DOMAIN_INSTANCE_ATTRIBUTE);
	TVarBasePtr domainVal = GET_ATTR(type, DOMAIN_ATTRIBUTE);
	if(!type.IsBase())
	{
		if((domainInstanceVal && IsTrue(*domainInstanceVal)) ||
			domainVal)
		{
			m_typeCode.insert(std::make_pair(type.GetAlias(), 
						DomainDeclGen(type)));
		}	
	}
}

std::string PostgresqlBackend::DomainDeclGen(TypeBase& type)
{
	std::ostringstream buffer;

	buffer << "CREATE DOMAIN " << type.GetAlias();
	buffer << " AS ";
	buffer << TypeSelector(type);

	TVarBasePtr defaultVal = GET_ATTR(type, DEFAULT_ATTRIBUTE);
	if(defaultVal)
		buffer << " DEFAULT " << defaultVal->ToSql();

	TVarBasePtr notNullVal = GET_ATTR(type, NOT_NULL_ATTRIBUTE);
	if(notNullVal && IsTrue(*notNullVal))
		buffer << " NOT NULL";
	else
		buffer << " NULL";
	
	TVarBasePtr domainVal = GET_ATTR(type, DOMAIN_ATTRIBUTE);
	if(domainVal)
	{
		buffer << " CONSTRAINT ck" << type.GetAlias();
		buffer << " CHECK (VALUE IN (" << domainVal->ToSql() << "))";
	}
	buffer << ";";

	return buffer.str();
}

std::string PostgresqlBackend::TypeInt(TypeBase& type)
{
	TVarBasePtr sequenceVal = GET_ATTR(type, SEQUENCE_ATTRIBUTE);
	if(sequenceVal && IsTrue(*sequenceVal))
		return "SERIAL";
	else
		return "INT";
}

std::string PostgresqlBackend::TypeCharacter(TypeBase& type)
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

std::string PostgresqlBackend::TypeDecimal(TypeBase& type)
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

std::string PostgresqlBackend::TypeDate(TypeBase& type)
{
	TVarBasePtr containsTimeVal = GET_ATTR(type, CONTAINS_TIME_ATTRIBUTE);
	if(containsTimeVal && IsTrue(*containsTimeVal))
		return "TIMESTAMP";
	else
		return "DATE";
}

std::string PostgresqlBackend::TypeFloat(TypeBase& type)
{
	TVarBasePtr precisionVal = GET_ATTR(type, PRECISION_ATTRIBUTE);

	if(precisionVal)
		return "FLOAT("+precisionVal->ToSql()+")";
	else
		return "FLOAT";
}

std::string PostgresqlBackend::TypeBlob(TypeBase& type)
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

void PostgresqlBackend::HandlePrimaryKey(ConstraintPrimaryKey& primaryKey)
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

	TVarBasePtr commentVal = GET_ATTR(primaryKey, COMMENT_ATTRIBUTE);
	if(commentVal)
	{
		m_backpatch << "COMMENT ON CONSTRAINT ";
		m_backpatch << primaryKey.GetName() << " ON " 
			<< primaryKey.GetOwnerTable()->GetName();
		m_backpatch << " IS " << commentVal->ToSql() << ";" << std::endl;
	}
}

void PostgresqlBackend::HandleUnique(ConstraintUnique& unique)
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
	
	TVarBasePtr commentVal = GET_ATTR(unique, COMMENT_ATTRIBUTE);
	if(commentVal)
	{
		m_backpatch << "COMMENT ON CONSTRAINT ";
		m_backpatch << unique.GetName() << " ON " 
			<< unique.GetOwnerTable()->GetName();
		m_backpatch << " IS " << commentVal->ToSql() << ";" << std::endl;
	}
}

void PostgresqlBackend::HandleForeignKey(ConstraintForeignKey& foreignKey)
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

	TVarBasePtr commentVal = GET_ATTR(foreignKey, COMMENT_ATTRIBUTE);
	if(commentVal)
	{
		m_backpatch << "COMMENT ON CONSTRAINT ";
		m_backpatch << foreignKey.GetName() << " ON " 
			<< foreignKey.GetOwnerTable()->GetName();
		m_backpatch << " IS " << commentVal->ToSql() << ";" << std::endl;
	}
}

void PostgresqlBackend::HandleNotNull(ConstraintNull& notNull)
{
	std::ostringstream buffer;
	
	if(notNull.NotNull())
	{
		buffer << ConstraintName(notNull);
		buffer << " NOT NULL";
		m_columnCode[notNull.GetPrimaryColumns().front()->GetName()] +=
			buffer.str();
	}
	
	TVarBasePtr commentVal = GET_ATTR(notNull, COMMENT_ATTRIBUTE);
	if(commentVal)
	{
		m_backpatch << "COMMENT ON CONSTRAINT ";
		m_backpatch << notNull.GetName() << " ON " 
			<< notNull.GetOwnerTable()->GetName();
		m_backpatch << " IS " << commentVal->ToSql() << ";" << std::endl;
	}
}

std::string PostgresqlBackend::TypeSelector(TypeBase& type)
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

void PostgresqlBackend::HandleColumn(Column& column)
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
		buffer << column.GetType()->GetAlias();
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

void PostgresqlBackend::HandleTable(Table& table)
{
	std::ostringstream buffer;

	buffer << OutputUsedTypes(table);
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

std::string PostgresqlBackend::OutputUsedTypes(Table& table)
{
	std::ostringstream buffer;
	TColumnPtrList::const_iterator iterator;

	for(iterator = table.GetColumns().begin(); 
			iterator != table.GetColumns().end(); ++iterator)
	{
		std::map<std::string, std::string>::iterator typeCodeIterator =
			m_typeCode.find((*iterator)->GetType()->GetAlias());
		if(typeCodeIterator != m_typeCode.end())
		{
			buffer << (*typeCodeIterator).second << std::endl;
			m_typeCode.erase(typeCodeIterator);
		}
	}

	return buffer.str();
}

std::string PostgresqlBackend::ConstraintName(ConstraintBase& constraint)
{
	if(constraint.GetName().size() > 0)
	{
		return " CONSTRAINT " + constraint.GetName();
	}
	else
		return "";
}

std::string PostgresqlBackend::GetOutput()
{
	return m_output.str();
}

std::string PostgresqlBackend::ForeignKeyOnDelete(
		ConstraintForeignKey& foreignKey)
{
	std::string toReturn;
	
	TVarBasePtr onDeleteVal = GET_ATTR(foreignKey, ON_DELETE_ATTRIBUTE);
	if(onDeleteVal && onDeleteVal->ToString() == "cascade")
		toReturn = " ON DELETE CASCADE";
	
	return toReturn;
}
