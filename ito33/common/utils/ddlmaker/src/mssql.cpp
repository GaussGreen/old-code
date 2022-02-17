/*************************************************************************** 
*  Name:        ddlmaker/src/mssql.cpp                              
*  Purpose:	Sql Server Backend                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: mssql.cpp,v 1.5 2005/08/29 16:45:56 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include <algorithm> 
#include <iterator>
#include <iostream>

#include "knownattributes.h"
#include "fundamentaltypes.h"
#include "mssql.h"
#include "variable.h"
#include "tools.h"
#include "column.h"
#include "constraint.h"

void MssqlBackend::HandleType(TypeBase& type)
{
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

std::string MssqlBackend::DomainDeclGen(TypeBase& type)
{
	std::ostringstream buffer;

	buffer << "sp_addtype " << type.GetAlias();
	buffer << ", " << "'" << TypeSelector(type) <<"'" << std::endl;
	buffer << "GO" << std::endl;

	TVarBasePtr defaultVal = GET_ATTR(type, DEFAULT_ATTRIBUTE);
	if(defaultVal)
	{
		buffer << "CREATE DEFAULT DF_" << type.GetAlias();
		buffer << " AS " << defaultVal->ToSql() << std::endl;
		buffer << "GO" << std::endl;

		buffer << "sp_bindefault " << "'DF_" << type.GetAlias() << "', ";
		buffer << "'" << type.GetAlias() << "'" << std::endl;
		buffer << "GO" << std::endl;
	}

	TVarBasePtr domainVal = GET_ATTR(type, DOMAIN_ATTRIBUTE);
	if(domainVal)
	{
		buffer << "CREATE RULE RL_" << type.GetAlias();
		buffer << " AS (@value in (" << domainVal->ToSql() << "))" << std::endl;
		buffer << "GO" << std::endl;
	}

	return buffer.str();
}

std::string MssqlBackend::TypeInt(TypeBase& type)
{
	TVarBasePtr sequenceVal = GET_ATTR(type, SEQUENCE_ATTRIBUTE);
	if(sequenceVal && IsTrue(*sequenceVal))
		return "INTEGER IDENTITY(1, 1)";
	else
		return "INTEGER";
}

std::string MssqlBackend::TypeCharacter(TypeBase& type)
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

std::string MssqlBackend::TypeDecimal(TypeBase& type)
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

std::string MssqlBackend::TypeDate(TypeBase& type)
{
	return "DATETIME";
}

std::string MssqlBackend::TypeFloat(TypeBase& type)
{
	TVarBasePtr precisionVal = GET_ATTR(type, PRECISION_ATTRIBUTE);

	if(precisionVal)
		return "FLOAT("+precisionVal->ToSql()+")";
	else
		return "FLOAT";
}

std::string MssqlBackend::TypeBlob(TypeBase& type)
{
	std::string toReturn;
	TVarBasePtr containsTextVal = GET_ATTR(type, CONTAINS_TEXT_ATTRIBUTE);

	if(containsTextVal && IsTrue(*containsTextVal))
		toReturn = "TEXT";
	else
	{
		toReturn = "IMAGE";
	}

	return toReturn;
}

void MssqlBackend::HandlePrimaryKey(ConstraintPrimaryKey& primaryKey)
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

	/* mark the constraint in the primary keys table */
	std::string primaryKeyTableName = primaryKey.GetOwnerTable()->GetName();
	if(m_pkeys.find(primaryKeyTableName) == m_pkeys.end())
		m_pkeys[primaryKeyTableName] = primaryKey; 
	
	TVarBasePtr commentVal = GET_ATTR(primaryKey, COMMENT_ATTRIBUTE);
	if(commentVal)
	{
		m_backpatch << "EXECUTE sp_addextendedproperty N'MS_Description', ";
		m_backpatch << "N" << commentVal->ToSql() << ", ";
		m_backpatch << "N'user', N'dbo', N'table', N'" << primaryKey.GetOwnerTable()->GetName() << "', ";
		m_backpatch << "N'constraint', ";
		m_backpatch << "N'" << primaryKey.GetName() << "'" << std::endl;
		m_backpatch << "GO" << std::endl;
	}
}

void MssqlBackend::HandleUnique(ConstraintUnique& unique)
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
		m_backpatch << "EXECUTE sp_addextendedproperty N'MS_Description', ";
		m_backpatch << "N" << commentVal->ToSql() << ", ";
		m_backpatch << "N'user', N'dbo', N'table', N'" << unique.GetOwnerTable()->GetName() 
			<< "', " << "N'constraint', ";
		m_backpatch << "N'" << unique.GetName() << "'" << std::endl;
		m_backpatch << "GO" << std::endl;
	}
}

void MssqlBackend::HandleForeignKey(ConstraintForeignKey& foreignKey)
{
	std::ostringstream buffer;

	/* 
	 * Stupid mssql barfs when more than one foreign key
	 * is in the column constraint line.
	 * We put all foreign keys as table constraints
	 */
	buffer << ConstraintName(foreignKey);
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

	/* we populate the m_fkeys */
	std::string targetTableName = 
		foreignKey.GetSecondaryColumns().front()->GetOwnerTable()->GetName();
	m_fkeys[targetTableName].push_back(foreignKey);
	
	TVarBasePtr commentVal = GET_ATTR(foreignKey, COMMENT_ATTRIBUTE);
	if(commentVal)
	{
		m_backpatch << "EXECUTE sp_addextendedproperty N'MS_Description', ";
		m_backpatch << "N" << commentVal->ToSql() << ", ";
		m_backpatch << "N'user', N'dbo', N'table', N'" << foreignKey.GetOwnerTable()->GetName()
			<< "', " << "N'constraint', ";
		m_backpatch << "N'" << foreignKey.GetName() << "'" << std::endl;
		m_backpatch << "GO" << std::endl;
	}
}

void MssqlBackend::HandleNotNull(ConstraintNull& notNull)
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
		m_backpatch << "EXECUTE sp_addextendedproperty N'MS_Description', ";
		m_backpatch << "N" << commentVal->ToSql() << ", ";
		m_backpatch << "N'user', N'dbo', N'table', N'" << notNull.GetOwnerTable()->GetName() 
			<< "', " << "N'constraint', ";
		m_backpatch << "N'" << notNull.GetName() << "'" << std::endl;
		m_backpatch << "GO" << std::endl;
	}
}

std::string MssqlBackend::TypeSelector(TypeBase& type)
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

void MssqlBackend::HandleColumn(Column& column)
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

	TVarBasePtr collateVal = GET_ATTR(*(column.GetType()), "collate");
	if(collateVal)
	{
		/* 
		 * we call ToString() instead of ToSql
		 * because we don't want the string to be 
		 * escaped inside single quotes
		 */
		buffer << " COLLATE " << collateVal->ToString();
	}

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
		m_backpatch << "EXECUTE sp_addextendedproperty N'MS_Description', ";
		m_backpatch << "N" << commentVal->ToSql() << ", ";
		m_backpatch << "N'user', N'dbo', N'table', ";
		m_backpatch << "N'" << column.GetOwnerTable()->GetName() << "', ";
		m_backpatch << "N'column', " << "N'" << column.GetName() << "'" << std::endl;
		m_backpatch << "GO" << std::endl;
	}
}

void MssqlBackend::HandleTable(Table& table)
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
	
	TVarBasePtr commentVal = GET_ATTR(table, COMMENT_ATTRIBUTE);
	if(commentVal)
	{
		m_backpatch << "EXECUTE sp_addextendedproperty N'MS_Description', ";
		m_backpatch << "N" << commentVal->ToSql() << ", ";
		m_backpatch << "N'user', N'dbo', N'table', ";
		m_backpatch << "N'" << table.GetName() << "', ";
		m_backpatch << "NULL, NULL;" << std::endl;
		m_backpatch << "GO" << std::endl;
	}

	buffer << std::endl << ")" << std::endl;
	buffer << "GO" << std::endl;
	
	/* dump everything output */
	m_output << buffer.str();
	m_output << m_backpatch.str();

	m_backpatch.str("");
	m_columnCode.clear();
	m_tableConstraintCode.clear();
}

std::string MssqlBackend::OutputUsedTypes(Table& table)
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

std::string MssqlBackend::ConstraintName(ConstraintBase& constraint)
{
	if(constraint.GetName().size() > 0)
	{
		return " CONSTRAINT " + constraint.GetName();
	}
	else
		return "";
}

std::string MssqlBackend::GetOutput()
{
	return m_output.str();
}

std::string MssqlBackend::ForeignKeyOnDelete(
		ConstraintForeignKey& foreignKey)
{
	std::string toReturn;
	
	TVarBasePtr onDeleteVal = GET_ATTR(foreignKey, ON_DELETE_ATTRIBUTE);
	if(onDeleteVal && onDeleteVal->ToString() == "cascade")
		toReturn = " ON DELETE CASCADE";
	
	return toReturn;
}

std::string MssqlBackend::GenerateOnDeleteTriggerWhereClause(
		TColumnPtrList& sonColumns, TColumnPtrList& fatherColumns)
{
	std::ostringstream buffer;

	TColumnPtrList::iterator sonColumnsIterator;
	TColumnPtrList::iterator fatherColumnsIterator;
	
	/* start generating the WHERE clause */
	for(sonColumnsIterator = sonColumns.begin(), 
		fatherColumnsIterator = fatherColumns.begin();
		sonColumnsIterator != sonColumns.end() &&
		fatherColumnsIterator != fatherColumns.end();
		++sonColumnsIterator, ++fatherColumnsIterator)
	{
		buffer << (*sonColumnsIterator)->GetName() << " IN ";
		buffer << "(SELECT " << (*fatherColumnsIterator)->GetName()
		<< " FROM DELETED) ";
					
		/* we can measure the distance for son or father iterator
		 * it should be the same thing
		 */
		if(std::distance(sonColumnsIterator, sonColumns.end()) > 1)
			buffer << " AND ";
	}

	return buffer.str();
}

std::string MssqlBackend::ConstructOnDeleteTriggers()
{
	std::map<std::string, std::list<ConstraintForeignKey> >::iterator iterator;
	std::list<ConstraintForeignKey>::iterator fkIterator;
	bool triggerCreated = false;
	std::ostringstream buffer;
	
	for(iterator = m_fkeys.begin(); iterator != m_fkeys.end(); 
			++iterator)
	{
		for(fkIterator = iterator->second.begin();
				fkIterator != iterator->second.end();
				++fkIterator)
		{
		        TVarBasePtr onDeleteVal = GET_ATTR(*fkIterator, ON_DELETE_ATTRIBUTE);
		        if(onDeleteVal && onDeleteVal->ToString() == "cascade")	
			{
				if(!triggerCreated)
				{
					buffer << "CREATE TRIGGER trigger_on_delete_cascade_"
						<< iterator->first << " ON " << iterator->first
						<< std::endl;
					buffer << "  INSTEAD OF DELETE" << std::endl;
					buffer << "AS" << std::endl;
					buffer << "BEGIN" << std::endl;
					
					TColumnPtrList pkeysList = m_pkeys[iterator->first].GetPrimaryColumns();
					buffer << "DELETE FROM " << iterator->first << " WHERE ";
					buffer << GenerateOnDeleteTriggerWhereClause(pkeysList, pkeysList);
					buffer << ";" << std::endl;
					
					triggerCreated = true;
				}
				TColumnPtrList sonColumns = fkIterator->GetPrimaryColumns();
				std::string sonTableName = sonColumns.front()->GetOwnerTable()->GetName();
				TColumnPtrList fatherColumns = 
					m_pkeys[iterator->first].GetPrimaryColumns();

				buffer << "DELETE FROM " << sonTableName 
					<< " WHERE ";
				buffer << GenerateOnDeleteTriggerWhereClause(sonColumns, fatherColumns);
				buffer << ";" << std::endl;
			}
		}

		if(triggerCreated)
		{
			buffer << "END" << std::endl;
			buffer << "GO" << std::endl;

			triggerCreated = false;
		}
	}

	return buffer.str();
}

void MssqlBackend::End()
{
	/* here we will dump all the code that goes at the
	 * end of the script (namely trigger creation
	 * to correctly handle ON DELETE CASCADE
	 * for mssql
	 */
	m_output << ConstructOnDeleteTriggers();
}
