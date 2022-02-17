/*************************************************************************** 
*  Name:        ddlmaker/src/type.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: type.cpp,v 1.2 2005/07/06 09:03:02 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include "type.h"

TTypeBasePtrList TypeBase::BuildTypeBaseTable()
{
	TTypeBasePtrList toReturn;

	toReturn.push_back(boost::shared_ptr<TypeBase>(
				new TypeBase(BASE_TYPE_INTEGER)));
	toReturn.push_back(boost::shared_ptr<TypeBase>(
				new TypeBase(BASE_TYPE_CHARACTER)));
	toReturn.push_back(boost::shared_ptr<TypeBase>(
				new TypeBase(BASE_TYPE_DECIMAL)));
	toReturn.push_back(boost::shared_ptr<TypeBase>(
				new TypeBase(BASE_TYPE_DATE)));
	toReturn.push_back(boost::shared_ptr<TypeBase>(
				new TypeBase(BASE_TYPE_FLOAT)));
	toReturn.push_back(boost::shared_ptr<TypeBase>(
				new TypeBase(BASE_TYPE_BLOB)));

	return toReturn;
}

std::string TypeBase::GetNameById(int id)
{
	switch(id)
	{
		case BASE_TYPE_INTEGER:
			return std::string("int");
		case BASE_TYPE_CHARACTER:
			return std::string("char");
		case BASE_TYPE_DECIMAL:
			return std::string("decimal");
		case BASE_TYPE_DATE:
			return std::string("date");
		case BASE_TYPE_FLOAT:
			return std::string("float");
		case BASE_TYPE_BLOB:
			return std::string("blob");
		default:
			return std::string("");
	}
}

/* 
 * return a list with all the names for base types
 */
std::list<std::string> TypeBaseNames()
{
	std::list<std::string> toReturn;
	
	toReturn.push_back(std::string("int"));
	toReturn.push_back(std::string("char"));
	toReturn.push_back(std::string("decimal"));
	toReturn.push_back(std::string("date"));
	toReturn.push_back(std::string("float"));
	toReturn.push_back(std::string("blob"));

	return toReturn;
}

int TypeBase::GetIdByName(const std::string& name)
{
	if(name == "int")
	{
		return BASE_TYPE_INTEGER;
	}
	else if(name == "char")
	{
		return BASE_TYPE_CHARACTER;
	}
	else if(name == "decimal")
	{
		return BASE_TYPE_DECIMAL;
	}
	else if(name == "date")
	{
		return BASE_TYPE_DATE;
	}
	else if(name == "float")
	{
		return BASE_TYPE_FLOAT;
	}
	else if(name == "blob")
	{
		return BASE_TYPE_BLOB;
	}
	else
		return -1;
}

bool TypeBase::IsFundamental(const std::string& name)
{
	if(name == "int" || name == "char" || name == "date" || 
		name == "float" || name == "blob" || name == "decimal")
		return true;
	else
		return false;
}

/*
 * return a type from a type symbol list
 */
TTypeBasePtr TypeBase::Search(TTypeBasePtrList& list, 
	const std::string& name)
{
	TTypeBasePtrList::const_iterator iterator;
	TTypeBasePtr toReturn;
	
	if(! TypeBase::IsFundamental(name))
	{
		/* this is an aliased type */
		for(iterator = list.begin(); iterator != list.end(); ++iterator)
		{
			if(! (*iterator)->IsBase() && (*iterator)->GetAlias() == name)
			{
				toReturn = *iterator;
				break;
			}
		}
	}
	else
	{
		/* 
		 * this is a fundamental type name, return the first found
		 * which is also a base
		 */
		int id = GetIdByName(name);
		
		for(iterator = list.begin(); iterator != list.end(); ++iterator)
		{
			if((*iterator)->IsBase() && (*iterator)->GetTypeId() == id)
			{
				toReturn = *iterator;
				break;
			}
		}
	}

	/* one needs to check the return value, maybe the type wasn't found */
	return toReturn;	
}

/*
 * return a type from a type symbol list
 */
TTypeBasePtr TypeBase::Search(const TTypeBasePtrList& list, 
	TypeBase& type)
{
	TTypeBasePtrList::const_iterator iterator;
	TTypeBasePtr toReturn;
	
	/* this is an aliased type */
	for(iterator = list.begin(); iterator != list.end(); ++iterator)
	{
		if(type == **iterator)
			toReturn = *iterator;
	}

	return toReturn;
}

void TypeBase::AddIfNotExists(TTypeBasePtrList& list,
	TTypeBasePtr type)
{
	if(! TypeBase::Search(list, *type))
	{
		list.push_back(type);
	}
}

void TypeBase::Add(const AttributeUnit& attribute)
{
	if(! IsBase())
	{
		Attributes::Add(attribute);
	}
	else
	{
		std::cerr << "Warning: you cannot add attributes to a base type. Attribute not bound!" 
			<< std::endl;
	}
}

void TypeBase::CodeGenerator(Backend& backend)
{
	backend.HandleType(*this);
}
