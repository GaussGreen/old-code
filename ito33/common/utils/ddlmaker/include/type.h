/*************************************************************************** 
*  Name:        ddlmaker/include/type.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  $Id: type.h,v 1.2 2005/07/06 09:28:30 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __TYPE_H__
#define __TYPE_H__

#include "fundamentaltypes.h"
#include "attribute.h"
#include "codegeneratorinterface.h"

/* 
 * A brief explanation:
 * The basetypes header contains the types that our frontend
 * will encounter in the parsed script which we want to transform
 * into a usable .sql file.
 * After that, every type has its attributes effectively
 * transforming it in another type.
 * A fundamental type doesn't have any attributes!
 * Attributes are only for types that are aliase, they derive
 * from fundamental types.
 */

class TypeBase;
typedef boost::shared_ptr<TypeBase> TTypeBasePtr;
typedef std::list<TTypeBasePtr> TTypeBasePtrList;

class TypeBase
:public Attributes,
public CodeGeneratorInterface
{
public:
	TypeBase(int typeId)
	:Attributes(), m_typeId(typeId), 
		m_isBase(true), m_alias("")
	{
	}

	TypeBase(int typeId, const std::string alias)
	:Attributes(), m_typeId(typeId), 
		m_isBase(false),
		m_alias(alias)
	{
	}

	TypeBase(const TypeBase& typeBase)
	:Attributes(), m_typeId(typeBase.GetTypeId()),
		m_isBase(typeBase.IsBase()),
		m_alias(typeBase.GetAlias())
	{
		*m_attributes = *(typeBase.GetAttributes());
	}
/*
 * TODO: add a method that would verify that a Var
 * is accepted by the type
 */

	/*
	 * Verify that a certain type is
	 * base or a derived type
	 */
	bool IsBase() const
	{
		return m_isBase;
	}

	std::string GetAlias() const
	{
		return m_alias; 
	} 

	void SetIsBase(bool isBase) 
	{
		m_isBase = isBase;
	}

	void SetAlias(const std::string alias)
	{
		m_isBase = false;
		m_alias = alias;
	}
	
	virtual void Add(const AttributeUnit& attribute);
	
	/*
	 * construct a list of pointer towards all the base types
	 */
	static TTypeBasePtrList BuildTypeBaseTable();

	/* 
	 * return a list with all the names for base types
	 */
	static std::list<std::string> TypeBaseNames();
	
	/*
	 * return a type from a type symbol list
	 */
	static TTypeBasePtr Search(TTypeBasePtrList& list, 
			const std::string& name);

	/*
	 * return a type from a type symbol list
	 */
	static TTypeBasePtr Search(const TTypeBasePtrList& list, 
			TypeBase& type);

	static void AddIfNotExists(TTypeBasePtrList& list,
			TTypeBasePtr type);

	/*
	 * return the name of a fundamental type given its id 
	 */
	static std::string GetNameById(int id);

	/*
	 * return the id of a fundamental type given its name
	 */
	static int GetIdByName(const std::string& name);
	
	/*
	 * Determine if a type name is fundamental or not
	 */
	static bool IsFundamental(const std::string& name);

	bool operator==(const TypeBase& type) const
	{
		return m_typeId == type.GetTypeId() &&
			Attributes::operator==(type) &&
			GetAlias() == type.GetAlias() &&
			IsBase() == type.IsBase();
	}
	
	virtual ~TypeBase()
	{
	}

	int GetTypeId() const
	{
		return m_typeId;
	}

	void CodeGenerator(Backend& backend);
	
private:
	/*
	 * Holds the id of the base type.
	 * The base types are described in an enum in fundamentaltypes.h
	 */
	int m_typeId;

	bool m_isBase;

	std::string m_alias;
};
#endif
