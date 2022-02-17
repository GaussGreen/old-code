/*************************************************************************** 
*  Name:        ddlmaker/include/attribute.h                                 
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: attribute.h,v 1.3 2005/08/24 08:21:52 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#ifndef __ATTRIBUTE_H__
#define __ATTRIBUTE_H__

#include <map>
#include <list>
#include <boost/shared_ptr.hpp>

#include "variable.h"

class AttributeUnit
{
public:
	AttributeUnit(const AttributeUnit& attribute)
	{
		m_key = attribute.GetKey();
		m_value = attribute.GetValue()->Clone();
	}
	
	AttributeUnit(std::string key, TVarBasePtr value)
	:m_key(key), m_value(value)
	{
	}

	std::string GetKey() const
	{
		return m_key;
	}

	TVarBasePtr GetValue() const
	{
		return m_value;
	}

	void SetValue(TVarBasePtr value)
	{
		m_value = value;
	}

	bool operator==(const AttributeUnit& attribute) const;
	
	AttributeUnit& operator=(const AttributeUnit& attribute);

	/*
	 * Use this for debugging purposes
	 */
	std::string ToString() const
	{
		return m_key + "=" + m_value->ToString();
	}
private:
	std::string m_key;
	TVarBasePtr m_value;
};

typedef std::list<AttributeUnit> TAttributeUnits;
typedef boost::shared_ptr<TAttributeUnits> TAttributeUnitsPtr;

/*
 * This class should be inherited by all
 * classes that can have attributes (i.e. the columns, types, columns..)
 */
class Attributes
{
public:
	Attributes()
	:m_attributes(new TAttributeUnits)
	{
	}

	virtual void Add(const AttributeUnit& attribute);
	
	TVarBasePtr GetByValue(const std::string& key);

	bool operator==(const Attributes& attributes) const;

	/*
	 * Use this for debuggin purposes
	 */
	std::string ToString() const;
	
	TAttributeUnitsPtr GetAttributes() const
	{
		return m_attributes;
	}

	Attributes& operator=(const Attributes& attributes);

	virtual ~Attributes()
	{
	}

protected:
	/* 
	 * this holds the pointer to the actual
	 * container of attributes
	 */
	TAttributeUnitsPtr m_attributes;
};

#define GET_ATTR(type, name_char) (type).GetByValue(std::string(name_char))

#endif
