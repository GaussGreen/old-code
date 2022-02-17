/*************************************************************************** 
*  Name:        ddlmaker/src/attribute.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: attribute.cpp,v 1.2 2005/07/06 09:03:02 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include "attribute.h"

class SameKeyAttribute
{
public:
	SameKeyAttribute(const std::string& key)
	:m_key(key)
	{
	}
	bool operator() (const AttributeUnit& attributeUnit)
	{
		return attributeUnit.GetKey() == m_key;
	}
private:
	std::string m_key;
};
	
void Attributes::Add(const AttributeUnit& attribute)
{
	TAttributeUnits::iterator pos;

	pos = find_if(m_attributes->begin(), m_attributes->end(),
			SameKeyAttribute(attribute.GetKey()));
	if(pos != m_attributes->end())
		pos->SetValue(attribute.GetValue());
	else
		m_attributes->push_back(attribute);
}

TVarBasePtr Attributes::GetByValue(const std::string& key)
{
	TVarBasePtr toReturn;
	TAttributeUnits::const_iterator iter = m_attributes->end();

	iter = find_if(m_attributes->begin(), m_attributes->end(),
			SameKeyAttribute(key));
	
	if(iter != m_attributes->end())
	{
		toReturn = iter->GetValue();
	}

	/* in case the requested attribute was not found
	 * it will return a newly constructed boost::shared_ptr
	 * which we can test */
	return toReturn;
}

std::string Attributes::ToString() const
{
	std::list<std::string>::const_iterator knownAIter;
	TAttributeUnits::const_iterator aIter;
	std::string toReturn;

	toReturn += "Attributes\n";
	for(aIter = m_attributes->begin();
		aIter != m_attributes->end();
		aIter++)
	{
		toReturn += aIter->ToString();
		toReturn += "\n";
	}
	toReturn += "\n";

	return toReturn;
}

bool AttributeUnit::operator==(const AttributeUnit& attribute) const
{
	return (m_key == attribute.GetKey() &&
			m_value == attribute.GetValue());
}

AttributeUnit& AttributeUnit::operator=(const AttributeUnit& attribute)
{
	if(!(*this == attribute))
	{
		this->m_key = attribute.GetKey();
		this->m_value = attribute.GetValue();
	}

	return *this;
}

bool Attributes::operator==(const Attributes& attributes) const
{
	return (*m_attributes) == (*attributes.GetAttributes());
}

Attributes& Attributes::operator=(const Attributes& attributes)
{
	TAttributeUnits::const_iterator iterator;

	if(!(*m_attributes == *(attributes.GetAttributes())))
	{
		m_attributes->clear();
		
		for(iterator = attributes.GetAttributes()->begin();
				iterator != attributes.GetAttributes()->end();
				++iterator)
		{
			AttributeUnit temp(*iterator);

			m_attributes->push_back(temp);
		}
	}

	return *this;
}
