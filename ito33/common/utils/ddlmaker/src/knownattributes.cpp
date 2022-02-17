/*************************************************************************** 
*  Name:        ddlmaker/src/knownattributes.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Cosmin Cremarenco                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: knownattributes.cpp,v 1.3 2005/08/24 08:21:59 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include <iostream>

#include "knownattributes.h"
#include "fundamentaltypes.h"

boost::shared_ptr<KnownAttributesSingleton>
KnownAttributesSingleton::m_instance;

boost::shared_ptr<KnownAttributesSingleton> 
KnownAttributesSingleton::Instance() 
{
	if(! m_instance)
	{
		m_instance = boost::shared_ptr<KnownAttributesSingleton>
			(new KnownAttributesSingleton);
	}
	return m_instance; 
}

void KnownAttributesSingleton::BuildColumnAttrs()
{
	tableColumnAttrs = boost::shared_ptr<KnownAttributes>(new KnownAttributes());

	tableColumnAttrs->push_back(std::string(DEFAULT_ATTRIBUTE));
	tableColumnAttrs->push_back(std::string(COMMENT_ATTRIBUTE));
}

void KnownAttributesSingleton::BuildTableAttrs()
{
	tableTableAttrs = boost::shared_ptr<KnownAttributes>(new KnownAttributes());

	tableColumnAttrs->push_back(std::string(COMMENT_ATTRIBUTE));
}

void KnownAttributesSingleton::BuildConstraintAttrs()
{
	tableConstraintAttrs = boost::shared_ptr<KnownAttributes>(new KnownAttributes());

	tableConstraintAttrs->push_back(std::string(COMMENT_ATTRIBUTE));
	tableConstraintAttrs->push_back(std::string(ON_DELETE_ATTRIBUTE));
}

void KnownAttributesSingleton::BuildIntegerTypeAttrs()
{
	tableIntegerTypeAttrs = boost::shared_ptr<KnownAttributes>
		(new KnownAttributes());

	tableIntegerTypeAttrs->push_back(std::string(COMMENT_ATTRIBUTE));
	tableIntegerTypeAttrs->push_back(std::string(DEFAULT_ATTRIBUTE));
	tableIntegerTypeAttrs->push_back(std::string(DOMAIN_ATTRIBUTE));
	tableIntegerTypeAttrs->push_back(std::string(SEQUENCE_ATTRIBUTE));
	tableIntegerTypeAttrs->push_back(std::string(DOMAIN_INSTANCE_ATTRIBUTE));
	tableIntegerTypeAttrs->push_back(std::string(NOT_NULL_ATTRIBUTE));
}

void KnownAttributesSingleton::BuildCharacterTypeAttrs()
{
	tableCharacterTypeAttrs = boost::shared_ptr<KnownAttributes>
		(new KnownAttributes());

	tableCharacterTypeAttrs->push_back(std::string(COMMENT_ATTRIBUTE));
	tableCharacterTypeAttrs->push_back(std::string(DEFAULT_ATTRIBUTE));
	tableCharacterTypeAttrs->push_back(std::string(VARYING_ATTRIBUTE));
	tableCharacterTypeAttrs->push_back(std::string(LENGTH_ATTRIBUTE));
	tableCharacterTypeAttrs->push_back(std::string(DOMAIN_ATTRIBUTE));
	tableCharacterTypeAttrs->push_back(std::string(DOMAIN_INSTANCE_ATTRIBUTE));
	tableCharacterTypeAttrs->push_back(std::string(NOT_NULL_ATTRIBUTE));
}

void KnownAttributesSingleton::BuildDecimalTypeAttrs()
{
	tableDecimalTypeAttrs = boost::shared_ptr<KnownAttributes>
		(new KnownAttributes());

	tableDecimalTypeAttrs->push_back(std::string(COMMENT_ATTRIBUTE));
	
	/* the number of numerics after the comma */
	tableDecimalTypeAttrs->push_back(std::string(SCALE_ATTRIBUTE));

	/* the total number of numerics */
	/* TODO: change this to something more meaningful */
	tableDecimalTypeAttrs->push_back(std::string(PRECISION_ATTRIBUTE));
	
	tableDecimalTypeAttrs->push_back(std::string(DEFAULT_ATTRIBUTE));
	tableDecimalTypeAttrs->push_back(std::string(DOMAIN_ATTRIBUTE));
	tableDecimalTypeAttrs->push_back(std::string(DOMAIN_INSTANCE_ATTRIBUTE));
	tableDecimalTypeAttrs->push_back(std::string(NOT_NULL_ATTRIBUTE));
}

void KnownAttributesSingleton::BuildDateTypeAttrs()
{
	tableDateTypeAttrs = boost::shared_ptr<KnownAttributes>
		(new KnownAttributes());

	tableDateTypeAttrs->push_back(std::string(COMMENT_ATTRIBUTE));
	tableDateTypeAttrs->push_back(std::string(DEFAULT_ATTRIBUTE));
	tableDateTypeAttrs->push_back(std::string(CONTAINS_TIME_ATTRIBUTE));
	tableDateTypeAttrs->push_back(std::string(DOMAIN_INSTANCE_ATTRIBUTE));
	tableDateTypeAttrs->push_back(std::string(NOT_NULL_ATTRIBUTE));
}

void KnownAttributesSingleton::BuildFloatTypeAttrs()
{
	tableFloatTypeAttrs = boost::shared_ptr<KnownAttributes>
		(new KnownAttributes());

	tableFloatTypeAttrs->push_back(std::string(COMMENT_ATTRIBUTE));
	tableFloatTypeAttrs->push_back(std::string(DEFAULT_ATTRIBUTE));
	tableFloatTypeAttrs->push_back(std::string(DOMAIN_ATTRIBUTE));
	tableFloatTypeAttrs->push_back(std::string(PRECISION_ATTRIBUTE));
	tableFloatTypeAttrs->push_back(std::string(DOMAIN_INSTANCE_ATTRIBUTE));
	tableFloatTypeAttrs->push_back(std::string(NOT_NULL_ATTRIBUTE));
}

void KnownAttributesSingleton::BuildBlobTypeAttrs()
{
	tableBlobTypeAttrs = boost::shared_ptr<KnownAttributes>
		(new KnownAttributes());

	tableBlobTypeAttrs->push_back(std::string(COMMENT_ATTRIBUTE));
	tableBlobTypeAttrs->push_back(std::string(DEFAULT_ATTRIBUTE));
	tableBlobTypeAttrs->push_back(std::string(DOMAIN_ATTRIBUTE));
	tableBlobTypeAttrs->push_back(std::string(CONTAINS_TEXT_ATTRIBUTE));
	tableBlobTypeAttrs->push_back(std::string(DOMAIN_INSTANCE_ATTRIBUTE));
	tableBlobTypeAttrs->push_back(std::string(NOT_NULL_ATTRIBUTE));
}

bool KnownAttributesSingleton::IsKnownTypeAttribute(
		int typeId, const std::string& attribute)
{
	switch(typeId)
	{
		case BASE_TYPE_INTEGER:
			return find(tableIntegerTypeAttrs->begin(),
					tableIntegerTypeAttrs->end(),
					attribute) != tableIntegerTypeAttrs->end();
		case BASE_TYPE_CHARACTER:
			return find(tableCharacterTypeAttrs->begin(),
					tableCharacterTypeAttrs->end(),
					attribute) != tableCharacterTypeAttrs->end();
		case BASE_TYPE_DECIMAL:
			return find(tableDecimalTypeAttrs->begin(),
					tableDecimalTypeAttrs->end(),
					attribute) != tableDecimalTypeAttrs->end();
		case BASE_TYPE_DATE:
			return find(tableDateTypeAttrs->begin(),
					tableDateTypeAttrs->end(),
					attribute) != tableDateTypeAttrs->end();
		case BASE_TYPE_FLOAT:
			return find(tableFloatTypeAttrs->begin(),
					tableFloatTypeAttrs->end(),
					attribute) != tableFloatTypeAttrs->end();
		case BASE_TYPE_BLOB:
			return find(tableBlobTypeAttrs->begin(),
					tableBlobTypeAttrs->end(),
					attribute) != tableBlobTypeAttrs->end();
		default:
			std::cerr << "Unrecognized type. Aborting!" << std::endl;
			exit(1);
	}
}

bool KnownAttributesSingleton::IsKnownColumnAttribute(
		const std::string& attribute)
{
	return find(tableColumnAttrs->begin(),
			tableColumnAttrs->end(),
			attribute) != tableColumnAttrs->end();
}

bool KnownAttributesSingleton::IsKnownTableAttribute(
		const std::string& attribute)
{
	return find(tableTableAttrs->begin(),
			tableTableAttrs->end(),
			attribute) != tableTableAttrs->end();
}

bool KnownAttributesSingleton::IsKnownConstraintAttribute(
		const std::string& attribute)
{
	return find(tableConstraintAttrs->begin(),
			tableConstraintAttrs->end(),
			attribute) != tableConstraintAttrs->end();
}

#if 0
BUILD_EMPTY_TABLE(Integer)
BUILD_EMPTY_TABLE(Character)
BUILD_EMPTY_TABLE(Decimal)
BUILD_EMPTY_TABLE(Date)
BUILD_EMPTY_TABLE(Float)
BUILD_EMPTY_TABLE(Blob)
#endif
