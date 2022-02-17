#ifndef __KNOWN_ATTRIBUTES_H__
#define __KNOWN_ATTRIBUTES_H__

#include <list>
#include <boost/shared_ptr.hpp>
#include <string>

typedef std::list<std::string> KnownAttributes;
typedef boost::shared_ptr<KnownAttributes> KnownAttributesPtr;

#define GET_ATTR_TYPE(type) \
	KnownAttributesPtr Get ## type ## Attributes() const \
	{ \
		return table ## type ## TypeAttrs; \
	}
/*
 * We use this class in order to build known attributes
 */
class KnownAttributesSingleton
{
public:
	static boost::shared_ptr<KnownAttributesSingleton>
		Instance();
	
	KnownAttributesPtr GetColumnAttributes() const
	{ 
		return tableColumnAttrs;
	}

	KnownAttributesPtr GetTableAttributes() const
	{ 
		return tableTableAttrs;
	}

	GET_ATTR_TYPE(Integer)
	GET_ATTR_TYPE(Character)
	GET_ATTR_TYPE(Decimal)
	GET_ATTR_TYPE(Float)
	GET_ATTR_TYPE(Date)
	GET_ATTR_TYPE(Blob)

	bool IsKnownTypeAttribute(int typeId, const std::string& attribute);
	
	bool IsKnownColumnAttribute(const std::string& attribute);
	
	bool IsKnownTableAttribute(const std::string& attribute);

	bool IsKnownConstraintAttribute(const std::string& attribute);
protected:
	KnownAttributesSingleton()
	{
		BuildColumnAttrs();
		BuildTableAttrs();
		BuildConstraintAttrs();
		BuildIntegerTypeAttrs();
		BuildCharacterTypeAttrs();
		BuildDecimalTypeAttrs();
		BuildFloatTypeAttrs();
		BuildDateTypeAttrs();
		BuildBlobTypeAttrs();
	}

private:
	void BuildColumnAttrs();
	void BuildTableAttrs();
	void BuildConstraintAttrs();
	
	void BuildIntegerTypeAttrs();
	void BuildCharacterTypeAttrs();
	void BuildDecimalTypeAttrs();
	void BuildFloatTypeAttrs();
	void BuildDateTypeAttrs();
	void BuildBlobTypeAttrs();
	
	KnownAttributesPtr tableColumnAttrs;
	KnownAttributesPtr tableTableAttrs;
	KnownAttributesPtr tableConstraintAttrs;

	/*
	 * Known attributes for fundamental types
	 */
	KnownAttributesPtr tableIntegerTypeAttrs;
	KnownAttributesPtr tableCharacterTypeAttrs;
	KnownAttributesPtr tableDecimalTypeAttrs;
	KnownAttributesPtr tableDateTypeAttrs;
	KnownAttributesPtr tableFloatTypeAttrs;
	KnownAttributesPtr tableBlobTypeAttrs;

	/*
	 * instance of this singleton
	 */
	static boost::shared_ptr<KnownAttributesSingleton> m_instance;

	/* make the singleton not copyable */
	KnownAttributesSingleton(const KnownAttributesSingleton&);
	KnownAttributesSingleton& operator=(const KnownAttributesSingleton&);
};

#define BUILD_EMPTY_TABLE(type) \
	void KnownAttributesSingleton::Build ## type ## TypeAttrs() \
	{ \
		table ## type ## TypeAttrs = boost::shared_ptr<KnownAttributes> \
			(new KnownAttributes()); \
	}

/*
 * The known attributes
 */
#define COMMENT_ATTRIBUTE "comment"
#define VARYING_ATTRIBUTE "varying"
#define LENGTH_ATTRIBUTE "length"
#define SEQUENCE_ATTRIBUTE "sequence"
#define DOMAIN_ATTRIBUTE "domain"
#define DEFAULT_ATTRIBUTE "default"
#define DECIMALS_ATTRIBUTE "decimals"
#define PRECISION_ATTRIBUTE "precision"
#define SCALE_ATTRIBUTE "scale"
#define NOT_NULL_ATTRIBUTE "notNull"
#define CONTAINS_TIME_ATTRIBUTE "containsTime"
#define CONTAINS_TEXT_ATTRIBUTE "containsText"
#define DOMAIN_INSTANCE_ATTRIBUTE "domainInstance"
#define ON_DELETE_ATTRIBUTE "onDelete"
#define NOT_NULL_ATTRIBUTE "notNull"

#define TYPE_ATTRIBUTE "type"

#endif
