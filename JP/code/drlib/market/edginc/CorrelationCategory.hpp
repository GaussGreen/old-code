//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrelationCategory.hpp 
//
//   Description : CorrelationCategory
//
//   Author      : Eva X Strasser
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_CORRELATIONCATEGORY_HPP
#define EDR_CORRELATIONCATEGORY_HPP
#include "edginc/MarketObject.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL CorrelationCategory: public MarketObject {
public:
    static CClassConstSP const TYPE;

    static const string SKIP_FWD_CORRELATION;

    friend class CorrelationCategoryHelper;

    virtual ~CorrelationCategory();

    virtual string getName() const; // mandatory

    virtual string getCategoryName() const;

	virtual string getRegionName() const;	

    CorrelationCategory();

private:
    string              assetName;
    string              categoryName;
	string				regionName;
};

typedef smartPtr<CorrelationCategory> CorrelationCategorySP;


DRLIB_END_NAMESPACE
#endif
