//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ICDOQuotesGenerator.hpp
//
//   Description : ICDOQuotesGenerator is an interface for
//                  generators of pseudo CDO quotes to be used
//                  afterwards in a calibration.
//
//   Author      : Sebastien Gay    
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------


#ifndef I_CDO_QUOTES_GENERATOR_HPP
#define I_CDO_QUOTES_GENERATOR_HPP

#include "edginc/IModel.hpp"
#include "edginc/MarketWrapper.hpp"

DRLIB_BEGIN_NAMESPACE

class CDOQuotes;
typedef smartConstPtr<CDOQuotes> CDOQuotesConstSP;

/**
 * Quote generator is a builder of CDOquotes
 * used whenever the market quotes are unavailable.
 * The CdoQuotes generator will generate the cdo quotes from 
 * the base correlation surface in the market environment
 * */
class PRODUCTS_DLL ICDOQuotesGenerator : public virtual IObject {
public:
	
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~ICDOQuotesGenerator();

    /** Main method : builds the CDOquotes object */
    virtual CDOQuotesConstSP buildCDOQuotes() const = 0;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

DECLARE(ICDOQuotesGenerator);

typedef MarketWrapper<ICDOQuotesGenerator> ICDOQuotesGeneratorWrapper;

DRLIB_END_NAMESPACE

#endif /* I_CDO_QUOTES_GENERATOR_HPP*/
