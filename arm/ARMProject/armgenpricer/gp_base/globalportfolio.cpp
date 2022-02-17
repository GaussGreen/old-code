
#include "gpbase/globalportfolio.h"
//#include "gpcalculators/pricerfactory.h"

/// gpbase
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/checkarg.h"
#include "gpbase/gpvector.h"

CC_USING_NS(ARM_Check,CheckSameArgSize);

CC_BEGIN_NAMESPACE( ARM )


ARM_GlobalPortfolio::ARM_PortfolioElem::ARM_PortfolioElem(const ARM_PortfolioElem& rhs)
:
	itsAsset( const_cast<ARM_Object*>(CreateClone(rhs.itsAsset))), 
	itsWeight(rhs.itsWeight), 
	itsPrice(rhs.itsPrice) 
{
}

ARM_GlobalPortfolio::ARM_PortfolioElem::ARM_PortfolioElem( const ARM_Object& asset, double weight, double price )
:	
	itsAsset( const_cast<ARM_Object*>(CreateClone(&asset))), 
	itsWeight(weight), 
	itsPrice(price) 
{
}

ARM_GlobalPortfolio::ARM_PortfolioElem::~ARM_PortfolioElem()
{
	delete itsAsset;
	itsAsset = NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_GlobalPortfolio
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_GlobalPortfolio::toString(const string& indent, const string& nextIndent) const
{
	return "ARM_GlobalPortfolio";
}


////////////////////////////////////////////////////
///	Class  : ARM_GlobalPortfolio
///	Routine: Constructor,Copy Constructor, Destructor, Assignment operator
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////

ARM_GlobalPortfolio::ARM_GlobalPortfolio( 
	const vector<ARM_Object*>& assets,
    const std::vector<double>& weights,
    const std::vector<double>& mktprices)
:	
	ARM_RootObject(),	
	itsElems( assets.size() )
{
	/// validation
	CheckSameArgSize(  assets, weights,		"assets", "weights"		);
	CheckSameArgSize(  assets, mktprices,	"assets", "mktprices"	);

	size_t nbElems = assets.size();
	for( size_t i=0; i<nbElems; ++i )
	{
		itsElems[i] = new ARM_PortfolioElem( *assets[i], weights[i], mktprices[i] );
	}
}

	////////////////////////////////////////////////////
///	Class  : ARM_GlobalPortfolio
///	Routine: copy sonstructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_GlobalPortfolio::ARM_GlobalPortfolio( const ARM_GlobalPortfolio& rhs )
:	
	ARM_RootObject(rhs), 
	itsElems( rhs.itsElems.size() )
{	
	DuplicateCloneablePointorVectorInPlace< ARM_PortfolioElem >( rhs.itsElems, itsElems );
}

ARM_GlobalPortfolio::~ARM_GlobalPortfolio()
{
	DeletePointorVector<ARM_GlobalPortfolio::ARM_PortfolioElem>( itsElems );
}

////////////////////////////////////////////////////
///	Class  : ARM_GlobalPortfolio
///	Routine: HasSameAssets
///	Returns: 
///	Action : tells whether the portfolio is composed of the same type of asset
////////////////////////////////////////////////////

bool ARM_GlobalPortfolio::HasSameAssets() const
{
	if( !itsElems.empty())
	{
		ARM_CLASS_NAME type = itsElems[0]->GetAsset()->GetRootName();
		for( size_t i=1; i<itsElems.size(); ++i )
			if( type != itsElems[i]->GetAsset()->GetRootName() )
				return false;
	}

	return true;
}

////////////////////////////////////////////////////
///	Class  : ARM_GlobalPortfolio
///	Routine: Price
///	Returns: double 
///	Action : computes the price of the portfolio
////////////////////////////////////////////////////
double ARM_GlobalPortfolio::Price(ARM_Object* model ) const
{
	double price=0.0;
	std::pair<bool,double> pricingResult(true, 0.0);
	for( const_iterator iter=begin(); iter!=end(); ++iter )
	{
		//pricingResult = ARM_PricerFactory::Price((*iter)->GetAsset(), model);

		if(pricingResult.first)
			price += pricingResult.second * (*iter)->GetWeight();
	}

	return price;
}

////////////////////////////////////////////////////
///	Class  : ARM_GlobalPortfolio
///	Routine: GetMktPrice
///	Returns: double 
///	Action : computes the mkt price of the portfolio
////////////////////////////////////////////////////

double ARM_GlobalPortfolio::GetMktPrice() const 
{
	double mktPrice=0.0;
	for( const_iterator iter=begin(); iter!=end(); ++iter )
	{
		mktPrice += (*iter)->GetPrice() * (*iter)->GetWeight();
	}

	return mktPrice;
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

