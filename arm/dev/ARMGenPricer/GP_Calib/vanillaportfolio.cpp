/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillaportfolio.cpp
 *
 *  \brief vanilla portfolio
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */



#include "gpcalib/vanillaportfolio.h"
#include "gpcalib/vanillaarg.h"

/// gpbase
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/checkarg.h"
CC_USING_NS(ARM_Check,CheckSameArgSize);

CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_VanillaPortfolio
///	Routine: GetMktPrice and Clone
///	Returns: double
///	Action : return the mkt price of a given security
////////////////////////////////////////////////////

double ARM_VanillaPortfolio::ARM_PortfolioElem::GetMktPrice( ) const
{
	return itsAsset->GetMktPrice();
}


ARM_VanillaPortfolio::ARM_PortfolioElem* ARM_VanillaPortfolio::ARM_PortfolioElem::Clone() const
{
	return new ARM_PortfolioElem( ARM_VanillaArgPtr( static_cast<ARM_VanillaArg*>(  itsAsset->Clone() ) ), itsWeight, itsPrecision );
}



////////////////////////////////////////////////////
///	Class  : ARM_VanillaPortfolio
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

string ARM_VanillaPortfolio::toString(const string& indent, const string& nextIndent) const
{
	return "ARM_VanillaPortfolio";
}


////////////////////////////////////////////////////
///	Class  : ARM_VanillaPortfolio
///	Routine: Constructor,Copy Constructor, Destructor, Assignment operator
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////

ARM_VanillaPortfolio::ARM_VanillaPortfolio( 
	const vector<ARM_VanillaArgPtr>& assets,
    const ARM_GP_Vector& weights,
    const ARM_GP_Vector& mktprices,
    const ARM_GP_Vector& precisions )
:	
	ARM_VanillaArg(),	
	itsElems( assets.size() )
{
	/// validation
	CheckSameArgSize(  assets, weights,		"assets", "weights"		);
	CheckSameArgSize(  assets, mktprices,	"assets", "mktprices"	);
	CheckSameArgSize(  assets, precisions,	"assets", "precisions"	);

	size_t nbElems = assets.size();
	for( size_t i=0; i<nbElems; ++i )
	{
		assets[i]->SetMktPrice( mktprices[i] );
		itsElems[i] = new ARM_PortfolioElem( assets[i], weights[i], precisions[i] );
	}
}


ARM_VanillaPortfolio::ARM_VanillaPortfolio( const ARM_VanillaPortfolio& rhs )
:	ARM_VanillaArg(rhs), itsElems( rhs.itsElems.size() )
{	DuplicateCloneablePointorVectorInPlace< ARM_PortfolioElem >( rhs.itsElems, itsElems ); }

ARM_VanillaPortfolio& ARM_VanillaPortfolio::operator=( const ARM_VanillaPortfolio& rhs )
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator =( rhs );
		DeletePointorVector<ARM_VanillaPortfolio::ARM_PortfolioElem>( itsElems );
		DuplicateCloneablePointorVectorInPlace<ARM_VanillaPortfolio::ARM_PortfolioElem>( rhs.itsElems, itsElems );
	}
	return *this;
}

ARM_VanillaPortfolio::~ARM_VanillaPortfolio()
{
	DeletePointorVector<ARM_VanillaPortfolio::ARM_PortfolioElem>( itsElems );
}



////////////////////////////////////////////////////
///	Class  : ARM_VanillaPortfolio
///	Routine: HasSameAssets
///	Returns: 
///	Action : tells whether the portfolio is composed of the same type of asset
////////////////////////////////////////////////////

bool ARM_VanillaPortfolio::HasSameAssets() const
{
	if( !itsElems.empty())
	{
		ARM_VanillaArg::VanillaArgType type = itsElems[0]->GetAsset()->GetType();
		for( size_t i=1; i<itsElems.size(); ++i )
			if( type != itsElems[i]->GetAsset()->GetType() )
				return false;
	}

	return true;
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaPortfolio
///	Routine: AreAssetsGrowingByExpiry
///	Returns: 
///	Action : tells whether assets are growing by expiry
////////////////////////////////////////////////////

bool ARM_VanillaPortfolio::AreAssetsGrowingByExpiry() const
{
	if( !itsElems.empty())
	{
		double firstExpiry = itsElems[0]->GetAsset()->GetExpiry();
		for( size_t i=1; i<itsElems.size(); ++i )
			if( firstExpiry != itsElems[i]->GetAsset()->GetExpiry() )
				return false;
	}
	return true;
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaPortfolio
///	Routine: GetSubPortfolio
///	Returns: ARM_VanillaPortfolio*
///	Action : get the sub portfolio defined by two indexes!
////////////////////////////////////////////////////

ARM_VanillaPortfolio* ARM_VanillaPortfolio::GetSubPortfolio( size_t begin, size_t end) const
{
	if( end<begin)
		CC_NS(std,swap)(begin,end);

	if( end>itsElems.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "out of range!" );

	vector<ARM_PortfolioElem*> pfElems( end-begin);

	for( size_t i=begin, j=0; i<end; ++i, ++j )
		pfElems[j] = itsElems[i]->Clone();

	return new ARM_VanillaPortfolio( pfElems );
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaPortfolio
///	Routine: BuiltPortfolioByMaturity
///	Returns: vector<ARM_VanillaPortfolio*>
///	Action : creates the vector of portfolio split by maturity
///			WARNING ASSUMES THAT THE PORTFOLIO IS SORTED BY EXPIRY!
////////////////////////////////////////////////////

vector<ARM_VanillaPortfolio*> ARM_VanillaPortfolio::SplitPortfolioByMaturity() const
{
	/// first collect the various maturities
    vector<int> IndexChanges;

    size_t i,size1,totalSize = itsElems.size();
    for(i=1, size1=1; i<totalSize;++i )
    {
        double expiry = itsElems[i-1]->GetAsset()->GetExpiry();
        while(i<totalSize && itsElems[i]->GetAsset()->GetExpiry()< expiry+ARM_LAG_THRESHOLD )
        {
            ++i;
            ++size1;
        }
        IndexChanges.push_back(size1);
    }

    size_t PfNb = IndexChanges.size();
    vector<ARM_VanillaPortfolio*> pfvector(PfNb);

    size_t begin = 0, end;

    for(i=0, end = IndexChanges[i]; i<PfNb; ++i, begin = end )
        pfvector[i] = GetSubPortfolio(begin,end);

    return pfvector;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaPortfolio
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaPortfolio::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaPortfolio
///	Routine: Price
///	Returns: double 
///	Action : computes the price of the portfolio
////////////////////////////////////////////////////

double ARM_VanillaPortfolio::Price(ARM_PricingModel* model ) const
{
	double price=0.0;
	for( const_iterator iter=begin(); iter!=end(); ++iter )
		price += (*iter)->GetAsset()->Price(model) * (*iter)->GetWeight();
	return price;
}


////////////////////////////////////////////////////
///	Class  : ARM_VanillaPortfolio
///	Routine: GetMktPrice
///	Returns: double 
///	Action : computes the mkt price of the portfolio
////////////////////////////////////////////////////

double ARM_VanillaPortfolio::GetMktPrice() const 
{
	double mktPrice=0.0;
	for( const_iterator iter=begin(); iter!=end(); ++iter )
		mktPrice += (*iter)->GetAsset()->GetMktPrice() * (*iter)->GetWeight();
	return mktPrice;
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

