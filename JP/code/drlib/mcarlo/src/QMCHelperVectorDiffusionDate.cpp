//
// C++ Implementation: QMCHelperVectorDiffusionDate
//
// Description:
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
//
//

#include "edginc/config.hpp"
#include "edginc/QMCHelperVectorDiffusionDate.hpp"

DRLIB_BEGIN_NAMESPACE


void QMCHelperVectorDiffusionDate::updateDependent( const DateTime& maturity ) const
{ // call updateMaturity on dependent
    for ( size_t i = 0; i < dependent.size(); ++i )
        if ( dependent[ i ] )
            dependent[ i ] ->updateMaturity( maturity );
}

QMCHelperVectorDiffusionDate::QMCHelperVectorDiffusionDate() : dependent()
{
} // zero dependent : IR

QMCHelperVectorDiffusionDate::QMCHelperVectorDiffusionDate( IQMCHelperMaxDiffusionDate *p ) : dependent()
{
    if ( p )
        dependent.push_back( p );
} // 1 dependent CR, EQ

QMCHelperVectorDiffusionDate::QMCHelperVectorDiffusionDate( IQMCHelperMaxDiffusionDate *p1, IQMCHelperMaxDiffusionDate *p2 ) : dependent()
{ // 2 dependent FX
    if ( p1 )
        dependent.push_back( p1 );
    if ( p2 )
        dependent.push_back( p2 );
}
/** append a new dependent to the list of dates to be notified */
void QMCHelperVectorDiffusionDate::add
    ( IQMCHelperMaxDiffusionDate *dep )
{
    if ( dep != NULL && find( dependent.begin(), dependent.end(), dep ) == dependent.end() ) {
        dep->updateMaturity( currentMaturity() );
        dependent.push_back( dep );
    }
}

DRLIB_END_NAMESPACE
