//
// C++ Implementation: QMCRNGManager
//
// Description:
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
//
//
#include "edginc/config.hpp"
#include "edginc/QMCRNGManager.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/IQMCAssetRNG.hpp"
#include "edginc/QMCAssetRNG.hpp"
#include "edginc/ran2.h" // SC_ran2Gen

DRLIB_BEGIN_NAMESPACE

QMCRNGManager::QMCRNGManager( MCRandomSP gen, size_t assetsNum ) :
        randoms( NULL ),
        randomGen( gen ),
        uniGen( SC_ran2Gen::create( -1 ) )
{}

void QMCRNGManager::seekToPath( size_t pathIdx )
{
    randomGen->generate( pathIdx );
    randoms = & ( randomGen->getRandomNumbers() ); // returns const DoubleMatrix&

    if ( uniGen )
        uniGen->seekToPath( pathIdx ); // FIXME: simple hack to initialize generators to diff. seeds
}
IQMCAssetRNGSP QMCRNGManager::getAssetRNG( size_t randomIdx )
{
    // FIXME: implement per-asset extra RNG
    return IQMCAssetRNGSP( new QMCAssetRNG( randomGen->getRandomNumbers(), randomIdx, ISuperCubeRNGGenSP() ) );
}
const double * QMCRNGManager::getCorrelatedRandoms( size_t assetIndex, size_t factor )
{
    return & ( ( *randoms ) [ assetIndex + factor ][ 0 ] );
}

IUniformRNGGenSP QMCRNGManager::getSharedGen()
{ ///< returns extra generator of uniform numbers shared between all assets
    return uniGen;
}


DRLIB_END_NAMESPACE
