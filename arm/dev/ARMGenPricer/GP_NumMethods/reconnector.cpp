/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file reconnector.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou
 *	\version 1.0
 *	\date November 2004
 */

#include "gpnummethods/reconnector.h"

#include "gpbase/utilityport.h"  /// for CC_Max

#include "gpnummethods/slice.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class  : ARM_ReconnnectorVar
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_ReconnnectorVar
///	Routine: ReconnectOnTheHedgeDown,ReconnectOnTheHedgeUp,ReconnectOutOfBoundDown,ReconnectOutOfBoundUp
///	Returns: void
///	Action : reconnects the node if necessary
////////////////////////////////////////////////////
void ARM_ReconnectorVar::ReconnectOnTheHedgeDown( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    int downIndex   = nextNodeIdx-1;
    double probaMid = (varRel+meanRel*meanRel - downIndex*downIndex)/(downIndex+nextNodeIdx);
    probaMid        = CC_Max(probaMid,0.0);
    probaDown       = 1.0 - probaMid;
    probaUp         = 0.0;
}

void ARM_ReconnectorVar::ReconnectOnTheHedgeUp( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    int upIndex = nextNodeIdx+1;
    double probaMid = (upIndex*upIndex-varRel-meanRel*meanRel)/(upIndex+nextNodeIdx);
    probaMid        = CC_Max(probaMid,0.0);
    probaUp         = 1.0 - probaMid;
    probaDown       = 0.0;
}

void ARM_ReconnectorVar::ReconnectOutOfBoundDown( int& nextNodeIdx, int minIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    /// Go straight on the edge : mean and variance are not fulfilled
    nextNodeIdx = minIdx;
    probaUp     = 0.0;
    probaDown   = 0.0;
}

void ARM_ReconnectorVar::ReconnectOutOfBoundUp( int& nextNodeIdx, int maxIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    /// Go straight on the edge : mean and variance are not fulfilled
    nextNodeIdx = maxIdx;
    probaUp     = 0.0;
    probaDown   = 0.0;
}

void ARM_ReconnectorVar::ReconnectInside( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    double meanErr = meanRel - nextNodeIdx;
    if(meanErr < 0)
    {
        /// Reconnect as on the lower hedge
        ReconnectOnTheHedgeDown(nextNodeIdx,meanRel,varRel,probaUp,probaDown);
    }
    else
    {
        /// Reconnect as on the upper hedge
        ReconnectOnTheHedgeUp(nextNodeIdx,meanRel,varRel,probaUp,probaDown);
    }

    /// Alternative strategy with only connections to upper & lower nodes
/****
    int downIndex   = nextNodeIdx-1;
    probaUp         = (varRel+meanRel*meanRel - downIndex*downIndex)/(4*nextNodeIdx);
    probaDown       = 1.0 - probaUp;
****/

}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class  : ARM_ReconnectorMean
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_ReconnnectorMean
///	Routine: ReconnectOnTheHedgeDown,ReconnectOnTheHedgeUp,ReconnectOutOfBoundDown,ReconnectOutOfBoundUp
///	Returns: void
///	Action : reconnects the node if necessary
////////////////////////////////////////////////////
void ARM_ReconnectorMean::ReconnectOnTheHedgeDown( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    double meanErr = meanRel - nextNodeIdx;
    if(meanErr > 0)
    {
        /// Strict truncation is possible
        probaUp     = meanErr;
        probaDown   = 0.0;
    }
    else
    {
        /// Lower connection is still used....
        if(probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0)
        {
            ///... but only mean is fulfilled
            probaDown   = -meanErr;
            probaUp     = 0.0;
        }
    }
}

void ARM_ReconnectorMean::ReconnectOnTheHedgeUp( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    double meanErr = meanRel - nextNodeIdx;
    if(meanErr < 0)
    {
        /// Strict truncation is possible
        probaDown   = -meanErr;
        probaUp     = 0.0;
    }
    else
    {
        /// Upper connection is still used....
        if(probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0)
        {
            ///... but only mean is fulfilled
            probaUp     = meanErr;
            probaDown   = 0.0;
        }
    }
}

void ARM_ReconnectorMean::ReconnectOutOfBoundDown( int& nextNodeIdx, int minIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    /// Go straight on the edge : mean and variance are not fulfilled
    nextNodeIdx = minIdx;
    probaUp     = 0.0;
    probaDown   = 0.0;
}

void ARM_ReconnectorMean::ReconnectOutOfBoundUp( int& nextNodeIdx, int maxIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    /// Go straight on the edge : mean and variance are not fulfilled
    nextNodeIdx = maxIdx;
    probaUp     = 0.0;
    probaDown   = 0.0;
}

void ARM_ReconnectorMean::ReconnectInside( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const
{
    double meanErr = meanRel - nextNodeIdx;
    if(meanErr < 0)
    {
        /// Reconnect as on the upper hedge
        probaDown   = -meanErr;
        probaUp     = 0.0;
    }
    else
    {
        /// Reconnect as on the lower hedge
        probaUp     = meanErr;
        probaDown   = 0.0;
    }

    /// Alternative strategy with only connections to upper & lower nodes
/****
    probaUp     = 0.5*(1+meanErr);
    probaDown   = 1 - probaUp;
****/
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/