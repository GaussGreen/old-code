/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file reconnector.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_RECONNECTOR_H
#define _INGPNUMMETHODS_RECONNECTOR_H

#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/env.h"
#include "gpbase/typedef.h"


CC_BEGIN_NAMESPACE( ARM )

struct ARM_ReconnectorBase : public ARM_RootObject
{
    enum ReconnectorType
	{
		DoNothing = 0,
        Mean,
		Variance
    };

    virtual void ReconnectOutOfBoundDown( int& nextNodeIdx, int minIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const = 0;
    virtual void ReconnectOutOfBoundUp( int& nextNodeIdx, int maxIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const = 0;
    virtual void ReconnectOnTheHedgeDown( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const = 0;
    virtual void ReconnectOnTheHedgeUp( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const = 0;
    virtual void ReconnectInside( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const = 0;
};

struct ARM_ReconnectorDoNothing : public ARM_ReconnectorBase
{
    virtual void ReconnectOutOfBoundDown( int& nextNodeIdx, int minIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const {};
    virtual void ReconnectOutOfBoundUp( int& nextNodeIdx, int maxIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const {};
    virtual void ReconnectOnTheHedgeDown( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const {};
    virtual void ReconnectOnTheHedgeUp( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const {};
    virtual void ReconnectInside( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const {};

    /// ARM_Object support
    virtual ARM_Object* Clone() const  { return new ARM_ReconnectorDoNothing(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_ReconnectorDoNothing"; }
};

struct ARM_ReconnectorVar : public ARM_ReconnectorBase
{
    virtual void ReconnectOutOfBoundDown( int& nextNodeIdx, int minIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;
    virtual void ReconnectOutOfBoundUp( int& nextNodeIdx, int maxIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;
    virtual void ReconnectOnTheHedgeDown( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;
    virtual void ReconnectOnTheHedgeUp( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;
    virtual void ReconnectInside( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;

    /// ARM_Object support
    virtual ARM_Object* Clone() const  { return new ARM_ReconnectorVar(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_ReconnectorVar"; }
};


struct ARM_ReconnectorMean : public ARM_ReconnectorBase
{
    virtual void ReconnectOutOfBoundDown( int& nextNodeIdx, int minIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;
    virtual void ReconnectOutOfBoundUp( int& nextNodeIdx, int maxIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;
    virtual void ReconnectOnTheHedgeDown( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;
    virtual void ReconnectOnTheHedgeUp( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;
    virtual void ReconnectInside( int& nextNodeIdx, double meanRel, double varRel, double& probaUp, double& probaDown ) const;

    /// ARM_Object support
    virtual ARM_Object* Clone() const  { return new ARM_ReconnectorMean(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_ReconnectorMean"; }
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

