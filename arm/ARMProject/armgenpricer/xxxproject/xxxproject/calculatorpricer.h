/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file vanillapricer.h
 *  \brief file the pricer of any vanilla product
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#ifndef INXXXPROJECT_VANILLAPRICER
#define INXXXPROJECT_VANILLAPRICER

#include "xxxproject/pricer.h"

class ARM_Model;
class ARM_Security;
class ARM_CapFloor;
class ARM_Digital;
class ARM_Swaption;
class ARM_CorridorLeg;
class ARM_Option;
class ARM_SpreadOption;
class ARM_Object;
class ARM_GenCalculator;

CC_BEGIN_NAMESPACE( ARM )


// All the ARM vanilla products pricer inherit from
// this pricer
class ARM_VanillaPricer : public ARM_Pricer
{
public:

	ARM_VanillaPricer(ARM_Security* security) : itsSecurity(security) {};

	double				PriceWithModel	( ARM_Model* model		);

	ARM_Security*		GetSecurity		( )		const	{ return	itsSecurity;	}

	ARM_Observator*		GetObservator	( )		const	{ return	itsObservator;	}

protected:

	ARM_Security*		itsSecurity;
	ARM_Model*			itsModel;
	ARM_Observator*		itsObservator;
};

// The Cap pricer
class ARM_CAP_Pricer : public ARM_VanillaPricer
{
public:

	ARM_CAP_Pricer				(	ARM_CapFloor* );
	void				SetMkt	(	ARM_MktData * );
	virtual double		Price	(	);
};

// The Osw pricer
class ARM_OSW_Pricer : public ARM_VanillaPricer
{
public:

	ARM_OSW_Pricer				(	ARM_Swaption* );
	void				SetMkt	(	ARM_MktData * );
	virtual double		Price	(	);
};

// The FX option pricer
class ARM_FX_Pricer : public ARM_VanillaPricer
{

public:
	ARM_FX_Pricer				(	ARM_Option*	  );
	void				SetMkt	(	ARM_MktData * );
	virtual double		Price	(	);
};

// The TARNFX calculator pricer
class ARM_TARNFX_Pricer : public ARM_VanillaPricer
{

public:
	ARM_TARNFX_Pricer			(	ARM_GenCalculator*	);
	void				SetMkt	(	ARM_MktData *		);
	virtual double		Price	(	);
};

/* The Swaption pricer
class ARM_SwaptionPricer : public ARM_VanillaPricer
{
public:
	ARM_SwaptionPricer	(		ARM_Swaption* swaption, const ARM_MktData & mktData);
	virtual double Price( );
};

// The Corridor pricer
class ARM_CorridorPricer : public ARM_VanillaPricer
{
public:
	ARM_CorridorPricer(ARM_CorridorLeg* corridorLeg);
	virtual double Price(const ARM_MktData& mktData);
};

// The Digital pricer
class ARM_DigitalPricer : public ARM_VanillaPricer
{
public:
	ARM_DigitalPricer(ARM_Digital* digital);
	virtual double Price(const ARM_MktData& mktData);
};

// The spread option pricer
class ARM_SpreadOptionPricer : public ARM_VanillaPricer
{
public:
	ARM_SpreadOptionPricer(ARM_SpreadOption* spreadOption);
	virtual double Price(const ARM_MktData& mktData);
};
*/
CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/