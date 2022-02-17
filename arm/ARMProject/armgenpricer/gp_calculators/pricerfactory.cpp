/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file pricerfactory.h
 *	\brief pricer factory
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/pricerfactory.h"
#include "gpcalculators/gencalculator.h"

/// gpbase
#include "gpbase/autocleaner.h"


/// gpcalib
#include "gpcalib/vanillapricer.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/calibmethod.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/pricingadviser.h"


/// ARM Kernel
//#include <inst/security.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct   : ARM_PricerFactory
///	Name     : Price (static function)
///	Action   : price instruments
/// Arguments: ARM_Object* sec, ARM_Object* mod
////////////////////////////////////////////////////

CC_NS(std,pair)<bool,double> ARM_PricerFactory::Price( ARM_Object* sec, ARM_Object* mod	)
{
	double price	= 0.0;
	bool canPrice	=false;

	/// Is it a generic model?
    if (( mod != NULL ) && ( ARM_PRICINGMODEL == mod->GetRootName() )) 
	{
		ARM_PricingModel* pmod=dynamic_cast<ARM_PricingModel*>(mod);
		if( !pmod )
			ARM_THROW( ERR_INVALID_ARGUMENT, " dynamic cast into pricing model failed!" );

		/// case of the generic security
		if( ARM_GENSECURITY == sec->GetRootName() )
		{
			ARM_GenSecurity*	genSec  = dynamic_cast<ARM_GenSecurity*>( sec);
			if( !sec )
				ARM_THROW( ERR_INVALID_ARGUMENT, " dynamic cast into generic security failed!" );

			/// creates locally a genpricer
			/// and uses it to price
			ARM_GenPricer* genpricer = new ARM_GenPricer(genSec,pmod);

			/// put an autocleaner on genpricer to ensure proper deletion after use!
			/// very simple smart pointor!
			CC_NS(ARM, ARM_AutoCleaner)<ARM_GenPricer> Hold(genpricer);
			price	= genpricer->Price();
			canPrice= true;
		}
		/// case of a calculator
		else if( ARM_GENCALCULATOR == sec->GetRootName() )
        {
			ARM_GenCalculator* genCalc=dynamic_cast<ARM_GenCalculator*>(sec);
			if( !genCalc)
				ARM_THROW( ERR_INVALID_ARGUMENT, " dynamic cast into generic calculator failed!" );

            /// force the model of the calculator then price
            genCalc->SetPricingModel( ARM_PricingModelPtr((ARM_PricingModel*)pmod->Clone()) );
            price	= genCalc->PriceAndTimeIt();
			canPrice= true;
        }
		/// other case!
        else
		{
			if( ARM_SECURITY == sec->GetRootName()  )
			{
				ARM_Security* security = NULL;//dynamic_cast<ARM_Security*>(sec);
				if( !security )
					ARM_THROW( ERR_INVALID_ARGUMENT, " dynamic cast into security failed!" );

				price	= ARM_VanillaPricer::Price( security, pmod );
				canPrice= true;
			}
			else if( ARM_GP_VANILLAARG == sec->GetRootName() )
			{
				ARM_VanillaArg* vanillaArg = dynamic_cast<ARM_VanillaArg*>(sec);
				if( !vanillaArg )
					ARM_THROW( ERR_INVALID_ARGUMENT, " dynamic cast into vanilla arg failed!" );

				price	= ARM_VanillaPricer::Price( *vanillaArg, pmod );
				canPrice= true;
			}
		}
	}
	else
	{
        if( ARM_GENCALCULATOR == sec->GetRootName() )
        {
			ARM_GenCalculator* genCalc=dynamic_cast<ARM_GenCalculator*>(sec);
			if( !genCalc )
				ARM_THROW( ERR_INVALID_ARGUMENT, " dynamic cast into generic calculator failed!" );
			price = genCalc->PriceAndTimeIt();
			canPrice= true;
		}
		else if( ARM_GENSECURITY == sec->GetRootName()  )
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, " a generic security can only be priced with a generic pricing model!" );
		}
	}

	return CC_NS(std,pair)<bool,double>(canPrice,price);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

