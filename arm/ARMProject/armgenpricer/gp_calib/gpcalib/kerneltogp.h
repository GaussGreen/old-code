/*§
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file kerneltogp.h
 *
 *  \brief conversion routines from kernel securities to gp securities
 *
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_KERNELTOGP_H
#define _INGPCALIB_KERNELTOGP_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"

/// forward declaration
/// kernel object in the global namespace!
class ARM_Security;
class ARM_Swaption;
class ARM_CapFloor;
class ARM_Digital;
class ARM_CorridorLeg;
class ARM_Option;
class ARM_SpreadOption;
class ARM_SwapLeg;
class ARM_Portfolio;
class ARM_SumOpt;
class ARM_SmiledSwaption;
class ARM_OptionPortfolio;
class ARM_StripOption;
class ARM_StripDigitalOption;

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_InfLeg;
class ARM_InfCapFloor;

struct ARM_VanillaArg;
struct ARM_VanillaSwaptionArg;
struct ARM_VanillaCapArg;
struct ARM_VanillaDigitalArg;
struct ARM_VanillaCorridorLegArg;
struct ARM_VanillaEqOption;
struct ARM_VanillaFxOption;
struct ARM_VanillaSpreadOptionArg;
struct ARM_VanillaIRSwaplegArg;
struct ARM_VanillaPortfolio;
struct ARM_VanillaInfSwaplegArg;
struct ARM_VanillaInfCapArg;
struct ARM_VanillaSumOptArg;
struct ARM_VanillaSmiledSwaption;
struct ARM_VanillaIrFxSwaption;
struct ARM_VanillaStripArg;

struct ARM_ConverterFromKernel
{
	/// general function to convert an arm kernel function to an ARM_VanillaArg
	static ARM_VanillaArg* ConvertSecuritytoArgObject( ARM_Security* Security, double asOfDate, const string& modelName="" );
	static ARM_VanillaPortfolio* ConvertPortoflio( ARM_Portfolio* port, double asOfDate, const string& modelName="" );
	static void TransferToSecurity(ARM_VanillaArg* vanillaArg, ARM_Security* sec);
	
private:
	friend class ARM_SFRM;
	/// --------- specific functions
	static ARM_VanillaSwaptionArg* ConvertVanillaSwaption       ( ARM_Swaption*    swaption,   double asOfDate, const string& modelName="");
	static ARM_VanillaCapArg* ConvertVanillaCapFloor			( ARM_CapFloor*    capFloor,   double asOfDate, const string& modelName="" );
    static ARM_VanillaDigitalArg* ConvertVanillaDigital			( ARM_Digital*     digital,    double asOfDate, const string& modelName="" );
    static ARM_VanillaCorridorLegArg* ConvertVanillaCorridorLeg ( ARM_CorridorLeg* corridorleg,double asOfDate, const string& modelName="" );
	static ARM_VanillaEqOption* ConvertToEqOption				( ARM_Option* option, double asOfDate, const string& modelName = "" );
	static ARM_VanillaFxOption* ConvertToFXOption				( ARM_Option* option, double asOfDate);
	static ARM_VanillaStripArg* ConvertToStripOption			( ARM_StripOption* option, double asOfDate, const string& modelName = "");
	static ARM_VanillaStripArg* ConvertToStripDigitalOption		( ARM_StripDigitalOption* option, double asOfDate, const string& modelName = "");
	static ARM_VanillaIRSwaplegArg* ConvertVanillaIRSwapleg		( ARM_SwapLeg*    swapleg,   double asOfDate, const string& modelName="" );
	static ARM_VanillaSpreadOptionArg* ConvertVanillaSpreadOption(	ARM_SpreadOption*	spreadoption,   double asOfDate, const string& modelName="");
	static ARM_VanillaInfSwaplegArg* ConvertVanillaInfSwapLeg	( ARM_InfLeg* infleg, double asOfDate, const string& modelName="");
	static ARM_VanillaInfCapArg* ConvertVanillaInfCapFloor		( ARM_InfCapFloor* infCapFloor, double asOfDate, const string& modelName="");
	static ARM_VanillaSumOptArg* ConvertVanillaSumOpt			( ARM_SumOpt* sumOpt, double asOfDate, const string& modelName="");
	static ARM_VanillaSmiledSwaption* ConvertVanillaSmiledSwaption	( ARM_SmiledSwaption* opt,double asOfDate, const string& modelName="");
	static ARM_VanillaIrFxSwaption* ConvertIrFxSwaption			(ARM_OptionPortfolio* irFxSwaption,double asOfDate, const string& modelName="");

	static void TranferToSumOptSec(ARM_VanillaSumOptArg* sumOptArg, ARM_SumOpt* sumOptSec);
};




CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
