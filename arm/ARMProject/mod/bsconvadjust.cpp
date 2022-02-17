/*
 *
 * Copyright (c) CDC IXIS CM September 2003 Paris
 *
 *
 * $Log: bsconvadjust.cpp,v $
 * Revision 1.4  2004/04/21 16:01:06  mcampet
 * MC case Start>Pay improved
 *
 * Revision 1.3  2004/03/30 10:10:36  rguillemot
 * Replication bug fix
 *
 * Revision 1.2  2004/03/24 13:44:38  rguillemot
 * SUMMIT Payment Lag
 *
 * Revision 1.28  2004/02/25 13:31:51  emezzine
 * Correct a bug (ComputeConxAdhVolatility).
 *
 * Revision 1.27  2004/02/25 11:05:25  mcampet
 * MC Coorect Pay<End Case
 *
 * Revision 1.26  2004/02/09 08:52:48  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.25  2004/02/05 08:51:19  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.24  2004/02/04 15:16:13  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.23  2004/02/02 17:10:16  rguillemot
 * Adjust libor bug fix
 *
 * Revision 1.22  2004/01/23 15:08:41  mcampet
 * MC correct R1tstp in CMS adj
 *
 * Revision 1.21  2004/01/23 14:16:07  mcampet
 * MC Modif TimeLagAdjustCMS
 *
 * Revision 1.20  2004/01/22 08:24:37  rguillemot
 * Replication Convexity Adjustment (Merge BS_Model)
 *
 * Revision 1.19  2004/01/21 11:38:44  emezzine
 * Improvement and adde exception in ConvadjustLibor
 *
 * Revision 1.16  2004/01/15 16:12:28  mcampet
 * MC modif ConvAdjustLibor
 *
 * Revision 1.15  2004/01/15 12:30:23  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.14  2004/01/15 09:07:01  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.13  2004/01/13 16:27:46  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.12  2003/12/29 10:38:16  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.11  2003/12/29 07:38:50  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.10  2003/12/23 16:55:20  rguillemot
 * Conv Adj Merge
 *
 * Revision 1.9  2003/12/22 09:52:11  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.8  2003/12/18 08:45:16  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.7  2003/12/11 14:59:20  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.6  2003/12/10 15:08:22  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.5  2003/12/09 08:30:15  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.4  2003/12/04 14:42:10  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.3  2003/11/28 07:04:51  rguillemot
 * Transfert des ajustements de convexite inflation
 *
 * Revision 1.2  2003/11/03 17:27:09  emezzine
 *  conv2unix
 *
 * Revision 1.2  2003/10/22 13:07:47  emezzine
 * First version
 *
 */
/*----------------------------------------------------------------------------*/
/*! \file bsconvadjust.cpp
 *
 *  \brief B&S Convexity adjustment 
 *
 *	\author  El Mostafa EZZINE
 *	\version 1.0
 *	\date October 2003
 */
/*----------------------------------------------------------------------------*/

/*! \class   ARM_BSConvAdjust
 *	\brief  object that manage convexity and payment lag for alla model
 *	
 *	\author  El Mostafa EZZINE
 *	\version 1.0
 */

#include "bsconvadjust.h"
#include "bsmodel.h"
#include "correlmanager.h"
#include "irindex.h"
#include "fromto.h"




////////////////////////////////////////////////////
///	Class  : ARM_BSConvAdjust
///	Routine: Init Method
///	Returns: 
///	Action : Init
////////////////////////////////////////////////////

void ARM_BSConvAdjust::Init(void)
{
	SetName(ARM_BSCONVADJUST);

    itsUsedModel          = NULL;

    itsSUMMITFormulaeUsed = 0; 
    itsUseSabrCMS         = 0;
}



////////////////////////////////////////////////////
///	Class  : ARM_BSConvAdjust
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_BSConvAdjust::ARM_BSConvAdjust(int SUMMITFormulaeUsed, int UseSabrCMS)
{
	Init();

	itsSUMMITFormulaeUsed = SUMMITFormulaeUsed;

    itsUseSabrCMS         = UseSabrCMS;
}



////////////////////////////////////////////////////
///	Class  : ARM_BSConvAdjust
///	Routine: copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_BSConvAdjust::ARM_BSConvAdjust(const ARM_BSConvAdjust& rhs)
{
	Init();

	itsSUMMITFormulaeUsed = rhs.itsSUMMITFormulaeUsed;

    itsUseSabrCMS         = rhs.itsUseSabrCMS;

	itsUsedModel          = rhs.itsUsedModel;
}



////////////////////////////////////////////////////
///	Class  : ARM_BSConvAdjust
///	Routine: operator =
///	Returns: 
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////

ARM_BSConvAdjust& ARM_BSConvAdjust::operator= (const ARM_BSConvAdjust& rhs)
{
    (*this).ARM_Object::operator = (rhs);

    BitwiseCopy(&rhs);   

    return(*this);
}



////////////////////////////////////////////////////
///	Class  : ARM_BSConvAdjust
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_BSConvAdjust::~ARM_BSConvAdjust()
{
}



////////////////////////////////////////////////////
///	Class   : ARM_BSConvAdjust
///	Routines: BitwiseCopy,Copy,View,Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////

ARM_Object* ARM_BSConvAdjust::Clone(void)
{
    return new ARM_BSConvAdjust(*this);
}



void ARM_BSConvAdjust::View(char* id, FILE* ficOut)
{
	FILE* fOut;
    char  fOutName[200];
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
 
       (void) unlink(fOutName);
       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

    fprintf(fOut, "\n ======> BS Convexity Adjustment Manager\n");

	if ( ficOut == NULL )
    {
        fclose(fOut);
    }
}



                      /*****************************/
                      /*                           */
                      /*     LIBOR CASE            */
                      /*                           */
                      /*****************************/

// return the convexity adjustment of LIBOR paid in Arrears or in advance
// Correspond to TIMELAG formalae in the case of SUMMIT
// N.B; 
// Whatever is the flag : SUMREP,EXP,SUMANA,SUMEXP in this cases we call the
// FOLLOWING function 
double ARM_BSConvAdjust::ConvAdjustLibor(ARM_Model* Model, const AdjustLiborData& Input)
{
	if (!Model)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                             " ARM_BSConvAdjust::Model not available in ConvAdjustLibor adjustment");
    
    if (Input.Start-Input.Pay > PRECISION)
	{
		if( (Input.Reset - Input.Pay) > PRECISION )
		{
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
                             " payment date is less than start date, this case valid only if payment date >= reset date");
		}
	}

	try
    {
		double asof     = Model->GetStartDate().GetJulian();
		int dayCount    = Input.Daycount;
        
		double resetLag = (Input.Reset-asof)/K_YEAR_LEN;
		
        // MA+MW : because of the correction of resetDate in :
        // double ARM_BSModel::ExpectedFwdYield(resetDate calc PB)
        // resetLag could be < 0 SO!:

        if ( resetLag < 0 )
        {
           resetLag = 0.0;
        }

        double startLag = (Input.Start-asof)/K_YEAR_LEN;

		double rawFwd   = Input.RawFwd;
		int YieldDecompFreq = Input.YieldDecompFreq;
		double Margin = Input.Margin;
        double paymentlag;
		double VolFwd;

		if ( itsSUMMITFormulaeUsed != K_CONV_ADJ_LIBORORG ) // ( TimeLag() in Summit )
		{
			// Whatever is the flag : SUMREP,EXP,SUMANA,SUMEXP

			double Term_pe = (Input.End-Input.Pay)/K_YEAR_LEN;
			
			double Term_se = (Input.End - Input.Start)/K_YEAR_LEN;
			
			double Term_rp = (Input.Pay - Input.Reset)/K_YEAR_LEN;
			
			double Term_re = (Input.End - Input.Reset)/K_YEAR_LEN;

			int indexFreq = Input.IndexFreq;
			
			double tenor = (indexFreq>0)? (1.0/indexFreq) : Term_se;
			VolFwd   = Model->ComputeCvxAdjVol(resetLag, tenor, 0.0, 0.0)/100.0;	//volatility for forward that we will adjust

			/*!
			 * natural payement p/r payment date
			 */
			if ( fabs(Term_pe) < PRECISION )
			{
			   paymentlag = 1.0;
			}
			else
			{	
				/*!
				 * standard payement In Advance p/r payment date
				 */
				double RTsTp, RTsTe ;
				double Term_sp;

				if (Input.Start+PRECISION<Input.Pay)
				{
					Term_sp   = (Input.Pay - Input.Start)/K_YEAR_LEN; 
					if( Term_sp<=1.0 )
						RTsTp = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),ARM_Date(Input.Pay),
																	K_MM_RATE)/100.;//K_MM_RATE,K_ADJUSTED)/100.;
					else
					{
						RTsTp = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),ARM_Date(Input.Pay),YieldDecompFreq)/100.0;// taux actuariel//RTsTp = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),ARM_Date(Input.Pay),YieldDecompFreq,K_ADJUSTED)/100.0;// taux actuariel
					}

					RTsTe = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),ARM_Date(Input.End),
																	K_MM_RATE)/100.;//K_MM_RATE,K_ADJUSTED)/100.;
				}

				// We use the same adjustement formulae for any payment date between reset and start
				if ( (Input.Reset<=Input.Pay+PRECISION) && (Input.Pay<=Input.Start+PRECISION) )
				{
				   paymentlag= (100.0+Term_pe*rawFwd*exp(VolFwd*VolFwd*resetLag))/(100.0 + Term_pe*rawFwd);
				}
				else
				{
					// To have a payment lag equal to 1.0 when the payment date is close to the
					// end date.
					double Rho_sp    = 1.0;		// correlation between rawFwd and RTsTp
					double Rho_se    = 1.0;		// correlation between rawFwd and RTsTe
					ARM_CorrelManager* correls = Model->GetCorrelManager();
					
					if (correls)
					{						
						double indextype = (indexFreq>0) ? (1.0/indexFreq) : (1.0/ROUND(1.0/Term_se));
						string yearTerm(YearTermToStringMatu(indextype));
						ARM_IRIndex IRIndex((char*) (yearTerm.c_str()),Model->GetZeroCurve()->GetCurrencyUnit());  

						// mode "EURIB_EURIB" in Summit
						string ccy	= IRIndex.GetCurrencyUnit()->GetCcyName();
						string index;
						if( ccy == "EUR")
							index = "EURIBOR";
						else
							index = "LIBOR";

						string intraMktTag	= ccy < index ? ccy + "_" + index : index + "_"+ ccy;
						
						double TermRhosp = (indexFreq>0) ? (1.0/indexFreq) : Term_sp;
						double TermRhose = (indexFreq>0) ? (1.0/indexFreq) : Term_se;

						try
						{
							Rho_sp = correls->ComputeCorrelData( "IR/IR", intraMktTag, Term_rp, TermRhosp)
																	/ ARM_Constants::correlBase;
							Rho_se = correls->ComputeCorrelData( "IR/IR", intraMktTag, Term_re, TermRhose )
																	/ ARM_Constants::correlBase;
						}
						// temporaire with old key
						catch(...)
						{
							string indexName	= IRIndex.GetIndexName();
							string intraMktTagTmp = ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;
							Rho_sp = correls->ComputeCorrelData( "IR/IR", intraMktTagTmp, Term_rp, TermRhosp)
																	/ ARM_Constants::correlBase;
							Rho_se = correls->ComputeCorrelData( "IR/IR", intraMktTagTmp, Term_re, TermRhose )
																	/ ARM_Constants::correlBase;
						}
					}
        
					double RspVol = Model->ComputeCvxAdjVol(resetLag, Term_sp, 0.0, 0.0)/100.0;	// vol ATM de tx Rsp pour expiry resetLag
					double RseVol = Model->ComputeCvxAdjVol(resetLag, Term_se, 0.0, 0.0)/100.0;	// vol ATM de tx Rse pour expiry resetLag
					// ajustment for the RspVol when Term_sp>1.0
					if ( Term_sp > 1.0 )
					{
						ARM_Date oneYAfterStart = ARM_Date(Input.Start);
						oneYAfterStart.AddYears(1.0);

						double maxTerm = YieldDecompFreq - Term_sp>0.0 ? (YieldDecompFreq - Term_sp):0.0;
						double vol1y = Model->ComputeCvxAdjVol(resetLag, 1.0, 0.0, 0.0)/100.0;
						double tx1yLin = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),oneYAfterStart,K_MM_RATE,K_ADJUSTED)/100.0;// taux linear
						double tx1yAct = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),oneYAfterStart,1,K_ADJUSTED)/100.0;// taux actuariel
						RspVol = RspVol + ( pow(1+RTsTp/YieldDecompFreq,1-YieldDecompFreq)*tx1yLin/tx1yAct -1 )*maxTerm*vol1y;
					}

					double deltaTse = Term_se;
					double deltaTsp = Term_sp;
					
					if ( Term_pe < -PRECISION )
					{
						if ( Term_sp<=1.0 )
						{
							deltaTsp = Term_sp;
						}
						else
						   deltaTsp = 1/YieldDecompFreq;
					}

					double frac_se = (Term_se*RTsTe)/(1.0+deltaTse*RTsTe);
					double term1 = frac_se*Rho_se*VolFwd*RseVol*resetLag;
					double frac_sp = (Term_sp*RTsTp)/(1.0+deltaTsp*RTsTp);
					double term2 = frac_sp*Rho_sp*VolFwd*RspVol*resetLag;

					/*!
					 * payementIn Arreas p/r payment date 
					 */
					if (Term_pe < -PRECISION)
					{
						paymentlag = exp(term1-term2);
					}
             
					/*!
					 * payement In Advance p/r payment date 
					 */          
					else if (Term_pe > PRECISION)
					{
						double RTpTe = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Pay),ARM_Date(Input.End),
																						K_MM_RATE,K_ADJUSTED)/100.;    

						double frac_pe=(Term_pe*RTpTe)/(1.0+Term_pe*RTpTe); 
						double exp1=1.0/frac_pe* term1;
						double exp2=1.0/frac_pe* term2;
						paymentlag=(1.0/(1.0+Term_pe*RTpTe)+ frac_pe*exp(exp1-exp2));
					}
				}
			}
		}
		else // ---> original ARM native FORMULAE LIBORORG
		{
			double Term_pe = Input.End>Input.Pay ? CountYears(dayCount, Input.Pay, Input.End) :
											 -CountYears(dayCount, Input.End, Input.Pay);
			double Term_se  = CountYears(dayCount, Input.Start, Input.End);
			VolFwd   = Model->ComputeCvxAdjVol(resetLag, Term_se, 0.0, 0.0)/100.0;//volatility for forward that we will adjust

			/*!
			 * natural payement p/r payment date
			 */
			if ( fabs(Term_pe) < PRECISION )
			{
			   paymentlag = 1.0;
			}
			else
			{	
				/*!
				 * standard payement In Advance p/r payment date
				 */
				double RTsTp ;
				double Term_sp;

				if (Input.Start+PRECISION<Input.Pay)
				{
					Term_sp   = CountYears(dayCount, Input.Start, Input.Pay); //CountYears(dayCount, Input.Start, Input.Pay);
					RTsTp     = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),ARM_Date(Input.Pay),
																	K_MM_RATE,K_ADJUSTED)/100.;
				}
            
				if ( fabs(Input.Start-Input.Pay)<PRECISION )
				{
				   paymentlag= (100.0+Term_pe*rawFwd*exp(VolFwd*VolFwd*resetLag))/(100.0 + Term_pe*rawFwd);
				}
				else
				{
					// To have a payment lag equal to 1.0 when the payment date is close to the
					// end date.
					double Rho       = 1.0;		// correlation between rawFwd and RTsTp
					ARM_CorrelManager* correls = Model->GetCorrelManager();
					if (correls)
					{
						double indextype = 1.0/ROUND(1.0/Term_se);
						string yearTerm(YearTermToStringMatu(indextype));
						ARM_IRIndex IRIndex((char*) (yearTerm.c_str()),Model->GetZeroCurve()->GetCurrencyUnit());  

						// mode "CORR" in Summit
						string ccy = IRIndex.GetCurrencyUnit()->GetCcyName();

						string index = "CORR";
						
						string intraMktTag	= ccy < index ? ccy + "_" + index : index + "_"+ ccy;

						try
						{
							Rho = correls->ComputeCorrelData( "IR/IR", intraMktTag, startLag, Term_sp )
																/ ARM_Constants::correlBase;
						}
						// temporaire with old key
						catch(...)
						{
							string indexName	= IRIndex.GetIndexName();
							string intraMktTagTmp = ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;
							Rho = correls->ComputeCorrelData( "IR/IR", intraMktTagTmp, startLag, Term_sp )
																/ ARM_Constants::correlBase;
						}

					}
        
   					double RVol = Model->ComputeCvxAdjVol(resetLag, Term_sp, 0.0, 0.0)/100.0;	

					double frac_se = (Term_se*rawFwd)/(100.0+Term_se*rawFwd);
					double term1 = frac_se*VolFwd*VolFwd*resetLag;
					double frac_sp = (Term_sp*RTsTp)/(1.0+Term_sp*RTsTp);
					double term2 = frac_sp*Rho*VolFwd*RVol*resetLag;

					/*!
					 * payementIn Arreas p/r payment date 
					 */
					if ( Term_pe < -PRECISION )
					{
						paymentlag = exp(term1-term2);
					}
             
					/*!
					 * payement In Advance p/r payment date 
					 */           
					else if ( Term_pe > PRECISION )
					{
						double RTpTe = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Pay),ARM_Date(Input.End),
																						K_MM_RATE,K_ADJUSTED)/100.;    

						double frac_pe=(Term_pe*RTpTe)/(1.0+Term_pe*RTpTe); 
						double exp1=1.0/frac_pe* term1;
						double exp2=1.0/frac_pe* term2;
						paymentlag=(1.0/(1.0+Term_pe*RTpTe)+ frac_pe*exp(exp1-exp2));
					}
				}
			}
		}

		double fwdAdj = rawFwd*paymentlag;
		double xi = paymentlag;

		if ( YieldDecompFreq > 1 )
		{
			xi = 0.5*pow(VolFwd*fwdAdj,2.0);
			xi *=(1.0/YieldDecompFreq-1)*pow(1.0+fwdAdj+Margin,1.0/YieldDecompFreq-2.0)*resetLag;
			xi /= (YieldDecompFreq*(pow(1.0+fwdAdj+Margin,1.0/YieldDecompFreq-2.0)-1.0));
		}

		return xi;
    }

	catch (Exception& e)
	{
		throw e;
	}

	catch (...)
	{
		throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
                        "ARM_BSConvAdjust:: error in ConvAdjustLibor");
	}

    return(0.0);
}


                      /*****************************/
                      /*                           */
                      /*       END LIBOR CASE      */
                      /*                           */
                      /*****************************/






// compute the fwd volatility used in the libor adjustment
double ARM_BSConvAdjust::ComputeFwdVolatility(ARM_Model* Model, 
                                              const AdjustLiborData& Input)
{
	if (!Model)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                             " ARM_BSConvAdjust::Model not available in ConvAdjustLibor adjustment");
    
    if (Input.Start-Input.Pay > 2)
	{
		if(Input.Reset - Input.Pay > PRECISION )
		{
			throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
                             " payment date is less than start date, this case valid only if payment date >= reset date");
		}
	}
			
	try
    {
		double asof     = Model->GetStartDate().GetJulian();
		int dayCount    = Input.Daycount;		
		double resetLag = (Input.Reset-asof)/K_YEAR_LEN;

        // MA+MW

        if ( resetLag < 0 )
        {
           resetLag = 0;
        }

		double rawFwd   = Input.RawFwd;

		double Term_se  = CountYears(dayCount, Input.Start, Input.End);
		double VolFwd   = Model->ComputeCvxAdjVol(resetLag, Term_se, 0.0, 0.0)/100.0;// volatility for forward that we will adjust

		return VolFwd;
	}

	catch (...)
	{
		throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
                        "ARM_BSConvAdjust:: error in ConvAdjustLibor");
	}

    return(0.0);
}



//return the natural convexity adjustment of CMS 
double ARM_BSConvAdjust::NaturalAdjstCMS(ARM_Model* Model, const NaturalAdjustData& Input, StoreFwdRateInfo* StoreInfo)
{
	if (!Model)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                             " ARM_BSConvAdjust::Model not available in NaturalAdjstCMS adjustment");
	}
	int dayCount= Input.Daycount;
	double asof= Model->GetStartDate().GetJulian();

	double resetLag = (Input.Reset-asof)/K_YEAR_LEN;
	double Term_se  = (Input.End - Input.Start)/K_YEAR_LEN;
 	double theta	= 1.0/Input.Frequency;
	double N		= ROUND((Input.End-Input.Start)/K_YEAR_LEN)*Input.Frequency;
	double rawSwap	= Input.RawSwap/100.0;
    double swapVol	= Model->ComputeCvxAdjVol(resetLag, ROUND(Term_se), 0.0,0.0)/100.;
   	int YieldDecompFreq = Input.YieldDecompFreq;
	double Margin		= Input.Margin/100.0;
	double payLagAdj	= Input.PayLagAdj;


    double lambda=1.0-(theta*N*rawSwap)/((1.0+theta*rawSwap)*(pow((1.0+theta*rawSwap),N)-1.0));
    double exp1=exp(swapVol*swapVol*resetLag)-1.0;
    double ki = lambda*exp1;

	if (YieldDecompFreq > 1 || Margin != 0.0)
	{
		double convlag = (1+ki);
        
		if(itsSUMMITFormulaeUsed==K_CONV_ADJ_EXP)// 0
        {
            convlag*=payLagAdj;
        }

		double cmsadust = rawSwap * convlag;
		double swapRateDcpSpr = FromRateToRate((rawSwap+Margin)*100.0, 1.0,
                                      K_COMP_ANNUAL, YieldDecompFreq)/100.0;

		double swapRateDcpSprAdjust = FromRateToRate((rawSwap*convlag+Margin)*100.0, 1.0,
                                      K_COMP_ANNUAL, YieldDecompFreq)/100.0;

		double xi = 0.0;
		if (YieldDecompFreq > 1)
		{
			xi = 0.5*pow(swapVol*(cmsadust+Margin),2.0);
			xi *=(1.0/YieldDecompFreq-1)*pow(1.0+cmsadust+Margin,1.0/YieldDecompFreq-2.0)*resetLag;
			xi /= (YieldDecompFreq*(pow(1.0+cmsadust+Margin,1.0/YieldDecompFreq)-1.0));
		}

		ki = swapRateDcpSprAdjust/swapRateDcpSpr*exp(xi)/payLagAdj-1.0;
	}

	return ki+1.0;

}



// compute the cms volatility used for the payment lag analytic formulae
double ARM_BSConvAdjust::ComputeCMSVolatility(ARM_Model* Model, const NaturalAdjustData& Input, StoreFwdRateInfo* StoreInfo)
{
	if (!Model)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                             " ARM_BSConvAdjust::Model not available in NaturalAdjstCMS adjustment");
	}
	int dayCount= Input.Daycount;
	double asof= Model->GetStartDate().GetJulian();

	double resetLag = (Input.Reset-asof)/K_YEAR_LEN;
	double Term_se  = (Input.End - Input.Start)/K_YEAR_LEN;
	double rawSwap	= Input.RawSwap/100.0;
    double swapVol	= Model->ComputeCvxAdjVol(resetLag, ROUND(Term_se), 0.0,0.0)/100.;
	double fwdAdj	= Input.Fwd;

	double cmsVol	= swapVol;

	if (Input.OptionType == K_CAP || Input.OptionType == K_FLOOR) // Compute the value and the BS vol implicit of an oplet
	{
		double espOplet = Model->EuroCaplet(0.,
                               (Input.Reset-asof)/K_YEAR_LEN,
                               (Input.Start-asof)/K_YEAR_LEN, 
                               (Input.End-asof)/K_YEAR_LEN,
                               (Input.Pay-asof)/K_YEAR_LEN,
                               Input.OptionType, 
                               fwdAdj, 
                               Input.Strike, 
                               Input.DomOrFrg,
                               Input.CompMeth, 
                               Input.Daycount,
                               true,
                               false,
                               Input.Frequency,
                               ROUND(Term_se),
                               Input.YieldDecompFreq,
                               Input.Margin,
							   Model->GetZeroCurve()->GetCurrencyUnit());
/*
		double espOplet = EuroCapletWithoutTimelag(
			Model,
			asof, 
			(Input.Reset-asof)/K_YEAR_LEN,
			(Input.Start-asof)/K_YEAR_LEN, 
			(Input.End-asof)/K_YEAR_LEN,
			(Input.Pay-asof)/K_YEAR_LEN,
			Input.OptionType, 
			Input.RawSwap,
			fwdAdj, 
			Input.Strike, 
			Input.DomOrFrg,
			Input.CompMeth, 
			Input.Daycount,
			true, 
			N,
			Input.Frequency,
			Input.YieldDecompFreq, 
			Input.Margin);
*/

		 // To prevent of bsVolImp with unreachable price

		double minPrice = 0., maxPrice = 0.;
		if (Input.OptionType == K_CALL)
		{
			minPrice = fwdAdj-Input.Strike;
			if (minPrice > 0)
				minPrice = 0.0;   
			maxPrice = fwdAdj;
		}
		else if (Input.OptionType == K_PUT)
		{
			minPrice = Input.Strike-fwdAdj;
			if (minPrice > 0)
				minPrice = 0.0;
			maxPrice = Input.Strike;
		}

		double volCMSCaplet = 0.0;

		if ((minPrice < espOplet) && (espOplet < maxPrice))
		{
			volCMSCaplet = bsVolImp(espOplet, fwdAdj, Input.Strike, 0, 0, (Input.Reset-asof)/K_YEAR_LEN, Input.OptionType, 1e-8, swapVol);
		}
		else
		{
			volCMSCaplet = -100.0;
		}

		cmsVol = volCMSCaplet;
	}

	return cmsVol;
}



//return the convexity adjustment of CMS paid in Arrears or in advance
double ARM_BSConvAdjust::TimeLagAdjustCMS(ARM_Model* Model, const TimeLagAdjustData& Input, StoreFwdRateInfo* StoreInfo)
{
	if (!Model)
	{
	   throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
             " ARM_BSConvAdjust::Model not available in TimeLagAdjustCMS adjustment");
	}

    double paymentlagAdjust = 1.0;

    int dayCount= Input.Daycount;
    double asof= Model->GetStartDate().GetJulian();
    double resetLag =  (Input.Reset-asof)/K_YEAR_LEN;

    double theta   = (Input.Pay - Input.Reset)/360;

    double Term_se   = (Input.End-Input.Start)/K_YEAR_LEN;
    double rawSwap	= Input.RawSwap/100.0;
    double YieldDecomp = Input.YieldDecompFreq;
    double Margin = Input.Margin/100.0;
    
    // Time lag formulae implemented in SUMMIT
	if (itsSUMMITFormulaeUsed==K_CONV_ADJ_SUMREP)//==1
    {   
        double Term_rp   = (Input.Pay - Input.Reset)/K_YEAR_LEN;

	    if (Input.Pay - Input.Start != 0) //(Term_rp > 0)
	    {
		    double	rho = CorrelDataDefault::defaultCorrel;

		    ARM_CorrelManager* correls = Model->GetCorrelManager();
		    if (correls)
		    {
				// mode "EURIB_EURIB" in Summit
				string ccy	= Model->GetZeroCurve()->GetCurrencyUnit()->GetCcyName();

				string index = "CORR";
						
				string intraMktTag	= ccy < index ? ccy + "_" + index : index + "_"+ ccy;

				try
				{
					rho = Model->GetCorrelManager()->ComputeCorrelData( "IR/IR", intraMktTag, Term_rp, floor(Term_se))
												    / ARM_Constants::correlBase;
				}
				// temporaire with old key
				catch(...)
				{
					int NbfixFlows= ROUND(Term_se);
					char iindexName[20]; 
					sprintf(iindexName,"CMS%d",NbfixFlows);
					string indexName	= iindexName;
					string intraMktTagTmp = ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;
					rho = Model->GetCorrelManager()->ComputeCorrelData( "IR/IR", intraMktTagTmp, Term_rp, floor(Term_se))
												    / ARM_Constants::correlBase;
				}
		    }

	        double R1TrTp    = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Reset),ARM_Date(Input.Pay),K_MM_RATE,K_ADJUSTED)/100.0;

		    NaturalAdjustData naturalAdjutData;

		    naturalAdjutData.RawSwap			= Input.RawSwap;
		    naturalAdjutData.Reset				= Input.Reset;
		    naturalAdjutData.Start				= Input.Start;
		    naturalAdjutData.End				= Input.End;
		    naturalAdjutData.Frequency			= Input.Frequency;
		    naturalAdjutData.Daycount			= Input.Daycount;
		    naturalAdjutData.YieldDecompFreq	= K_COMP_PROP;
		    naturalAdjutData.Margin				= 0.0;

			naturalAdjutData.PayLagAdj	= 1.0;

			naturalAdjutData.OptionType			= 0;
			naturalAdjutData.Strike				= 0.0;

		    double CMSVol    = Model->GetConvAdjustManager()->ComputeCMSVolatility(Model, naturalAdjutData);
            double adjustNoDecap = Model->GetConvAdjustManager()->NaturalAdjstCMS(Model,naturalAdjutData);

            double fwdAdjNoDecap = rawSwap*adjustNoDecap;

            naturalAdjutData.YieldDecompFreq=Input.YieldDecompFreq;
            naturalAdjutData.Margin=Input.Margin;

            double adjustDecap = Model->GetConvAdjustManager()->NaturalAdjstCMS(Model,naturalAdjutData);

            // double fwdAdjDecap = FromRateToRate((rawSwap+Margin)*100.0, 1.0,K_COMP_ANNUAL, Input.YieldDecompFreq)*adjustDecap/100.0;
			double fwdAdjDecap = rawSwap*adjustDecap;

		    double RVol = 0.0;

		    if (resetLag > 0)
		    {
                double xi = 0.0;
                double volAdj = 1.0;

                if (YieldDecomp > 1)
		        {                                      
                    xi = 0.5*pow(CMSVol*fwdAdjNoDecap,2.0);
			        xi *=(1.0/YieldDecomp-1)*pow(1.0+fwdAdjNoDecap,1.0/YieldDecomp-2.0)*resetLag;
			        xi /= (YieldDecomp*(pow(1.0+fwdAdjNoDecap,1.0/YieldDecomp)-1.0));

                    volAdj = pow(1.0+fwdAdjNoDecap,1.0/YieldDecomp-1.0);
			        volAdj *= fwdAdjNoDecap;
			        volAdj /= FromRateToRate(fwdAdjNoDecap*100.0, 1.0,
											K_COMP_ANNUAL, YieldDecomp)/100.0;
			    }

                CMSVol *= volAdj;

			    // Bad code to manage caplet & swaption volatility with extended BS model.
			    ARM_BSModel* BSModel = dynamic_cast<ARM_BSModel*>(Model);

			    if (BSModel)
			    {
				   RVol = BSModel->ComputeCvxAdjVol(resetLag,Term_rp,0.0,0.0)/100.;
			    }
			    else
                    throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                                 " Model::ComputeCvxAdjVol not available in TimeLagAdjustCMS adjustment");

			    double varR = R1TrTp*(1+R1TrTp*theta/(1+theta*R1TrTp)*(exp(RVol*RVol*resetLag)-1))*sqrt(exp(RVol*RVol*resetLag)-1);

                double fwdAdjDecapSynth = FromRateToRate((fwdAdjNoDecap)*100.0, 1.0,K_COMP_ANNUAL, Input.YieldDecompFreq)*exp(xi)/100.0;

                double varCMS = fwdAdjDecapSynth*sqrt(exp(CMSVol*CMSVol*resetLag)-1);
                // varCMS /= fwdAdjDecap;

                double cov = rho*varR*varCMS;

                paymentlagAdjust=1-theta/(1+theta*R1TrTp)*cov/fwdAdjDecap;
		    }
	    }
    }
	else if(itsSUMMITFormulaeUsed==K_CONV_ADJ_EXP)//==0
    {
        double Term_sp   = (Input.Pay - Input.Start)/K_YEAR_LEN;

	    if (Term_sp > 0)
	    {
		    double	rho = CorrelDataDefault::defaultCorrel;

		    ARM_CorrelManager* correls = Model->GetCorrelManager();
		    if (correls)
		    {
				// mode "EURIB_EURIB" in Summit
				string ccy	= Model->GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
				string index = "CORR";				
				string intraMktTag	= ccy < index ? ccy + "_" + index : index + "_"+ ccy;

				try
				{
					rho = Model->GetCorrelManager()->ComputeCorrelData( "IR/IR", intraMktTag, Term_sp, floor(Term_se))
												    / ARM_Constants::correlBase;
				}
				// temporaire with old key
				catch(...)
				{
					int NbfixFlows= ROUND(Term_se);
					char iindexName[20]; 
					sprintf(iindexName,"CMS%d",NbfixFlows);
					string indexName	= iindexName;
					string intraMktTagTmp	= ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;
					rho = Model->GetCorrelManager()->ComputeCorrelData( "IR/IR", intraMktTagTmp, Term_sp, floor(Term_se))
														/ ARM_Constants::correlBase;
				}
		    }

	        double R1TsTp    = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),ARM_Date(Input.Pay),K_MM_RATE,K_ADJUSTED)/100.;

		    NaturalAdjustData naturalAdjutData;

		    naturalAdjutData.RawSwap			= Input.RawSwap;
		    naturalAdjutData.Reset				= Input.Reset;
		    naturalAdjutData.Start				= Input.Start;
		    naturalAdjutData.End				= Input.End;
		    naturalAdjutData.Frequency			= Input.Frequency;
		    naturalAdjutData.Daycount			= Input.Daycount;
		    naturalAdjutData.YieldDecompFreq	= Input.YieldDecompFreq;
		    naturalAdjutData.Margin				= Input.Margin;

			naturalAdjutData.OptionType			= 0;
			naturalAdjutData.Strike				= 0.0;

		    double CMSVol    = Model->GetConvAdjustManager()->ComputeCMSVolatility(Model, naturalAdjutData, StoreInfo);

		    double RVol = 0.0;

		    if ( resetLag > 0 )
		    {
			    // Bad code to manage caplet & swaption volatility with extended BS model.
			    ARM_BSModel* BSModel = dynamic_cast<ARM_BSModel*>(Model);

			    if(BSModel)
			    {
				    RVol = BSModel->ComputeCvxAdjVol(resetLag,Term_sp,0.0,0.0)/100.;
			    }
			    else
                    throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                                 " Model::ComputeCvxAdjVol not available in TimeLagAdjustCMS adjustment");

			    double nominator=Term_sp*R1TsTp*rho*CMSVol*RVol*resetLag;
			    double donominator=1.0+Term_sp*R1TsTp;
			    paymentlagAdjust=exp(-nominator/donominator);
		    }
        }
    }
	else if(itsSUMMITFormulaeUsed==K_CONV_ADJ_SUMANA)//CMSE approch
	{
        double Term_rp   = (Input.Pay - Input.Reset)/K_YEAR_LEN;
	    if (Term_rp > 0)
	    {
		    double	rho = CorrelDataDefault::defaultCorrel;

		    ARM_CorrelManager* correls = Model->GetCorrelManager();
		    if (correls)
		    {
				// mode "EURIB_EURIB" in Summit
				string ccy	= Model->GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
				string index = "CORR";				
				string intraMktTag	= ccy < index ? ccy + "_" + index : index + "_"+ ccy;
				try
				{
					rho = Model->GetCorrelManager()->ComputeCorrelData( "IR/IR", intraMktTag, Term_rp, floor(Term_se))
												    / ARM_Constants::correlBase;
				}
				// temporaire with old key
				catch(...)
				{
					int NbfixFlows= ROUND(Term_se);
					char iindexName[20]; 
					sprintf(iindexName,"CMS%d",NbfixFlows);
					string indexName	= iindexName;
					string intraMktTagTmp	= ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;
					rho = Model->GetCorrelManager()->ComputeCorrelData( "IR/IR", intraMktTagTmp, Term_rp, floor(Term_se))
														/ ARM_Constants::correlBase;					
				}
		    }

	        double R1TrTp   = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Reset),ARM_Date(Input.Pay),K_MM_RATE,K_ADJUSTED)/100.0;

		    NaturalAdjustData naturalAdjutData;

		    naturalAdjutData.RawSwap			= Input.RawSwap;
		    naturalAdjutData.Reset				= Input.Reset;
		    naturalAdjutData.Start				= Input.Start;
		    naturalAdjutData.End				= Input.End;
		    naturalAdjutData.Frequency			= Input.Frequency;
		    naturalAdjutData.Daycount			= Input.Daycount;
		    naturalAdjutData.YieldDecompFreq	= K_COMP_PROP;
		    naturalAdjutData.Margin				= 0.0;

			naturalAdjutData.PayLagAdj	= 1.0;

			naturalAdjutData.OptionType			= 0;
			naturalAdjutData.Strike				= 0.0;

		    double CMSVol    = Model->GetConvAdjustManager()->ComputeCMSVolatility(Model, naturalAdjutData);

            naturalAdjutData.YieldDecompFreq=Input.YieldDecompFreq;
            naturalAdjutData.Margin=Input.Margin;

		    double RVol = 0.0;

		    if (resetLag > 0)
		    {
                double xi = 0.0;
                double volAdj = 1.0;

			    ARM_BSModel* BSModel = dynamic_cast<ARM_BSModel*>(Model);

			    if (BSModel)
			    {
				    RVol = BSModel->ComputeCvxAdjVol(resetLag,Term_rp,0.0,0.0)/100.;
			    }
			    else
                    throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                                 " Model::ComputeCvxAdjVol not available in TimeLagAdjustCMS adjustment");
				//theta = Term_rp;
				paymentlagAdjust = 1 - rho*sqrt(exp(pow(theta*RVol*R1TrTp*sqrt(resetLag)/(1+theta*R1TrTp),2))-1)*sqrt(exp(CMSVol*CMSVol*resetLag)-1);

		    }
	    }		
	}
	else if( itsSUMMITFormulaeUsed == K_CONV_ADJ_SUMEXP)  // CMSA approach
	{
		double Term_sp = (Input.Pay - Input.Start)/K_YEAR_LEN; // en base ACT/365
		double Term_rp = (Input.Pay - Input.Reset)/K_YEAR_LEN; // en base ACT/365
		int Freq = Input.Frequency;
		
		if ( Term_sp > 0.0)
		{
			//---- Correlation entre le CMS et Rsp pour expiry=resetLag ----// à voir
			double rho = CorrelDataDefault::defaultCorrel;
		    ARM_CorrelManager* correls = Model->GetCorrelManager();
		    if (correls)
		    {
				// mode "EURIB_EURIB" in Summit
				string ccy	= Model->GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
				/*
				string index = "CORR";				
				*/

				string index;
				if( ccy == "EUR")
					index = "EURIBOR";
				else
					index = "LIBOR";

				string intraMktTag	= ccy < index ? ccy + "_" + index : index + "_"+ ccy;
				try
				{
					rho = Model->GetCorrelManager()->ComputeCorrelData( "IR/IR", intraMktTag, Term_sp, floor(Term_se))/ ARM_Constants::correlBase;
				}
				// temporaire with old key
				catch(...)
				{
					int NbfixFlows= ROUND(Term_se);
					char iindexName[20]; 
					sprintf(iindexName,"CMS%d",NbfixFlows);
					string indexName	= iindexName;
					string intraMktTagTmp	= ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;
					rho = Model->GetCorrelManager()->ComputeCorrelData( "IR/IR", intraMktTagTmp, Term_sp, floor(Term_se))/ ARM_Constants::correlBase;				
				}
		    }
			
			//---- Vol ATM du CMS pour expiry=resetLag ----//
			NaturalAdjustData naturalAdjutData;
		    naturalAdjutData.RawSwap			= Input.RawSwap;
		    naturalAdjutData.Reset				= Input.Reset;
		    naturalAdjutData.Start				= Input.Start;
		    naturalAdjutData.End				= Input.End;
		    naturalAdjutData.Frequency			= Input.Frequency;
		    naturalAdjutData.Daycount			= Input.Daycount;
		    naturalAdjutData.YieldDecompFreq	= Input.YieldDecompFreq;
		    naturalAdjutData.Margin				= Input.Margin;

			naturalAdjutData.OptionType			= 0;
			naturalAdjutData.Strike				= 0.0;

			double CMSVol = Model->GetConvAdjustManager()->ComputeCMSVolatility(Model, naturalAdjutData);
			
			double RspVol = 0.0;
			double RTsTp = 0.0;	

			if( resetLag > 0.0 )
			{
				//---- Vol ATM du taux Rsp pour expiry=resetLag ----//
			    ARM_BSModel* BSModel = dynamic_cast<ARM_BSModel*>(Model);

				if(BSModel)
				{
					RspVol = BSModel->ComputeCvxAdjVol(resetLag,Term_sp,0.0,0.0)/100.;
					
					// ajustment for the RspVol when Term_sp>1.0
					if( Term_sp > 1.0 )
					{
						ARM_Date oneYAfterStart = ARM_Date(Input.Start);
						oneYAfterStart.AddYears(1.0);

						double maxTerm = Freq - Term_sp>0.0 ? (Freq - Term_sp):0.0;
						double vol1y = Model->ComputeCvxAdjVol(resetLag, 1.0, 0.0, 0.0)/100.0;
						double tx1yLin = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),oneYAfterStart,K_MM_RATE,K_ADJUSTED)/100.0;// taux linear
						double tx1yAct = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),oneYAfterStart,1,K_ADJUSTED)/100.0;// taux actuariel
						RspVol = RspVol + ( pow(1+RTsTp/Freq,1-Freq)*tx1yLin/tx1yAct -1 )*maxTerm*vol1y;
					}
				}
				else
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
							" Model::ComputeCvxAdjVol not available in TimeLagAdjustCMS adjustment");
				}
				
				//---- Taux de placement sur la période [Ts, Tp] ----//
				if( Term_sp <= 1.0 )
				{
					// taux monétaire
					RTsTp = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),ARM_Date(Input.Pay),K_MM_RATE,K_ADJUSTED)/100.0;
					paymentlagAdjust = exp(-Term_sp*RTsTp*rho*CMSVol*RspVol*resetLag/(1.0+RTsTp*Term_sp));
				}
				else 
				{
					// taux actuariel
					RTsTp = Model->GetZeroCurve()->ForwardYield(ARM_Date(Input.Start),ARM_Date(Input.Pay),Freq,K_ADJUSTED)/100.0;
					paymentlagAdjust = exp(-Term_sp*RTsTp*rho*CMSVol*RspVol*resetLag/(1.0+RTsTp/Freq)); 
				}
			}
		}
	}
	else
	{
	   throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
             " ARM_BSConvAdjust::Formulae unkown in TimeLagAdjustCMS adjustment. They should be : EXP, SUMEXP, SUMREP, SUMANA");

	}

    return paymentlagAdjust;
}


// return the forward CPI ratio convexity adjustment
double ARM_BSConvAdjust::FwdCPIRatioAdjust(ARM_Model* Model, const ARM_Vector* Input)
{
	double TjForMaturity	= (*Input)[0]; // CountYear(dayCount,modelAsOfDate,Tj)
	double TiForMaturity	= (*Input)[1]; // CountYear(dayCount,modelAsOfDate,Ti)
	double tenor			= (*Input)[2]; // CountYear(dayCount,Tj,Ti)
	double volTj			= (*Input)[3]; // InfZCVol->ComputeVolatility(SpotDate,Tj)
	double volTi			= (*Input)[4]; // InfZCVol->ComputeVolatility(SpotDate,Ti)
	double volYtY			= (*Input)[5]; // InfZCVol->ComputeVolatility(Tj,Ti)
	double forward			= (*Input)[6]; // ZC->ForwardYield(Tj,Ti)
	double volIR			= (*Input)[7];// ZCVol->ComputeVolatility(Tj,Tenor)
	double RhoInfIR			= (*Input)[8]; //InfIRCorr->ComputeCorrelDate(Tj,Ti)

	double tenorPayment		= (*Input)[9];
	double forwardPayment	= (*Input)[10];
	double volIRPayment		= (*Input)[11];
	double RhoInfIRPaymentj	= (*Input)[12];
	double RhoInfIRPaymenti	= (*Input)[13];

	double volBond = tenor * forward * volIR / ( 1.0 + tenor * forward );
	double volBondPayment = tenorPayment * forwardPayment * volIRPayment / ( 1.0 + tenorPayment * forwardPayment );

	/// computation of the adjustment
	/// part on the inflation only
	double logAdjInf	=  0.5 * ( volYtY * volYtY * tenor + volTj * volTj * TjForMaturity - volTi * volTi * TiForMaturity );
	/// part on the inflation vs the interest rates
	double logAdjInfbond=  TjForMaturity * volTj * RhoInfIR * volBond
		//// payment lag!
		+ TjForMaturity * volTj * RhoInfIRPaymentj * volBondPayment
		- TiForMaturity * volTi * RhoInfIRPaymenti * volBondPayment;

	double adj = exp(logAdjInfbond+logAdjInf);

	return adj;
}



/*----------------------------------------------------------------------------*/
/*---- End of file ----*/