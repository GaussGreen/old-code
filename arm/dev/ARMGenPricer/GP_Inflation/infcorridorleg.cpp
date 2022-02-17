
/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	 /
 *																							 *
 *			Class Inf Corridor Leg															 *
 *																							 *
 *			This class builds a corridor leg from swap legs	and inherits from Inf Swap Leg	 *									 *
 *																							 *
 *			Author	: Mathieu Bernardo														 *
 *			Date	: October, 25th 2006													 *																											 *
 *			Copyright (c) NatIxis															 *
 *			Version	: 1.0																	 *
 *																							 *
 *	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
#pragma warning(disable : 4786)

#include <gpbase\assignop.h>
#include "gpinflation/infcorridorleg.h"
#include "gpinflation/infleg.h"			
#include "gpinflation/infdata.h"		
#include "gpinflation/infcurvmodel.h"
#include "gpinflation/infbsmodel.h"
#include <gpbase/typedef.h>

#include "gpbase/stringconvert.h"
#include "gpbase/gpvector.h"
#include <gpbase/curve.h>
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/curvetypedef.h"
#include "gpcalib\densityfunctors.h"


#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/normal.h"

CC_BEGIN_NAMESPACE( ARM )


ARM_InfCorridorLeg::ARM_InfCorridorLeg(	map<string,ARM_Date>	&	mDate, 
										map<string,string>		&	mString,
										map<string,ARM_Curve*>	&	mCurve,
										map<string,int>			&	mInt,
										map<string,double>		&	mDouble){

	ARM_IRIndex *	tmpIrIndex		= NULL;
	char *			tmpInfResetCal	= "INF";
	double			tmpSpread		= 0.0;
	double			tmpMultiple		= 1.0;
	double			tmpComultiple	= 1.0;
	ARM_Currency*	tmpCcy			= new ARM_Currency	( mString["currency"].c_str() );
	ARM_INDEX_TYPE  tmpIndex		= (ARM_INDEX_TYPE)	mInt["irIndex"];
	char *			tmpResetCal		= const_cast< char *> (mString["resetCal"].c_str() );
	char *			tmpPayCal		= const_cast< char *> (mString["payCal"].c_str() );
	char *			tmpRefDate		= const_cast< char *> (mString["refDate"].c_str() );;

	itsEpsilon						= mDouble[	"epsilon"			];
	itsNbGaussLeg					= mDouble[	"nbGaussLeg"		];
	itsPtGaussLeg					= mDouble[	"ptGaussLeg"		];

	itsAsOfDate						= mDate	[	"asOfDate"			];
	itsStartDate					= mDate	[	"startDate"			];
	itsEndDate						= mDate	[	"endDate"			];
	itsRecOrPay						= mInt	[	"recOrPay"			];
	isIrCritera						= mInt	[	"irCritera"			];
	isModulable						= mInt	[	"isModulable"		];
	isComputed						= false;

	tmpIrIndex	= new ARM_IRIndex(	mInt	[	"cpnDayCount"		], 
									mInt	[	"resetFreq"			], 
									mInt	[	"payFreq"			], 
									mInt	[	"maturity"			], 
									mInt	[	"compMeth"			], 
									mInt	[	"fwdRule"			], 
									mInt	[	"irIntRule"			], 
									mInt	[	"resetTiming"		],  
									mInt	[	"resetGap"			],
									mInt	[	"payTiming"			],
									mInt	[	"payGap"			], 
									tmpCcy,
									tmpIndex,
									mInt	[	"decompFreq"		]);


	itsIrLeg	= ARM_CountedPtr<ARM_SwapLeg> (new ARM_SwapLeg(	
									mDate	[	"startDate"			],			
									mDate	[	"endDate"			],	
									tmpIrIndex,
									mInt	[	"recOrPay"			],
									tmpSpread ,
									mInt	[	"stubRule"			],
									mInt	[	"decompFreq"		],  
									tmpCcy,
									mInt	[	"irDayCount"		],
									10000,					//	to be in agreement with the schedule of Cap Floor Pricer // resetGap, 
									tmpResetCal,
									tmpPayCal,
									mInt	[	"decompPriceFlag"	],
									mInt	[	"finNotioType"		],
									tmpRefDate,
									mInt	[	"adjFirstRule"		]) );

	itsInfLeg	= ARM_CountedPtr<ARM_InfLeg>(  new ARM_InfLeg(	
									mDate	[	"startDate"			],
									mDate	[	"endDate"			],
									mString [	"infIndex"			],
									K_YEARTOYEAR_LEG,
									mInt	[	"recOrPay"			],
									mInt	[	"infInterType"		],
									tmpMultiple,
									tmpComultiple,
									tmpSpread, 
									mInt	[	"resetFreq"			],
									mInt	[	"infDayCount"		],
									tmpInfResetCal,
									mInt	[	"fwdRule"			],
									mInt	[	"infIntRule"		],
									mInt	[	"stubRule"			],
									mInt	[	"resetNumGap"		],
									mInt	[	"resetDemGap"		],
									mInt	[	"payFreq"			],
									mInt	[	"payGap"			],	
									tmpPayCal,
									mInt	[	"adjFirstRule"		],
									mInt	[	"finNotioType"		],
									mInt	[	"firstReset"		],
									tmpCcy) );

	BuildSchedule( );
	InitSchedule (	mCurve[	"notional"		],	
					mCurve[	"irLeverage"	],
					mCurve[	"infLeverage"	],
					mCurve[	"constant"		],
					mCurve[	"multipleUp"	],	
					mCurve[	"multipleDown"	],	
					mCurve[	"rangeUp"		],	
					mCurve[	"rangeDown"		] );

	if ( tmpIrIndex){ delete tmpIrIndex;	tmpIrIndex=NULL; }
	if ( tmpCcy   ) { delete tmpCcy;		tmpCcy=NULL;	 }

}


void ARM_InfCorridorLeg::BuildSchedule(){

	itsNbFlows = itsIrLeg->GetFlowStartDates()->GetSize();

	itsFlowStartDates	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );	
	itsFlowEndDates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );	;		
	itsPaymentDates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );			 
	itsInterestDays		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		
	itsInterestTerms	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		
	itsDiscFactor		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );			

	itsFwdRateStartDates= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );	
	itsFwdRateEndDates	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );	
	itsIrResetDates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		 
	itsIrResetTerms		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		
	itsIrInterTerms		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsFwdIrRates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );			
	itsIrAdjConv		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );				

	itsNumResetDates	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		
	itsDemResetDates	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		
	itsInfResetDates	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		 
	itsInfResetTerms	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		
	itsInfInterTerms	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsNumCPIRates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		
	itsDemCPIRates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsNumPublishDates	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsDemPublishDates	= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsFwdInfRates		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );			
	itsInfAdjConv		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );				

	itsIrLeverage		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );
	itsInfLeverage		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );	
	itsConstant			= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );			
	itsUpMultiple		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );			
	itsDownMultiple		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );		
	itsRangeUp			= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );				
	itsRangeDown		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );				
	itsNotional			= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );				
	
	itsUpPrice			= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );	
	itsDownPrice		= ARM_GP_VectorPtr(new ARM_GP_Vector( itsNbFlows ) );	

	for ( int i =0; i< itsNbFlows; i++){

		itsFlowStartDates	->Elt(i) =	itsIrLeg	->GetFlowStartDates()	->Elt(i);
		itsFlowEndDates		->Elt(i) =  itsIrLeg	->GetFlowEndDates()		->Elt(i);
		itsPaymentDates		->Elt(i) =  itsIrLeg	->GetPaymentDates()		->Elt(i);
		itsFwdRateStartDates->Elt(i) =  itsIrLeg	->GetFwdRateStartDates()->Elt(i);
		itsFwdRateEndDates	->Elt(i) =  itsIrLeg	->GetFwdRateEndDates()	->Elt(i);
		itsIrResetDates		->Elt(i) =  itsIrLeg	->GetResetDates()		->Elt(i);
		itsInfResetDates	->Elt(i) =  itsInfLeg	->GetResetDates()		->Elt(i);
		itsNumResetDates	->Elt(i) =  itsInfLeg	->GetNumResetDates()	->Elt(i);
		itsDemResetDates	->Elt(i) =	itsInfLeg	->GetDenomResetDates()	->Elt(i);
		
	}
}


void ARM_InfCorridorLeg::InitSchedule	(	ARM_Curve	*	notional,
											ARM_Curve	*	irLeverage,
											ARM_Curve	*	infLeverage,
											ARM_Curve	*	constant,						
											ARM_Curve	*	multipleUp,						
											ARM_Curve	*	multipleDown,						
											ARM_Curve	*	rangeUp,						
											ARM_Curve	*	rangeDown){
	double	lag;
	double	tmpIr;
	double	tmpInf;
	double	tmp;

	for(int i =0; i<itsNbFlows; i++){
		lag	= itsIrLeg ->GetResetDates() ->Elt(i) - itsAsOfDate.GetJulian();
		itsNotional		->Elt(i)	= notional		->Interpolate(lag);
		itsIrLeverage	->Elt(i)	= irLeverage	->Interpolate(lag);
		itsInfLeverage	->Elt(i)	= infLeverage	->Interpolate(lag);
		itsConstant		->Elt(i)	= constant		->Interpolate(lag);  
		itsRangeUp		->Elt(i)	= rangeUp		->Interpolate(lag);   
		itsRangeDown	->Elt(i)	= rangeDown		->Interpolate(lag);
		itsUpMultiple	->Elt(i)	= multipleUp	->Interpolate(lag);
		itsDownMultiple	->Elt(i)	= multipleDown	->Interpolate(lag); 

		itsIrInterTerms	->Elt(i)	= itsIrLeg -> GetInterestTerms() ->Elt(i);
		itsInfInterTerms->Elt(i)	= itsInfLeg-> GetInterestTerms() ->Elt(i);

		itsInterestTerms->Elt(i)	= itsIrLeg	-> GetInterestTerms()->Elt(i);
		itsInterestDays	->Elt(i)	= itsIrLeg	-> GetInterestDays() ->Elt(i);

		itsConstant		->Elt(i)	= itsConstant->Elt(i) - RateBase*itsInfLeverage->Elt(i);

		if( isModulable){
			tmpIr	=	itsIrInterTerms	->Elt(i);
			tmpInf	=	itsInfInterTerms->Elt(i);
			if (isIrCritera){
				tmp = tmpInf/tmpIr;
				itsUpMultiple	->Elt(i) = itsUpMultiple	->Elt(i) * tmp;
				itsDownMultiple	->Elt(i) = itsDownMultiple	->Elt(i) * tmp;
				itsRangeUp		->Elt(i) = itsRangeUp		->Elt(i) * tmp;;
				itsRangeDown	->Elt(i) = itsRangeDown		->Elt(i) * tmp;
			}
			else{
				tmp = tmpIr/tmpInf;
				itsUpMultiple	->Elt(i) = itsUpMultiple	->Elt(i) * tmp;
				itsDownMultiple	->Elt(i) = itsDownMultiple	->Elt(i) * tmp;
				itsRangeUp		->Elt(i) = itsRangeUp		->Elt(i) * tmp;
				itsRangeDown	->Elt(i) = itsRangeDown		->Elt(i) * tmp;
			}
		}
	}
}

ARM_InfCorridorLeg::ARM_InfCorridorLeg( const ARM_InfCorridorLeg& 	leg	){


	itsInfLeg			= CreateClonedPtr ( &*leg.GetInfLeg()			); 
	itsIrLeg			= CreateClonedPtr ( &*leg.GetIrLeg()			);

	itsFlowStartDates	= CreateClonedPtr ( &*	leg.itsFlowStartDates	);
	itsFlowEndDates		= CreateClonedPtr ( &*	leg.itsFlowEndDates		);
	itsPaymentDates		= CreateClonedPtr ( &*	leg.itsPaymentDates		);
	itsInterestDays		= CreateClonedPtr ( &*	leg.itsInterestDays		);
	itsInterestTerms	= CreateClonedPtr ( &*	leg.itsInterestTerms	);
	itsDiscFactor		= CreateClonedPtr ( &*	leg.itsDiscFactor		);

	itsFwdRateStartDates= CreateClonedPtr ( &*	leg.itsFwdRateStartDates);
	itsFwdRateEndDates	= CreateClonedPtr ( &*	leg.itsFwdRateEndDates	);
	itsFwdIrRates		= CreateClonedPtr ( &*	leg.itsFwdIrRates		);
	itsIrResetDates		= CreateClonedPtr ( &*	leg.itsIrResetDates		);
	itsIrResetTerms		= CreateClonedPtr ( &*	leg.itsIrResetTerms		);
	itsIrInterTerms		= CreateClonedPtr ( &*	leg.itsIrInterTerms		);
	itsIrAdjConv		= CreateClonedPtr ( &*	leg.itsIrAdjConv		);
	
	itsNumResetDates	= CreateClonedPtr ( &*	leg.itsNumResetDates	);
	itsDemResetDates	= CreateClonedPtr ( &*	leg.itsDemResetDates	);
	itsNumPublishDates	= CreateClonedPtr ( &*	leg.itsNumPublishDates	);
	itsDemPublishDates	= CreateClonedPtr ( &*	leg.itsDemPublishDates	);
	itsInfResetDates	= CreateClonedPtr ( &*	leg.itsInfResetDates	);
	itsInfResetTerms	= CreateClonedPtr ( &*	leg.itsInfResetTerms	);
	itsInfInterTerms	= CreateClonedPtr ( &*	leg.itsInfInterTerms	);
	itsNumCPIRates		= CreateClonedPtr ( &*	leg.itsNumCPIRates		);
	itsDemCPIRates		= CreateClonedPtr ( &*	leg.itsDemCPIRates		);
	itsFwdInfRates		= CreateClonedPtr ( &*	leg.itsFwdInfRates		);
	itsInfAdjConv		= CreateClonedPtr ( &*	leg.itsInfAdjConv		);

	itsIrLeverage		= CreateClonedPtr ( &*	leg.itsIrLeverage		);
	itsInfLeverage		= CreateClonedPtr ( &*	leg.itsInfLeverage		);
	itsConstant			= CreateClonedPtr ( &*	leg.itsConstant			);	
	itsUpMultiple		= CreateClonedPtr ( &*	leg.itsUpMultiple		);
	itsDownMultiple		= CreateClonedPtr ( &*	leg.itsDownMultiple		);	   
	itsRangeUp			= CreateClonedPtr ( &*	leg.itsRangeUp			);
	itsRangeDown		= CreateClonedPtr ( &*	leg.itsRangeDown		);
	itsNotional			= CreateClonedPtr ( &*	leg.itsNotional			);

	itsUpPrice			= CreateClonedPtr ( &*	leg.itsUpPrice			);
	itsDownPrice		= CreateClonedPtr ( &*	leg.itsDownPrice		);

	itsRecOrPay			= leg.itsRecOrPay;
	itsNbFlows			= leg.itsNbFlows;

	itsEpsilon			= leg.itsEpsilon;
	itsNbGaussLeg		= leg.itsNbGaussLeg;
	itsPtGaussLeg		= leg.itsPtGaussLeg;

}

		/********************************************************************************/
		/*																				*/
		/*									VIEWER										*/
		/*																				*/
		/********************************************************************************/
void	ARM_InfCorridorLeg::View		( char* id , FILE* ficOut ){
    FILE* fOut;
    char fOutName[200];
	
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

	CC_Ostringstream	os;
	Ostream<> ss;
	Ostream<string,int> is;
	Ostream<string,double> ds;

	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(55)<<"INFLATION CORRIDOR LEG"<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n";	
	os	<<"INDEX DESCRIPTION\n";
	os	<<"================="<<"\n";

	os	<< ss.ToStream("Index Currency",	itsIrLeg->GetIRIndex()->GetCurrencyUnit()->GetCcyName() );
	os	<< ss.ToStream("Index IR",			ARM_ParamView::GetMappingName(	S_INDEX_TYPES, itsIrLeg->GetIRIndex()->GetIndexType() ));
	os	<< ss.ToStream("Index INF",			ARM_ParamView::GetMappingName(	S_INDEX_TYPES, itsInfLeg->GetIRIndex()->GetIndexType() ));
	os	<< ss.ToStream("INF Leg",			ARM_ParamView::GetMappingName(	S_LEGTYPE, itsInfLeg->GetSwapType() ));   
	os	<< ss.ToStream("INF Interpol",		ARM_ParamView::GetMappingName(	S_INTERPOL_METHOD, itsInfLeg->GetInterpType() ));
	if ( isIrCritera)
		os	<< ss.ToStream("Critera On",	ARM_ParamView::GetMappingName(	S_INDEX_TYPES, itsIrLeg->GetIRIndex()->GetIndexType() ));
	else
		os	<< ss.ToStream("Critera On",	ARM_ParamView::GetMappingName(	S_INDEX_TYPES, itsInfLeg->GetIRIndex()->GetIndexType() ));

	os	<< "\n";
	os	<<"SCHEDULE DESCRIPTION\n";
	os	<<"===================="<<"\n";

	char d[20];

	itsAsOfDate.JulianToStrDate(d);
	os	<< ss.ToStream("As of Date",(string) d);
	itsStartDate.JulianToStrDate(d);
	os	<< ss.ToStream("Start Date",(string) d);
	itsEndDate.JulianToStrDate(d);
	os	<< ss.ToStream("End Date",	(string) d);
	os	<< "\n";
	
	os	<< ss.ToStream("Reset Calendar",	itsIrLeg->GetResetCalName() );
	os	<< ss.ToStream("Reset Frequency",	ARM_ParamView::GetMappingName(	S_FREQUENCY,	itsIrLeg->GetIRIndex()->GetResetFrequency() ));
	os	<< ss.ToStream("Reset Timing",		ARM_ParamView::GetMappingName(	S_TIMING_MOD,	itsIrLeg->GetIRIndex()->GetResetTiming() ));    
	os	<< is.ToStream("Reset Gap",			itsIrLeg->GetIRIndex()->GetResetGap() );
	os	<< is.ToStream("Num Reset Gap",		itsInfLeg->GetNumResetGap() );
	os	<< is.ToStream("Dem Reset Gap",		itsInfLeg->GetDemResetGap() );
  	os	<< "\n";
	os	<< ss.ToStream("Pay Calendar",		itsIrLeg->GetPayCalName() );
	os	<< ss.ToStream("Pay Frequency",		ARM_ParamView::GetMappingName( S_FREQUENCY,		itsIrLeg->GetIRIndex()->GetPayFrequency() ));
	os	<< ss.ToStream("Reset Timing",		ARM_ParamView::GetMappingName( S_TIMING_MOD,	itsIrLeg->GetIRIndex()->GetPayTiming() ));    
	os	<< is.ToStream("Reset Gap",			itsIrLeg->GetIRIndex()->GetPayGap() );
  	os	<< "\n";
	os	<< ss.ToStream("Forward Rule",		ARM_ParamView::GetMappingName( S_FORWARD_RULES,	itsIrLeg->GetIRIndex()->GetFwdRule() ));
	os	<< ss.ToStream("Int Rule",			ARM_ParamView::GetMappingName( S_INTEREST_RULES,itsIrLeg->GetIRIndex()->GetIntRule() )); 
	os	<< ss.ToStream("Stub Rule",			ARM_ParamView::GetMappingName( S_STUB_RULES,	itsIrLeg->GetStubMeth() )); 
	os	<< ss.ToStream("Adj First Rule",	ARM_ParamView::GetMappingName( S_INTEREST_RULES,itsIrLeg->GetAdjStartDateFlag() )); 

	if ( isIrCritera)
		os	<< ss.ToStream("Coupon Day Count",	ARM_ParamView::GetMappingName( S_DAYCOUNT,	itsIrLeg->GetDayCount()) );
	else
		os	<< ss.ToStream("Coupon Day Count",	ARM_ParamView::GetMappingName( S_DAYCOUNT,	itsInfLeg->GetDayCount()) );

	os	<< ss.ToStream("Ir Day Count",		ARM_ParamView::GetMappingName( S_DAYCOUNT,		itsIrLeg->GetIRIndex()->GetDayCount() ));
	os	<< ss.ToStream("Inf Day Count",		ARM_ParamView::GetMappingName( S_DAYCOUNT,		itsInfLeg->GetIRIndex()->GetDayCount() ));
	os	<< ss.ToStream("Receive or Pay",	ARM_ParamView::GetMappingName( S_RECEIVE_PAY,	itsRecOrPay ) );

	os	<< "\n";
	os	<<"MODEL DESCRIPTION\n";
	os	<<"===================="<<"\n";

	os	<< ss.ToStream("Copula Model", ViewModelFeatures() );
	os	<< ds.ToStream("Epsilon Digital Parameter", itsEpsilon );
	os	<< is.ToStream("Nb Point Gauss Legendre", itsNbGaussLeg );
	os	<< ds.ToStream("Integration Symetric Bound", itsPtGaussLeg );

	if ( isComputed ){
	  	os	<< "\n";	
		os	<<"SCHEDULE MANAGER\n";
		os	<<"================"<<"\n";
	
		os  << ViewGeneralFeatures	( );
		os  << ViewInflationFeatures( );
		os  << ViewRateFeatures		( );
		os  << ViewDealFeatures		( );
		os  << ViewPricerFeatures	( );
	}
	fprintf(fOut, "%s" ,(os.str()).c_str() );
	if ( ficOut == NULL )
		fclose(fOut);
}
		/********************************************************************************/
		/*																				*/
		/*								COMPUTE PRICE									*/
		/*																				*/
		/********************************************************************************/

double ARM_InfCorridorLeg::ComputePrice(int ){

	double price = 0.0;

	CptExpFwdRates();

	InitModel();

	CptExpFwdCorrel();
	CptExpFwdInfVol();
	CptExpFwdIrVol();


	price = 0.0;
	for ( int i =0 ;i <itsNbFlows; i++)		StripCompute(i);

	for (  i =0 ;i <itsNbFlows; i++)
		price += (itsDownPrice->Elt(i)-itsUpPrice->Elt(i))*itsNotional->Elt(i)*itsDiscFactor->Elt(i)*itsInterestTerms->Elt(i);

	isComputed = true;
	return itsRecOrPay*price/100;
} 

		/********************************************************************************/
		/*																				*/
		/*						COMPUTE RATE & INFLATION FWD 							*/
		/*																				*/
		/********************************************************************************/

void	ARM_InfCorridorLeg::CptExpFwdRates(){

	ARM_GP_Vector CPIRatio;

	ARM_Model*		model	= GetModel();
	ARM_IRIndex*	tmp		= dynamic_cast<ARM_IRIndex*>(itsInfLeg->GetIRIndex()->Clone() );

	itsIrLeg->SetModel(model);
	itsIrLeg->CptExpectedFwdRates();
	if( dynamic_cast<ARM_InfBSModel*>(model) ){

		ARM_InfBSModel*	pricingMod	= dynamic_cast<ARM_InfBSModel*	>(model);
		ARM_InfIdx*		infIdx		= (ARM_InfIdx*	) (tmp);

		for(int i=0; i<itsNbFlows; ++i )	{

			CPIRatio = pricingMod->CptFwdCPIRatio(	itsNumResetDates->Elt(i), 
													itsDemResetDates->Elt(i),
													itsPaymentDates->Elt(i),
													itsInfLeg->GetInterpType(), 
													itsInfLeg->GetFirstReset(), 
													infIdx);

			itsFwdInfRates	->Elt(i)	= CPIRatio[0] * RateBase ;
			itsNumCPIRates	->Elt(i)	= CPIRatio[1];
			itsDemCPIRates	->Elt(i)	= CPIRatio[2];
			itsInfAdjConv	->Elt(i)	= (1+CPIRatio[0])*CPIRatio[2]/CPIRatio[1];
			itsFwdIrRates	->Elt(i)	= itsIrLeg -> GetFwdRates() -> Elt(i);
			itsDiscFactor   ->Elt(i)	= model->ZeroPrice( 0.0, itsIrLeg->GetYearTerms()->Elt(i), itsInfLeg->GetDiscountingYC() );  
		}

		itsIrLeg->CptExpectedFwdRatesNoAdj();
		for( i=0; i<itsNbFlows; ++i )	
			itsIrAdjConv ->Elt(i)	= itsFwdIrRates	->Elt(i) /itsIrLeg -> GetFwdRates() -> Elt(i);
		itsIrLeg->CptExpectedFwdRates();

		if(tmp) {delete tmp; tmp=NULL;}
	}
	else 
		ARMTHROW(ERR_INVALID_ARGUMENT,"model is of wrong type.. is supposed to price CPI fwd but failed to do so because it does  not inherit from InfFwdModel");	
}



		/********************************************************************************/
		/*																				*/
		/*									VIEW FEATURES								*/
		/*																				*/
		/********************************************************************************/

string	ARM_InfCorridorLeg::ViewGeneralFeatures	(  ){
   
	CC_Ostringstream	os;
	vector<string>		vc(3);
	vector<double>		vd(3);
	char d[20];

	os	<< "\n" <<"General Features :"<<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<" ";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"StartDates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"EndDates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"PayDates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"IntDays";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"IntTerms";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"DiscFactor";
	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"==";
	os << CC_NS(std,setfill('=')) << CC_NS(std,setw)(6*LAG+2)<<" \n";
	
	for ( int i = 0; i < itsNbFlows; ++i )
	{
		( (ARM_Date) itsFlowStartDates->Elt(i)).JulianToStrDate(d);		vc[0] = (string) d;
		( (ARM_Date) itsFlowEndDates->Elt(i)).JulianToStrDate(d);		vc[1] = (string) d;
		( (ARM_Date) itsPaymentDates->Elt(i)).JulianToStrDate(d);		vc[2] = (string) d;
	
		vd[0] =	itsInterestDays	->Elt(i);
		vd[1] =	itsInterestTerms->Elt(i);
		vd[2] =	itsDiscFactor	->Elt(i);
	
	
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<i;
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[0];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[1];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[2];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(0)<<vd[0];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[1];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[2];
		os <<"\n";
	}
	return os.str();
}

  
string	ARM_InfCorridorLeg::ViewInflationFeatures(  ){
  
	CC_Ostringstream	os;
	vector<string>		vc(5);
	vector<double>		vd(6);
	char d[20];

	os	<< "\n" <<"Infla Features :"<<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<" ";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"DemDates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"NumDates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"DemPubDates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"NumPubDates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"ResetDates";	
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"ResetTerms";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"InterTerms";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"FwdDemCPI";	
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"FwdNumCPI";	
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"FwdRates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"FwdAdjConv";
	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed)<< CC_NS(std,setw)(LAG)<<"==";
	os << CC_NS(std,setfill('=')) << CC_NS(std,setw)(11*LAG+2)<<" \n";
	
	for (int i = 0; i < itsNbFlows; ++i )
	{
		( (ARM_Date) itsDemResetDates	->Elt(i)).JulianToStrDate(d);		vc[0] = (string) d;
		( (ARM_Date) itsNumResetDates	->Elt(i)).JulianToStrDate(d);		vc[1] = (string) d;
		( (ARM_Date) itsDemPublishDates	->Elt(i)).JulianToStrDate(d);		vc[2] = (string) d;
		( (ARM_Date) itsNumPublishDates	->Elt(i)).JulianToStrDate(d);		vc[3] = (string) d;
		( (ARM_Date) itsInfResetDates	->Elt(i)).JulianToStrDate(d);		vc[4] = (string) d;
	
		vd[0] =	itsInfResetTerms->Elt(i);
		vd[1] =	itsInfInterTerms->Elt(i);
		vd[2] =	itsDemCPIRates	->Elt(i);
		vd[3] =	itsNumCPIRates	->Elt(i);
		vd[4] =	itsFwdInfRates	->Elt(i);
		vd[5] =	itsInfAdjConv	->Elt(i);
	
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<i;
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[0];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[1];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[2];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[3];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[4];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[0];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[1];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[2];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[3];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[4];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[5];
		os <<"\n";
	}
	return os.str( );
}

string	ARM_InfCorridorLeg::ViewRateFeatures( ){
   
	CC_Ostringstream	os;
	vector<string>		vc(3);
	vector<double>		vd(4);
	char d[20];

	os	<< "\n" <<"Rate Features :"<<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<" ";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"StartDates";	
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"EndDates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"ResetDates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"ResetTerms";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"InterTerms";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"FwdRates";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"FwdAdjConv";
	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed)<< CC_NS(std,setw)(LAG)<<"==";
	os << CC_NS(std,setfill('=')) << CC_NS(std,setw)(7*LAG+2)<<" \n";

	for (int i = 0; i < itsNbFlows; ++i )
	{
		( (ARM_Date) itsFwdRateStartDates->Elt(i)).JulianToStrDate(d);	vc[0] = (string) d;
		( (ARM_Date) itsFwdRateEndDates->Elt(i)).JulianToStrDate(d);	vc[1] = (string) d;
		( (ARM_Date) itsIrResetDates->Elt(i)).JulianToStrDate(d);		vc[2] = (string) d;
	
		vd[0] =	itsIrResetTerms ->Elt(i);
		vd[1] =	itsIrInterTerms ->Elt(i);
		vd[2] =	itsFwdIrRates	->Elt(i);
		vd[3] =	itsIrAdjConv	->Elt(i);
	
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<i;
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[0];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[1];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<vc[2];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[0];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[1];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[2];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[3];
		os <<"\n";
	}
	return os.str( );
}

string	ARM_InfCorridorLeg::ViewDealFeatures(  ){
 

	CC_Ostringstream	os;
	vector<double>		vd(8);

	os	<< "\n" <<"Deal Features :"<<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<" ";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"IrLeverage";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"InfLeverage";	
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Constant";	
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"UpMult";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"DownMult";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"UpRange";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"DownRange";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<"Notional";
	os <<"\n";
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed)<< CC_NS(std,setw)(LAG)<<"==";
	os << CC_NS(std,setfill('=')) << CC_NS(std,setw)(8*LAG+2)<<" \n";
	
	for (int i = 0; i < itsNbFlows; ++i )
	{
		vd[0] =	itsIrLeverage	->Elt(i);
		vd[1] =	itsInfLeverage	->Elt(i);
		vd[2]=	itsConstant		->Elt(i);
		vd[3]=	itsUpMultiple	->Elt(i);
		vd[4]=	itsDownMultiple	->Elt(i);
		vd[5]=	itsRangeUp		->Elt(i);
		vd[6]=	itsRangeDown	->Elt(i);
		vd[7]=	itsNotional		->Elt(i);
	
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<i;
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[0];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[1];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[2];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[3];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[4];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[5];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[6];
		os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(LAG)<<CC_NS(std,setprecision)(4)<<vd[7];
		os <<"\n";
	}
	return os.str( );
}

		/********************************************************************************/
		/*																				*/
		/*										COPULA									*/
		/*																				*/
		/********************************************************************************/

void ARM_Copula::InitGaussLeg( ){
	GaussLegendre_Coefficients c( (int) itsNbGaussLeg, -itsPtGaussLeg, itsPtGaussLeg);
	itsPosit.resize(  (int) itsNbGaussLeg );
	itsWeigth.resize( (int) itsNbGaussLeg );

	for ( int i = 0; i< (int) itsNbGaussLeg; i++){
		itsPosit.Elt(i)	=	c.get_point(i);
		itsWeigth.Elt(i)=	c.get_weight(i)*exp(-0.5*itsPosit.Elt(i)*itsPosit.Elt(i))/sqrt(2.0*PI);
	}
}

ARM_Copula::ARM_Copula( ARM_VolCurvePtr vol1, ARM_VolCurvePtr vol2, map<string,double> & modelParam ){
	itsVol1			= CreateClonedPtr ( &*vol1 );
	itsVol2			= CreateClonedPtr ( &*vol2 );
	itsEpsilon		= modelParam["epsilon"		];
	itsNbGaussLeg	= modelParam["nbGaussLeg"	];
	itsPtGaussLeg	= modelParam["ptGaussLeg"	];

	InitGaussLeg( );
}

void	ARM_Copula::SetVolParam ( const ARM_GP_Vector & param1, const ARM_GP_Vector & param2 ){
	itsParam1	= param1;
	itsParam2	= param2;
}

void	ARM_Copula::SetMarginal(	double (*f1)(ARM_GP_Vector, ARM_VolCurvePtr, double), 
									double (*f2)(ARM_GP_Vector, ARM_VolCurvePtr, double) ) {
	CptVol1 = f1;
	CptVol2 = f2;
}

ARM_Copula::ARM_Copula( const ARM_Copula & cop ){


	itsParam1		= cop.itsParam1;
	itsVol1			= CreateClonedPtr ( &* cop.itsVol1 );
	CptVol1			= cop.CptVol1;

	itsParam2		= cop.itsParam2;
	itsVol2			= CreateClonedPtr ( &* cop.itsVol2 );
	CptVol2			= cop.CptVol2;

	itsPosit		= cop.itsPosit;
	itsWeigth		= cop.itsWeigth;

	itsFwd1			= cop.itsFwd1;
	itsTerms1		= cop.itsTerms1;
	itsStrike1		= cop.itsStrike1;
	itsLev1			= cop.itsLev1;

	itsFwd2			= cop.itsFwd2;
	itsTerms2		= cop.itsTerms2;
	itsStrike2		= cop.itsStrike2;
	itsLev2			= cop.itsLev2;

	itsCorrel		= cop.itsCorrel;
	itsStrike		= cop.itsStrike;
	itsMultiple		= cop.itsMultiple;
	itsConstant		= cop.itsConstant;

	itsEpsilon		= cop.itsEpsilon;
	itsNbGaussLeg	= cop.itsNbGaussLeg;
	itsPtGaussLeg	= cop.itsPtGaussLeg;
}


void	ARM_Copula::SetPayoff	( const ARM_GP_Vector & input){

	itsFwd1			= input.Elt(0);
	itsTerms1		= input.Elt(1);

	itsFwd2			= input.Elt(2);
	itsTerms2		= input.Elt(3);

	itsStrike		= input.Elt(4);
	itsMultiple		= input.Elt(5);

	itsCorrel		= input.Elt(6);
	itsLev1			= input.Elt(7);
	itsLev2			= input.Elt(8);
	itsConstant		= input.Elt(9);

	itsStrike1	= itsMultiple*itsFwd2 + itsStrike;
	itsStrike2	= itsMultiple!=0.0 ? (itsFwd1-itsStrike)/itsMultiple: K_NEW_DOUBLE_TOL;
}


int ARM_GaussCopula::nbStrike = 2;


double	ARM_GaussCopula::Compute (){

	double tmp;
	double weigth;

	double pos1, uni1, F1, W1;
	double pos2, uni2, F2, W2;


	ARM_GP_Vector	vol1(2*nbStrike+1);
	ARM_GP_Vector	mon1(2*nbStrike+1);
	ARM_GP_Vector	vol2(2*nbStrike+1);
	ARM_GP_Vector	mon2(2*nbStrike+1);

	double volMoney1	=	(*CptVol1)(itsParam1, itsVol1, itsFwd1)*sqrt(itsTerms1);
	double volMoney2	=	(*CptVol2)(itsParam2, itsVol2, itsFwd2)*sqrt(itsTerms2);

	for( int i=0; i<2*nbStrike+1;i++){
		mon1.Elt(i)  =	(i-nbStrike);
		vol1.Elt(i)  =	(*CptVol1)( itsParam1, itsVol1,  itsFwd1 * exp( (i-nbStrike) * volMoney1 ) ) * sqrt(fabs(itsTerms1));

		mon2.Elt(i)  =	(i-nbStrike);
		vol2.Elt(i)  =  (*CptVol2)( itsParam2, itsVol2,  itsFwd2 * exp( (i-nbStrike) * volMoney2 ) ) * sqrt(fabs(itsTerms2));
	}

	ARM_SplineDensityFunctor	dist1( mon1, vol1, ARM_SmileViewer::BLACK);
	ARM_SplineDensityFunctor	dist2( mon2, vol2, ARM_SmileViewer::BLACK);


	double value	=	0.0;
	for (  i = 0; i < (int) itsNbGaussLeg; i++){
		pos1		= itsPosit.Elt(i);
		uni1		= NormalCDF(pos1);
		F1			= dist1.Quantile( uni1, itsFwd1 , 1.0 );
		if(F1 > 100 * itsFwd1)	F1		= 0.0;
		W1			= itsWeigth.Elt(i)*exp( 0.5 * pos1 * pos1 );

		if( itsMultiple == 0.0  && itsLev2 == 0.0){
			tmp		=	F1-(itsStrike)>0.0?F1-(itsStrike):0.0 ;
			tmp		-=	F1-(itsStrike+itsEpsilon)>0.0?F1-(itsStrike+itsEpsilon):0.0;
			tmp		/=  itsEpsilon;

			value	+=	itsWeigth.Elt(i)* ( itsLev1*F1 + itsConstant ) * tmp;
		}

		else{	
			for ( int j = 0; j< (int) itsNbGaussLeg; j++){
				pos2	= itsPosit.Elt(j);
				uni2	= NormalCDF(pos2);
				F2		= dist2.Quantile( uni2, itsFwd2 , 1.0 );
				if(	F2	>	100*itsFwd2)	F2	= 0.0;
				W2		= itsWeigth.Elt(j)*exp( 0.5 * pos2 * pos2 );
			
				tmp		=	F1-(itsMultiple * F2+itsStrike)>0.0?F1-(itsMultiple * F2+itsStrike):0.0 ;
				tmp		-=	F1-(itsMultiple * F2 + itsStrike + itsEpsilon)>0.0?F1-(itsMultiple * F2 + itsStrike + itsEpsilon):0.0;
				tmp		/=  itsEpsilon;

				weigth 	= -0.5*( pos1*pos1+pos2*pos2-2*pos1*pos2*itsCorrel)/(1-itsCorrel*itsCorrel); 
				weigth  =	W1 * W2 * exp( weigth )/sqrt(1-itsCorrel*itsCorrel);
		
				value	+= weigth* ( itsLev1*F1 + itsLev2*F2 + itsConstant ) * tmp;
			}		
		}
	}
	return value;
}





CC_END_NAMESPACE()	

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

