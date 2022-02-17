/* ===============================================================
   FILE_NAME:	utcurrency.c

   PURPOSE:     A few translation functions for currencies   
   =============================================================== */

#include        "utallhdr.h"

typedef struct {
	String	name;
	CcyCode code;
	void*	next;
} CcyCodeListItem;

static CcyCodeListItem* dynamic_ccy_code_list = 0;
static CcyCode next_available_code = DEFAULT_CCY + 1;


/* ---------------------------------------------------------------------------
                                CURRENCIES
   --------------------------------------------------------------------------- */

/*************************************************************************************/

void delete_dynamic_ccy_codes(void)
{
	CcyCodeListItem *item = dynamic_ccy_code_list,*prev;
	next_available_code = DEFAULT_CCY + 1;
	while(item)
	{
		srt_free(item->name);
		prev = item;
		item = item->next;
		srt_free(prev);
	}
	dynamic_ccy_code_list=0;
}

/*************************************************************************************/

Err ccy_get_or_create_code(String cs, CcyCode *ccy)
{
	strupper(cs);
	strip_white_space(cs);
	
	if(interp_default_ccy_string(cs,ccy))
	{
		/* not found - add a new one */
		CcyCodeListItem *item, *prev;

		if(dynamic_ccy_code_list)
		{
			item = dynamic_ccy_code_list;
			do {
				prev = item;
				if(strcmp(item->name,cs)==0)
				{
					/* Found the ccy */
					*ccy = item->code;
					return 0;
				}
				item = item->next;

			} while(item);
			/* Add a new element */
			item = (CcyCodeListItem *)srt_calloc(1,sizeof(CcyCodeListItem));
			item->name = (String) malloc(strlen(cs)+1);
			strcpy(item->name,cs);
			item->code = next_available_code;
			*ccy=next_available_code;
			next_available_code++;
			item->next = 0;
			prev->next = item;
		}
		else
		{
			/* Add a new element */
			item = (CcyCodeListItem *)srt_calloc(1,sizeof(CcyCodeListItem));
			item->name = (String) malloc(strlen(cs)+1);
			strcpy(item->name,cs);
			next_available_code = DEFAULT_CCY + 1;
			item->code = next_available_code;
			*ccy=next_available_code;
			next_available_code++;
			item->next = 0;
			dynamic_ccy_code_list = item;
		}
		return 0;

	}
	else
		return 0;
}

/*************************************************************************************/

Err interp_ccy_string(String cs, CcyCode *ccy)
{
	strupper(cs);
	strip_white_space(cs);

	if(interp_default_ccy_string(cs,ccy))
	{
		CcyCodeListItem *item = dynamic_ccy_code_list;
		while(item)
		{
			if(strcmp(cs,item->name)==0)
			{
				*ccy = item->code;
				return 0;
			}
			item = item->next;
		}
		/*
		Special case: if not recognised above but three-character code assume
		default currency...
		*/
		if (strlen(cs) == 3) {*ccy = DEFAULT_CCY; return 0;}

		return serror("unknown currency: %s",cs);
	}
	else
		return 0;
}

/*************************************************************************************/
/* Added some new currencies May 21, 2002 */

int interp_default_ccy_string(String cs, CcyCode *ccy)
{
	strupper(cs);
	strip_white_space(cs);

	if(!strcmp(cs,"AED")){*ccy = AED; return 0;}   
	if(!strcmp(cs,"ARS")){*ccy = ARS; return 0;}   
/*	if(!strcmp(cs,"ATS")){*ccy = ATS; return 0;} */
	if(!strcmp(cs,"AUD")){*ccy = AUD; return 0;}
/*	if(!strcmp(cs,"BEF")){*ccy = BEF; return 0;}*/
	if(!strcmp(cs,"BGN")){*ccy = BGN; return 0;}  
	if(!strcmp(cs,"BHD")){*ccy = BHD; return 0;}       
	if(!strcmp(cs,"BRL")){*ccy = BRL; return 0;}   
	if(!strcmp(cs,"BYR")){*ccy = BYR; return 0;}       
	if(!strcmp(cs,"CAD")){*ccy = CAD; return 0;}
	if(!strcmp(cs,"CHF")){*ccy = CHF; return 0;}
	if(!strcmp(cs,"CLP")){*ccy = CLP; return 0;}   
	if(!strcmp(cs,"CNY")){*ccy = CNY; return 0;}       
	if(!strcmp(cs,"COO")){*ccy = COO; return 0;}   
	if(!strcmp(cs,"COP")){*ccy = COP; return 0;}   
	if(!strcmp(cs,"CZK")){*ccy = CZK; return 0;}       
	if(!strcmp(cs,"DKK")){*ccy = DKK; return 0;}
/*	if(!strcmp(cs,"DM")) {*ccy = DEM; return 0;}
	if(!strcmp(cs,"DEM")){*ccy = DEM; return 0;}
	if(!strcmp(cs,"ECU")){*ccy = XEU; return 0;} */
	if(!strcmp(cs,"EEK")){*ccy = EEK; return 0;}       
	if(!strcmp(cs,"EGP")){*ccy = EGP; return 0;}       
/*	if(!strcmp(cs,"ESB")){*ccy = ESB; return 0;}
	if(!strcmp(cs,"ESP")){*ccy = ESP; return 0;}*/
	if(!strcmp(cs,"EUR")){*ccy = EUR; return 0;}
/*	if(!strcmp(cs,"FFR")){*ccy = FFR; return 0;}
	if(!strcmp(cs,"FIM")){*ccy = FIM; return 0;}       
	if(!strcmp(cs,"FRF")){*ccy = FFR; return 0;}*/
	if(!strcmp(cs,"GBP")){*ccy = GBP; return 0;}
/*	if(!strcmp(cs,"GRD")){*ccy = GRD; return 0;}       */
	if(!strcmp(cs,"HKD")){*ccy = HKD; return 0;}
	if(!strcmp(cs,"HRK")){*ccy = HRK; return 0;}
	if(!strcmp(cs,"HUF")){*ccy = HUF; return 0;}
	if(!strcmp(cs,"IDR")){*ccy = IDR; return 0;}       
/*	if(!strcmp(cs,"IEP")){*ccy = IEP; return 0;}       */
	if(!strcmp(cs,"ILS")){*ccy = ILS; return 0;}
	if(!strcmp(cs,"INO")){*ccy = INO; return 0;}    // P McCallum: 23.10.2003  
	if(!strcmp(cs,"INR")){*ccy = INR; return 0;}       
/*	if(!strcmp(cs,"ITL")){*ccy = ITL; return 0;}*/
	if(!strcmp(cs,"JPB")){*ccy = JPB; return 0;}       
	if(!strcmp(cs,"JPY")){*ccy = JPY; return 0;}
	if(!strcmp(cs,"KRO")){*ccy = KRO; return 0;}
	if(!strcmp(cs,"KRW")){*ccy = KRW; return 0;}
	if(!strcmp(cs,"KWD")){*ccy = KWD; return 0;}       
	if(!strcmp(cs,"KZT")){*ccy = KZT; return 0;}
	if(!strcmp(cs,"LBP")){*ccy = LBP; return 0;}       
	if(!strcmp(cs,"LTL")){*ccy = LTL; return 0;}       
/*	if(!strcmp(cs,"LUF")){*ccy = LUF; return 0;}*/
	if(!strcmp(cs,"LVL")){*ccy = LVL; return 0;}       
	if(!strcmp(cs,"MAD")){*ccy = MAD; return 0;}       
	if(!strcmp(cs,"MXN")){*ccy = MXN; return 0;}       
	if(!strcmp(cs,"MYR")){*ccy = MYR; return 0;}       
/*	if(!strcmp(cs,"NLG")){*ccy = NLG; return 0;}*/
	if(!strcmp(cs,"NOK")){*ccy = NOK; return 0;}       
	if(!strcmp(cs,"NZD")){*ccy = NZD; return 0;}
	if(!strcmp(cs,"OMR")){*ccy = OMR; return 0;}       
	if(!strcmp(cs,"PEN")){*ccy = PEN; return 0;}       
	if(!strcmp(cs,"PHP")){*ccy = PHP; return 0;}       
	if(!strcmp(cs,"PKR")){*ccy = PKR; return 0;}       
	if(!strcmp(cs,"PLN")){*ccy = PLN; return 0;}   
/*	if(!strcmp(cs,"PTE")){*ccy = PTE; return 0;}  */     
	if(!strcmp(cs,"QAR")){*ccy = QAR; return 0;}       
	if(!strcmp(cs,"ROL")){*ccy = ROL; return 0;}       
	if(!strcmp(cs,"RUB")){*ccy = RUB; return 0;}
	if(!strcmp(cs,"SAR")){*ccy = SAR; return 0;}       
	if(!strcmp(cs,"SEK")){*ccy = SEK; return 0;}
	if(!strcmp(cs,"SGD")){*ccy = SGD; return 0;}       
	if(!strcmp(cs,"SIT")){*ccy = SIT; return 0;}       
	if(!strcmp(cs,"SKK")){*ccy = SKK; return 0;}       
	if(!strcmp(cs,"THB")){*ccy = THB; return 0;}       
	if(!strcmp(cs,"THO")){*ccy = THO; return 0;}       
	if(!strcmp(cs,"TND")){*ccy = TND; return 0;}       
	if(!strcmp(cs,"TRL")){*ccy = TRL; return 0;}
	if(!strcmp(cs,"TWD")){*ccy = TWD; return 0;}       
	if(!strcmp(cs,"TWO")){*ccy = TWO; return 0;}       
	if(!strcmp(cs,"UAH")){*ccy = UAH; return 0;}
	if(!strcmp(cs,"USD")){*ccy = USD; return 0;}
	if(!strcmp(cs,"UTF")){*ccy = UTF; return 0;} // Stanley Mrose: 12.02.2003
	if(!strcmp(cs,"VEB")){*ccy = VEB; return 0;}   
	if(!strcmp(cs,"XAU")){*ccy = XAU; return 0;}
/*	if(!strcmp(cs,"XEU")){*ccy = XEU; return 0;}*/
	if(!strcmp(cs,"ZAR")){*ccy = ZAR; return 0;}       
	
	/* not a hard-coded value */
	return 1;
	
}

/*************************************************************************************/

/* Copy the output string - since if ccy dynamic it can be deleted */
Err translate_ccy(String *cs, CcyCode ccy)
{
	strupper( *cs );
	strip_white_space( *cs );

	if(translate_default_ccy(cs,ccy))
	{
		CcyCodeListItem *item = dynamic_ccy_code_list;
		while(item)
		{
			if(item->code==ccy)
			{
				*cs = item->name;
				return 0;
			}
			item = item->next;
		}

		/* Special case: default currency */
		if(ccy == DEFAULT_CCY) {*cs = "CCY"; return 0;}
		
		return "translate_ccy: unknown currency.";
	}
	else
		return 0;
}

/*************************************************************************************/

int translate_default_ccy(String *cs, CcyCode ccy)
{

	strupper( *cs);					     
	strip_white_space( *cs);


	if(ccy == AED){*cs = "AED"; return 0;}   
	if(ccy == ARS){*cs = "ARS"; return 0;}   
/*	if(ccy == ATS){*cs = "ATS"; return 0;} */
	if(ccy == AUD){*cs = "AUD"; return 0;}
/*	if(ccy == BEF){*cs = "BEF"; return 0;}*/
	if(ccy == BGN){*cs = "BGN"; return 0;}  
	if(ccy == BHD){*cs = "BHD"; return 0;}       
	if(ccy == BRL){*cs = "BRL"; return 0;}   
	if(ccy == BYR){*cs = "BYR"; return 0;}       
	if(ccy == CAD){*cs = "CAD"; return 0;}
	if(ccy == CHF){*cs = "CHF"; return 0;}
	if(ccy == CLP){*cs = "CLP"; return 0;}   
	if(ccy == CNY){*cs = "CNY"; return 0;}       
	if(ccy == COO){*cs = "COO"; return 0;}   
	if(ccy == COP){*cs = "COP"; return 0;}   
	if(ccy == CZK){*cs = "CZK"; return 0;}       
	if(ccy == DKK){*cs = "DKK"; return 0;}
/*	if(ccy == DM) {*sy "= "DEM; return 0;}
	if(ccy == DEM){*cs = "DEM"; return 0;}
	if(ccy == ECU){*cs = "XEU"; return 0;} */
	if(ccy == EEK){*cs = "EEK"; return 0;}       
	if(ccy == EGP){*cs = "EGP"; return 0;}       
/*	if(ccy == ESB){*cs = "ESB"; return 0;}
	if(ccy == ESP){*cs = "ESP"; return 0;}*/
	if(ccy == EUR){*cs = "EUR"; return 0;}
/*	if(ccy == FFR){*cs = "FFR"; return 0;}
	if(ccy == FIM){*cs = "FIM"; return 0;}       
	if(ccy == FRF){*cs = "FFR"; return 0;}*/
	if(ccy == GBP){*cs = "GBP"; return 0;}
/*	if(ccy == GRD){*cs = "GRD"; return 0;}       */
	if(ccy == HKD){*cs = "HKD"; return 0;}
	if(ccy == HRK){*cs = "HRK"; return 0;}
	if(ccy == HUF){*cs = "HUF"; return 0;}
	if(ccy == IDR){*cs = "IDR"; return 0;}       
/*	if(ccy == IEP){*cs = "IEP"; return 0;}       */
	if(ccy == ILS){*cs = "ILS"; return 0;}
	if(ccy == INR){*cs = "INR"; return 0;}       
/*	if(ccy == ITL){*cs = "ITL"; return 0;}*/
	if(ccy == JPB){*cs = "JPB"; return 0;}       
	if(ccy == JPY){*cs = "JPY"; return 0;}
	if(ccy == KRO){*cs = "KRO"; return 0;}
	if(ccy == KRW){*cs = "KRW"; return 0;}
	if(ccy == KWD){*cs = "KWD"; return 0;}       
	if(ccy == KZT){*cs = "KZT"; return 0;}
	if(ccy == LBP){*cs = "LBP"; return 0;}       
	if(ccy == LTL){*cs = "LTL"; return 0;}       
/*	if(ccy == LUF){*cs = "LUF"; return 0;}*/
	if(ccy == LVL){*cs = "LVL"; return 0;}       
	if(ccy == MAD){*cs = "MAD"; return 0;}       
	if(ccy == MXN){*cs = "MXN"; return 0;}       
	if(ccy == MYR){*cs = "MYR"; return 0;}       
/*	if(ccy == NLG){*cs = "NLG"; return 0;}*/
	if(ccy == NOK){*cs = "NOK"; return 0;}       
	if(ccy == NZD){*cs = "NZD"; return 0;}
	if(ccy == OMR){*cs = "OMR"; return 0;}       
	if(ccy == PEN){*cs = "PEN"; return 0;}       
	if(ccy == PHP){*cs = "PHP"; return 0;}       
	if(ccy == PKR){*cs = "PKR"; return 0;}       
	if(ccy == PLN){*cs = "PLN"; return 0;}   
/*	if(ccy == PTE){*cs = "PTE"; return 0;}  */     
	if(ccy == QAR){*cs = "QAR"; return 0;}       
	if(ccy == ROL){*cs = "ROL"; return 0;}       
	if(ccy == RUB){*cs = "RUB"; return 0;}
	if(ccy == SAR){*cs = "SAR"; return 0;}       
	if(ccy == SEK){*cs = "SEK"; return 0;}
	if(ccy == SGD){*cs = "SGD"; return 0;}       
	if(ccy == SIT){*cs = "SIT"; return 0;}       
	if(ccy == SKK){*cs = "SKK"; return 0;}       
	if(ccy == THB){*cs = "THB"; return 0;}       
	if(ccy == THO){*cs = "THO"; return 0;}       
	if(ccy == TND){*cs = "TND"; return 0;}       
	if(ccy == TRL){*cs = "TRL"; return 0;}
	if(ccy == TWD){*cs = "TWD"; return 0;}       
	if(ccy == TWO){*cs = "TWO"; return 0;}       
	if(ccy == UAH){*cs = "UAH"; return 0;}
	if(ccy == USD){*cs = "USD"; return 0;}
	if(ccy == UTF){*cs = "UTF"; return 0;} // Stanley Mrose: 12.02.2003
	if(ccy == VEB){*cs = "VEB"; return 0;}   
	if(ccy == XAU){*cs = "XAU"; return 0;}
/*	if(ccy == XEU){*cs = "XEU"; return 0;}*/
	if(ccy == ZAR){*cs = "ZAR"; return 0;}       


	
	return 1;
}

