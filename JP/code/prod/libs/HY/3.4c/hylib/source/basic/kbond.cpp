#include "kbond.h"
#include "kexception.h"


KBondType::KBondType(const KString &type)
{
	if( type == "ATS_T")
		_bondType = GTO_ATS;         /* Austrian Federal Government Bonds. */
	else if(type == "AUD_T")
		_bondType =GTO_AUD,         /* Australian CGBs. */
	else if(type == "BEF_T")
		_bondType = GTO_BEF,         /* Belgian OLOs and BGBs. */
	else if(type == "CAD_T")
		_bondType = GTO_CAD,         /* Canadian Government Bonds. */
	else if(type == "CHF_T")
		_bondType = GTO_CHF,         /* Swiss SGBs and SGNs. */
	else if(type == "DEM_T")
		_bondType = GTO_DEM,         /* German Bund, Schatze, Bobls, etc. */
	else if(type == "DKK_T")
		_bondType = GTO_DKK,         /* Danish Government Bonds. */
	else if(type == "ECU_T")
		_bondType = GTO_ECU,         /* European Currency Unit bonds. */
	else if(type == "ESP_T")
		_bondType = GTO_ESP,         /* Spanish Bonos and Obligaciones. */
	else if(type == "FIM_T")
		_bondType = GTO_FIM,         /* Finnish Markkas. */
	else if(type == "FRF_T")
		_bondType = GTO_FRF,         /* French BTANs and BTFs. */
	else if(type == "FRF_OAT")
		_bondType = GTO_FRF_OAT,     /* French OATs (accrual is rounded to 3 places). */
	else if(type == "GBP_T")
		_bondType = GTO_GBP,         /* U.K. Gilts, including perpetuals and partly-paids. */
	else if(type == "GBP_IL")
		_bondType = GTO_GBP_IL,      /* U.K. index-linked Gilts. */
	else if(type == "IEP_T")
		_bondType = GTO_IEP,         /* Irish Government Bonds. */
	else if(type == "ITL_T")
		_bondType = GTO_ITL,         /* Italian Government Bonds. */
	else if(type == "JPY_T")
		_bondType = GTO_JPY,         /* Japanese Government Bonds. */
	else if(type == "NLG_T")
		_bondType = GTO_NLG,         /* Dutch State Loans. */
	else if(type == "NZD_T")
		_bondType = GTO_NZD,         /* New Zealand Government Stock. */
	else if(type == "SEK_T")
		_bondType = GTO_SEK,         /* Swedish Government Bonds. */
	else if(type == "USD_T")
		_bondType = GTO_USD,         /* U.S. Treasuries. */
	else if(type == "USD_CORP")
		_bondType = GTO_USD_CORP,    /* U.S. Corporates. */
	else if(type == "USD_MUNI")
		_bondType = GTO_USD_MUNI,    /* U.S. Municipals. */
	else if(type == "ZAR_T")
		_bondType = GTO_ZAR,         /* South African Government Bonds. */
	else if(type == "FRN")
		_bondType = GTO_FRN,         /* A floating rate note that obeys ISMA conventions
								(day count is Act/360). */
	else if(type == "FRN_GBP")
		_bondType = GTO_FRN_GBP,     /* A Euro-sterling FRN (day count is Act/365L). */
	else if(type == "FRCD_GBP")
		_bondType = GTO_FRCD_GBP,    /* A sterling floating rate CD (day count is Act/365). */
	else if(type == "EUROBOND")
		_bondType = GTO_EUROBOND,    /* Standard corporate Eurobond */
	else if(type == "USD_IL")
		_bondType = GTO_USD_IL,      /* U.S. index-linked Treasuries. */
	else if(type == "USD_IL_NPSA")
		_bondType = GTO_USD_IL_NPSA, /* Same as GTO_USD_IL. But instead of using PSA's formula 
						to estimate unknown index value, use future inflation 
						curve */
	else if(type == "PLN_T")
		_bondType = GTO_PLN,          /* Polish Bond with NEW currency re-denomination */
	else if(type == "GRD_T")
		_bondType = GTO_GRD,          /* Greek Government Bonds.*/
	else if(type == "ILS_T")
		_bondType = GTO_ILS,          /* Israeli Government Bonds. */
	else if(type == "HUF_T")
		_bondType = GTO_HUF,          /* Hungarian Government Bonds */ 
	else if(type == "RUB_T")
		_bondType = GTO_RUB,          /* Russian Federation Bonds */
	else if(type == "CZK_T")
		_bondType = GTO_CZK,          /* Czech Republic Bonds */
	else if(type == "EUR")
		_bondType = GTO_EUR,          /* European Union Government Bond denominated in EUROs */
	else if(type == "ATS_EUR")
		_bondType = GTO_ATS_EUR,      /* Austria redenominating in Euro */
	else if(type == "BEF_EUR")
		_bondType = GTO_BEF_EUR,      /* Belgium redenominating in Euro */
	else if(type == "DEM_EUR")
		_bondType = GTO_DEM_EUR,      /* Germany redenominating in Euro */
	else if(type == "ESP_EUR")
		_bondType = GTO_ESP_EUR,      /* Spain redenominating in Euro */
	else if(type == "FIM_EUR")
		_bondType = GTO_FIM_EUR,      /* Finland redenominating in Euro */
	else if(type == "FRF_EUR")
		_bondType = GTO_FRF_EUR,      /* France redenominating in Euro */
	else if(type == "ITL_EUR")
		_bondType = GTO_ITL_EUR,      /* Italy redenominating in Euro */
	else if(type == "NLG_EUR")
		_bondType = GTO_NLG_EUR       /* Netherlands redenominating in Euro */
	else
		KException("Type is not supported!");
}

KBond::KBond()
{
	_tbond = GtoNewEmptyTBond();
	if(!_tbond)
		KException("Error is KBond constructor!");
}
KBond::KBond(const KBond &b)
{
	_tbond = GtoCopyTBond(b._tbond);
	if(!_tbond)
		KException("Error is KBond constructor!");
}
KBond& KBond::operator=(const KBond&b)
{
	if( this != &b)
	{
		_tbond = GtoCopyTBond(b._tbond);
		if(!_tbond)
			KException("Error is KBond constructor!");
	}
	return *this;
}
