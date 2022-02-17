/*------------------------------------------------------------------
C FILE:         kwrapcds.c

CREATED BY:     Neil Yang - Feb 2000

PURPOSE:        DR Wrapper to price credit default swap
---------------------------------------------------------------------- */
//#include "kwrapcds.h"

#include <string.h>

#include "cerror.h"
#include "cfileio.h"
#include "kexception.h"
#include "kwrapgen.h"
#include "readdeal.h"
#include "ddmap.h"
#include "kratecurve.h"
#include "genfunc.h"
#include "instruments.h"
#include "defswap.h"
#include "macros.h"


/*----------------------------------------------------------------------------
FUNCTION:       DrKapCDS

CREATED BY:     Neil Yang - Feb 2000

DESCRIPTION:    Valuation rutine for a credit default swap
----------------------------------------------------------------------------*/
double DrKapCDS()
{
    static char routine[]="DrKapCDS";
    int        status = FAILURE;

    
    char    discFile[] = "disczero.dat";
    char    equityDyn[] = "equity.dyn";
	char    equitySta[] = "equity.sta";
    char    dealFile[] = "hycds_t.dat";

	ErrLog_OnOff = HY_ERR_ON;

	TCurve  *zeroCurve = DrKapWrapGetZC(discFile);
	if(zeroCurve == NULL)
	{
		KException("Error reading disczero.date file.\n",routine);
	}

	TDate  valueDate=zeroCurve->fBaseDate;
	double stockPrice;
	double correlation; // not used
	TCurve *volCurve = NULL;

	if(DrKapWrapReadEquityDyn(equityDyn,
							  valueDate,
							  &stockPrice,
							  &correlation,
							  &volCurve) == FAILURE)
	{
		KException("Error reading equity.dyn file.\n",routine);
	}

	/* get divident info */
	long numPoints;
	TDate  *dividentDates = NULL;
	double *dividentRates = NULL;

	if(DrKapWrapReadEquitySta(equitySta,
							  &numPoints,
							  &dividentDates,
							  &dividentRates) == FAILURE)
	{
		KException("Error reading equity.sta file.\n",routine);
	}

	KDate  *kdividentDates = new KDate[numPoints];
	for(int i =0; i<numPoints;i++)
	{
		kdividentDates[i] = dividentDates[i];
	}

	/* get deal info */
	double  ppy;
	double  beta;
	double  dps;
	double  x;
	double  lim;
	double  lim1;
	double  lim2;
	double  recovery;
	double  notional;
	TDate   maturityDate;
	char   fee[128];

	if(DrKapWrapGetCreditDefaultSwap(dealFile,
									 &ppy,
									 &beta,
									 &dps,
									 &x,
									 &lim,
									 &lim1,
									 &lim2,
									 &recovery,
									 &notional,
									 &maturityDate,
									 fee) == FAILURE)
	{
		KException("Error reading deal.dat file.\n",routine);
	}


	double price;
	try
	{
		double repoValue = 0.0;
		double period = (maturityDate-valueDate)/365.0;
		long ppy2 = (long)(MAX(ppy,ppy*sqrt(period))/period);
		DDMap  *divident = new DDMap(numPoints, kdividentDates, dividentRates);
		KRateCurve  *krepoCurve = new KRateCurve(valueDate,&maturityDate,&repoValue,1,1);
		KRateCurve  *kzeroCurve = new KRateCurve(zeroCurve);
		KRateCurve  *kvolCurve = new KRateCurve(volCurve);
		KRateCurve  *strike1 = new KRateCurve(valueDate, &maturityDate,&lim1,1,1);
		KRateCurve  *strike2 = new KRateCurve(valueDate, &maturityDate,&lim2,1,1);
		LogNormalFunction  *assetProcess = new LogNormalFunction();
		AssetToStockMapping1 *ats = new AssetToStockMapping1(x,lim,dps);
		DefaultSwap   *cds = new DefaultSwap( strike1, strike2, recovery, fee);
		
		/* call pricer */
		KValarray<double> &p = GeneralPricer(stockPrice,
							  divident,
							  krepoCurve,
							  kzeroCurve,
							  kvolCurve,
							  cds,
							  ats,
							  assetProcess,
							  ppy2,
							  beta,
							  valueDate);
		price = p[0];

	}

	catch(KException &e)
	{
		HYErrLog(e.what());
	}

	return price*notional;
}

/*----------------------------------------------------------------------------
FUNCTION:       main
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
CREATED BY:     Neil Yang - Feb 2000

DESCRIPTION:    run valuation of credit default swap
----------------------------------------------------------------------------*/
int main()
{
    int status = FAILURE;
    char    priceFile[] = "price";
    double  pv;
    TFile   *fp   = NULL;   /* File pointer */

    GtoErrMsgOn();
    GtoErrMsgFilePointer(stdout);

	pv = DrKapCDS();
   
    fp = GtoFopen(  priceFile,
                    GTO_FWRITE);
    if(fp IS NULL)
        goto done;

    if(GtoFprintf( fp,"%.12f\n",pv) ISNT FAILURE)
        status = SUCCESS;

done:

    if(fp)GtoFclose( fp);

    return(status);
}


						
