/*
 * $Log: armrpc.h,v $
 * Revision 1.87  2003/12/16 17:59:46  arm
 * RPC_VERSION: commented
 *
 * Revision 1.86  2003/08/01 15:19:41  jpriaudel
 * Ajout de CalibrationFrmModel
 *
 * Revision 1.85  2003/08/01 13:56:38  mab
 * Added : RPC_CREATE_FRMMARKOVTREE RPC_SET_FRMMARKOVTREE
 * RPC_CREATE_BOOTSCALIB_FRMMODEL RPC_SET_BOOTSCALIB_FRMMODEL
 *
 * Revision 1.84  2003/05/09 09:58:14  mab
 * Added : RPC_CREATE_BUMPCURVE , RPC_SET_BUMPCURVE
 * y
 *
 * Revision 1.83  2003/03/28 16:40:44  mab
 * New RPC_CREATESET_FIXEDRATES RPC_SETSET_FIXEDRATES
 *
 * Revision 1.82  2003/02/11 17:31:43  mab
 * Added : RPC_TRIANGULAR_INTERPOL RPC_CREATE_HWSIGVARFROMANA
 * RPC_SET_HWSIGVARFROMANA RPC_CREATE_HWSIGVARCALIB RPC_SET_HWSIGVARCALIB
 * RPC_CREATE_HWSIGVAR_CALIBRATOR RPC_SET_HWSIGVAR_CALIBRATOR RPC_CREATE_CALIBRATE_BSSMILED
 * RPC_SET_CALIBRATE_BSSMILED RPC_CREATE_GET_SABR_SIG_RHO_NU RPC_SET_GET_SABR_SIG_RHO_NU
 *
 * Revision 1.81  2002/11/15 11:07:26  mab
 * Added :  RPC_GET_SABR_SIGMA RPC_CREATEPORTOPTION RPC_SETPORTOPTION
 * RPC_CUST_FIRST_PERIOD RPC_DISP_SCHED_DATES RPC_DISP_SCHED_VALUES
 *
 * Revision 1.80  2002/10/11 09:26:43  mab
 * Added :  RPC_BSTHETA RPC_BSGAMMA
 *
 * Revision 1.79  2002/10/09 09:50:53  mab
 * Added : RPC_CREATEBUMPVOLATILITY , RPC_SETBUMPVOLATILITY
 *
 * Revision 1.78  2002/09/24 13:29:03  mab
 * Added : RPC_CREATETOYZCSWAPINT , RPC_SETTOYZCSWAPINT
 *
 * Revision 1.77  2002/09/17 16:52:24  mab
 * Added : #define RPC_GETSPOTDAYSFROMCCY          512
 *
 * Revision 1.76  2002/07/29 13:36:56  mab
 * Added : RPC_BSOPTION RPC_BSDELTA RPC_BSVEGA
 *
 * Revision 1.75  2002/07/18 13:59:59  mab
 * RPC_GETUNDPRICE added
 *
 * Revision 1.74  2002/07/08 09:47:35  mab
 * Rajout de : RPC_CREATE_BSCORRMODEL, RPC_SET_BSCORRMODEL
 * RPC_CREATE_SPREADOPTION, RPC_SET_SPREADOPTION
 *
 * Revision 1.73  2002/07/02 13:42:41  mab
 * Added : RPC_CREATEFRMTREE_AUTO_B
 * and
 * RPC_SETFRMTREE_AUTO_B
 *
 * Revision 1.72  2002/04/26 10:25:50  mab
 * added : RPC_ADDPERIOD RPC_INTERPOL RPC_DISCOUNTPRICEREFVAL
 *
 * Revision 1.71  2002/03/28 13:24:45  mab
 * Rajout de : RPC_CREATEFRMTREE_AUTO_G
 * RPC_SETFRMTREE_AUTO_G
 * RPC_CREATEFRMLSMC_AUTO_G
 * RPC_SETFRMLSMC_AUTO_G
 *
 * Revision 1.70  2002/02/14 10:52:46  mab
 * Rajout : RPC_CREATEBSSMILEDMODEL ,
 * RPC_SETBSSMILEDMODEL
 * RPC_GEN_AMORT
 * RPC_SET_GEN_AMORT
 *
 * Revision 1.69  2002/02/04 10:03:41  mab
 * Rajout de : RPC_ADDMONTHS , RPC_ADDYEARS
 *
 * Revision 1.68  2002/01/23 16:52:07  mab
 * rajout de :
 * RPC_GETDEFIDXFROMCCY , RPC_GETPAYCALNAME , RPC_GETCCYNAME
 *
 * Revision 1.67  2001/12/04 13:00:37  mab
 * Rajout de : RPC_BETWEENDATES
 *
 * Revision 1.66  2001/11/21 14:21:06  mab
 * Rajout de : RPC_GETFXVOLFROMSUMMIT et RPC_SETFXVOLFROMSUMMIT
 *
 * Revision 1.65  2001/10/04 09:35:14  abizid
 * Ameliorations
 *
 * Revision 1.64  2001/09/28 13:12:10  mab
 * rajout de : RPC_BMCFRM2CR_CREATE , RPC_BMCFRM2CR_SET
 *
 * Revision 1.63  2001/09/28 12:31:48  mab
 * Rajout de : RPC_BMCFRM_CREATE RPC_BMCFRM_SET
 *
 * Revision 1.62  2001/08/13 14:35:10  mab
 * Rajout de : RPC_CREATEVOLCUBE, RPC_SETVOLCUBE
 *
 * Revision 1.61  2001/08/06 12:41:55  abizid
 * Ajout interface XBSFX
 *
 * Revision 1.60  2001/07/23 13:03:36  mab
 * Rajout de : Vol. Cube
 *
 * Revision 1.59  2001/06/21 12:27:47  smysona
 * nouveau constructeur de smoothcurve
 *
 * Revision 1.58  2001/06/18 10:58:26  sgasquet
 * Ajout RPC_CREATE_OPTIONALACCRUALZC et RPC_SET_OPTIONALACCRUALZC
 *
 * Revision 1.57  2001/05/07 13:12:10  mab
 * Rajout de : RPC_DISPATCH_PID
 *
 * Revision 1.56  2001/04/23 07:59:09  vberger
 * ajout des variables de dëfinition pour les "reverse calendar"
 * RPC_CREATE_REVERSE_CALENDAR et RPC_SET_REVERSE_CALENDAR
 *
 * Revision 1.55  2001/04/20 15:30:21  abizid
 * Ajout RPC_SMILEMCRN... et RPC_INSTLOGDEC...
 *
 * Revision 1.54  2001/03/19 11:09:25  abizid
 * Ajout PFINSTLOGDEC Fitting
 *
 * Revision 1.53  2001/03/19 10:58:59  smysona
 * Ajout des LSMC_AUTO,
 * CRF -> REVERSEFLOAT
 *
 * Revision 1.52  2001/03/12 17:17:55  mab
 * Rajout de : RPC_CREATE_VARFIX_SWOPT, RPC_SET_VARFIX_SWOPT
 * RPC_CREATE_CRF, RPC_SET_CRF
 *
 * Revision 1.51  2001/03/02 19:05:44  smysona
 * XCcy Set
 *
 * Revision 1.50  2001/02/22 20:02:15  smysona
 * ajout RPC_XCCYADJUST
 *
 * Revision 1.49  2001/02/02 08:21:21  sgasquet
 * Rajout de : RPC_FRMSHORTRATEVOLS
 *
 * Revision 1.48  2001/01/30 15:36:04  smysona
 * Ajout des FRMTREE
 *
 * Revision 1.47  2001/01/19 18:41:09  mab
 * *** empty log message ***
 *
 * Revision 1.46  2001/01/05 18:04:59  sgasquet
 * Ajout RPC_CREATE_MATCAPFLOOR et RPC_SET_MATCAPFLOOR
 *
 * Revision 1.45  2000/12/08 21:43:37  mab
 * Ajout Restrikable
 *
 * Revision 1.44  2000/11/23 16:56:59  mab
 * Ajout interface FIXED et CMS
 *
 * Revision 1.43  2000/11/10 17:42:42  smysona
 * FRM integration
 *
 * Revision 1.42  2000/11/08 14:21:14  mab
 * Ajout CorridorLeg
 *
 * Revision 1.41  2000/10/20 09:29:56  mab
 * Ajout Calibrage + Createur MC_FROM_ANA pour LogDec
 *
 * Revision 1.40  2000/10/13 19:26:19  mab
 * *** empty log message ***
 *
 * Revision 1.39  2000/09/28 14:02:18  mab
 * Rajout RPC_...SPREADLEG
 *
 * Revision 1.38  2000/07/21 16:29:10  mab
 * Integration fonction correl.
 *
 * Revision 1.37  2000/07/19 09:49:51  mab
 * Amelioration HistoVol
 *
 * Revision 1.36  2000/07/11 17:12:07  mab
 * Ajout constructeur Log Dec
 *
 * Revision 1.35  2000/06/22 15:51:59  mab
 * adding the Smooth Interpolation
 *
 * Revision 1.34  2000/04/13 12:16:34  ypilchen
 * rajout de RPC_PFHW2FSIGVARGLOBMIN et RPC_SETPFHW2FSIGVARGLOBMIN par jfg
 * >
 *
 * Revision 1.33  2000/02/29 13:40:28  sgasquet
 * ajout modele lognormal decale
 *
 * Revision 1.32  2000/02/16 15:01:46  mab
 * Rajout de Ratchet et Sticky
 *
 * Revision 1.31  2000/02/02 08:46:18  mab
 * Rajout de : RPC_ADJUSTBUSDATE
 *
 * Revision 1.30  2000/01/11 14:37:46  vberger
 * prise en compte d'un nouveau constructeur pour les reversecoupons
 *
 * Revision 1.29  1999/11/24 16:57:31  vberger
 * creation de variables pour les reverse reversecoupon et l'arbre quantoHW
 *
 * Revision 1.28  1999/11/15 08:22:56  sgasquet
 * Ajout fonctions modele analyt et mc hw2f sigvar
 *
 * Revision 1.27  1999/10/01 08:42:36  mab
 * Rajout de : RPC_SALESFILTER
 *
 * Revision 1.26  1999/08/31 10:31:24  mab
 * Rajout de RPC_SUMMITVOLMATRIX
 *
 * Revision 1.25  1999/08/16 13:29:48  mab
 * Rajout de : RPC_PFHW2FPENALFIT , RPC_SETPFHW2FPENALFIT
 *
 * Revision 1.24  1999/07/27 16:11:36  mab
 * Rajout de : RPC_GETFWDRATESHISTO
 *
 * Revision 1.23  1999/06/15 13:23:03  mab
 * Rajout de : RPC_GNPV_FILTER
 *
 * Revision 1.22  1999/05/21 12:02:35  mab
 * Rajout de : RPC_FXCONVERT
 *
 * Revision 1.21  1999/05/04 14:36:53  mab
 * Rajout de : define RPC_GETINITIALCURVEFROMSUMMIT   130
 *
 * Revision 1.20  1999/04/26 14:21:12  mab
 * Rajout de : #define RPC_IMPLIED_SPREAD_WITH2MODELS  461
 *
 * Revision 1.19  1999/03/29 11:36:34  vberger
 * Ajout de deux macros modele BS et une macro HW sigma variable pour le pricing des quantos
 *
 * Revision 1.18  1999/03/22 18:35:57  sgasquet
 * Ajout RPC_CREATE_CFEXOSWAPTION et RPC_SET_CFEXOSWAPTION
 *
 * Revision 1.17  1999/03/10 17:23:12  nicolasm
 * Ajout interface EXOFLEXCF
 *
 * Revision 1.16  1999/03/04 13:16:16  ypilchen
 * Rajout de RPC_GETFWDRATESMATRIX, RPC_GETEXPIRY
 *
 * Revision 1.15  1999/02/26 14:55:04  ypilchen
 * Rajout de #define RPC_CAPFLOOREXPIRY
 *
 * Revision 1.14  1999/02/25 11:53:04  ypilchen
 * Rajout de RPC_CREATEDFBSMODEL, RPC_SETDFBSMODEL
 *
 * Revision 1.13  1999/02/18 17:03:08  nicolasm
 * Ajout RPC_CREATEZCSWAPFUTINT et RPC_SETZCSWAPFUTINT
 *
 * Revision 1.12  1999/02/02 11:13:57  ypilchen
 * Rajout de C_NPV_FILTER
 *
 * Revision 1.11  1999/01/28 16:59:25  nicolasm
 * Ameliorations
 *
 * Revision 1.10  1999/01/19 14:01:40  nicolasm
 * Ajout C_CAPLETPRICE
 *
 * Revision 1.9  1999/01/15 18:00:18  nicolasm
 * Ajout MC HW sig var
 *
 * Revision 1.8  1998/12/23 17:49:25  ypilchen
 * Remontee HW analytique SigVar et sa Calibration
 *
 * Revision 1.7  1998/12/17 16:06:44  nicolasm
 * Ajout de #define ARM_FRF_CCY_OBJECT  -11112
 *
 * Revision 1.6  1998/12/17 15:49:07  ypilchen
 * Rajout de RPC_SETHWSIGVARANALYTIC, RPC_SETHWSIGVARANALYTIC
 *
 * Revision 1.5  1998/12/07 15:06:56  ypilchen
 * Rajout de RPC_GETDEFAULTCOUNTRY, RPC_SETDEFAULTCOUNTRY
 *
 * Revision 1.4  1998/11/30 15:10:10  nicolasm
 * Ajout RPC_CREATE_IASEC RPC_SET_IASEC RPC_CREATE_IA3LEVREFVAL RPC_SET_IA3LEL
 *
 * Revision 1.3  1998/11/27 14:55:48  ypilchen
 * Rajout de RPC_CREATEHWSIGVAR , RPC_SETHWSIGVAR
 *
 * Revision 1.2  1998/11/24 18:33:18  nicolasm
 * Ajout Modele MonteCarloFN
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE           : armrpc.h                                                  */
/*                                                                            */
/* DESCRIPTION    : Applix RPC Requests codes                                 */
/*                                                                            */
/* DATE           : Tue Aug 20 1996                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
 



#define ARM_NULL_OBJECT -11111
#define ARM_FRF_CCY_OBJECT  -11112


/*---- GLOBAL ----*/

#define RPC_DISPATCH_PID             8

#define RPC_FREE_ALL                 9

#define RPC_TEST                    10

#define RPC_CLEANMEMORY             11

#define RPC_SETPRICE                12

#define RPC_PRINT                   13

#define RPC_FREEOBJECT              14

#define RPC_EXIT                    15

#define RPC_PRICE                   16

#define RPC_FORWARD_PRICE           17

#define RPC_ACCRUED                 18

#define RPC_BSSPOT                  19

#define RPC_SENSITIVITY             20

#define RPC_CVSENSITIVITY           21

#define RPC_SETMARKETPRICE          22

#define RPC_VIEW                    23

#define RPC_SCREDIT                 24

#define RPC_SETAMOUNT               25

/* Pb on XLL PC
#define RPC_VERSION                 26
*/

#define RPC_SUMMIT_VALUE_TRADE      27

#define RPC_SUMMIT_VALUE_GREEKS     28

#define RPC_SUMMIT_SWAP_INFO        29

#define RPC_EXPROBA                 30

#define RPC_GETDEFAULTCOUNTRY       31

#define RPC_SETDEFAULTCOUNTRY       32

#define RPC_CPARALLELSHIFT          33

#define RPC_SPARALLELSHIFT          34

#define RPC_NPV_FILTER              35

#define RPC_GETFWDRATESMATRIX       36

#define RPC_GETEXPIRY               37

#define RPC_FXCONVERT               38

#define RPC_GNPV_FILTER             39

#define RPC_GETFWDRATESHISTO        40

#define RPC_EVALSUMMITASSET         41

#define RPC_SUMMITVOLMATRIX         42

#define RPC_SALESFILTER             43

#define RPC_ASOFVOLRATE             44

#define RPC_HISTORICALVOL           45

#define RPC_HISTOCORREL             46

#define RPC_XCCYADJUST              47

#define RPC_SETXCCYADJUST           48

#define RPC_BETWEENDATES            49

#define RPC_ADDPERIOD               50 

#define RPC_INTERPOL                51

#define RPC_DISCOUNTPRICEREFVAL     52

#define RPC_CREATECLONEDUPDNOTIONAL 53

#define RPC_SETCLONEDUPDNOTIONAL    54 

#define RPC_TRIANGULAR_INTERPOL     55


/*---- CURVES ----*/

#define RPC_CREATEZEROCURVELIN      100

#define RPC_SETZEROCURVELIN         101

#define RPC_CREATEZEROCURVESPLI     102

#define RPC_SETZEROCURVESPLI        103

#define RPC_CREATEZEROCURVEVSK      104

#define RPC_SETZEROCURVEVSK         105

#define RPC_DISCOUNTPRICE           106
    
#define RPC_DISCOUNTYIELD           107

#define RPC_FORWARDPRICE            108

#define RPC_FORWARDYIELD            109

#define RPC_GETZCFROMSUMMIT         110

#define RPC_SETCREATEDZEROCURVEFROMSUMMIT    111

#define RPC_CREATEZEROFLAT          112
 
#define RPC_SETZEROFLAT             113

#define RPC_CREATEZCSWAPINT         114 

#define RPC_SETZCSWAPINT            115 

#define RPC_CREATEZCCASHINT         116 

#define RPC_SETZCCASHINT            117

#define RPC_CREATESPREADCURVE       118

#define RPC_SETSPREADCURVE          119

#define RPC_CREATESPLICUBCURVE      120

#define RPC_SETSPLICUBCURVE         121

#define RPC_CREATEZCTAMINT          122

#define RPC_SETZCTAMINT             123

#define RPC_CREATECUBDIFFCURVE      124

#define RPC_SETCUBDIFFCURVE         125

#define RPC_CREATEZCSWAPCUBDIFF     126

#define RPC_SETZCSWAPCUBDIFF        127

#define RPC_CREATEZCSWAPFUTINT      128
 
#define RPC_SETZCSWAPFUTINT         129

#define RPC_GETINITIALCURVEFROMSUMMIT    130

#define RPC_CREATEZCINTSMOOTH       131
 
#define RPC_SETZCINTSMOOTH          132

#define RPC_CREATEZCSWAPINTSMOOTH   133 

#define RPC_SETZCSWAPINTSMOOTH      134 

#define RPC_CREATEZCSWAPFUTINTSMOOTH 135 

#define RPC_SETZCSWAPFUTINTSMOOTH   136 

#define RPC_CREATE_TRANS2SMOOTH     137

#define RPC_SET_TRANS2SMOOTH        138

#define RPC_CREATETOYZCSWAPINT      139

#define RPC_SETTOYZCSWAPINT         140

#define RPC_CREATE_BUMPCURVE        141

#define RPC_SET_BUMPCURVE           142



/*---- BONDS ----*/

#define RPC_CREATEBOND                  200

#define RPC_SETBOND                     201

#define RPC_CREATEBASISYIELDCURVEMODEL  202

#define RPC_YIELDTOPRICE                203

#define RPC_PRICETOYIELD                204

#define RPC_BONDFORWARDACTUARIELPRICE   205

#define RPC_BONDREPORATEFORFORWARDDATE  206

#define RPC_GET_BOND_FROM_SUMMIT        207

#define RPC_SET_CREATED_BOND_FROM_SUMMIT    208

#define RPC_YIELDTODURATION                 209

#define RPC_YIELDTOCONVEXITY                210

#define RPC_SETYIELD                        211


/*---- MODELS ----*/


#define RPC_CREATEHWTREEMODEL                   300

#define RPC_SETHWTREEMODEL                      301

#define RPC_CREATECRRTREEMODEL                  302

#define RPC_SETCRRTREEMODEL                     303
 
#define RPC_CREATEYCMODEL                       304
 
#define RPC_SETYCMODEL                          305

#define RPC_CREATEBSMODEL                       306
 
#define RPC_SETBSMODEL                          307

#define RPC_CREATEGYCMODEL                      308
 
#define RPC_SETGYCMODEL                         309

#define RPC_CREATEDFGYCMODEL                    310
 
#define RPC_SETDFGYCMODEL                       311
 
#define RPC_CREATEIR3DTHWMODEL                  312
 
#define RPC_SETIR3DTHWMODEL                     313

#define RPC_CREATEHW2FMODEL                     314
 
#define RPC_SETHW2FMODEL                        315
 
#define RPC_GETPARAMETER                        316

#define RPC_CREATEBKIRTREE                      317
 
#define RPC_SETBKIRTREE                         318
 
#define RPC_COVARTXFWD                          319

#define RPC_CORRTXFWD                           320

#define RPC_CREATEHW2FTREEMODEL                 321
 
#define RPC_SETHW2FTREEMODEL                    322

#define RPC_CREATEFNHWMONTECARLO                323

#define RPC_SETFNHWMONTECARLO                   324

#define RPC_CREATEGYCLSMODEL                    325

#define RPC_SETGYCLSMODEL                       326

#define RPC_CREATEG2YCMODEL                     327

#define RPC_SETG2YCMODEL                        328

#define RPC_CREATEBSSLMODEL                     329
 
#define RPC_SETBSSLMODEL                        330

#define RPC_CREATEHWSIGCST                      331

#define RPC_SETHWSIGCST                         332

#define RPC_CREATEHWSIGVAR                      333

#define RPC_SETHWSIGVAR                         334

#define RPC_CREATEHWSIGVARANALYTIC              335
 
#define RPC_SETHWSIGVARANALYTIC                 336

#define RPC_CREATEFNHWMONTECARLOSV              337
 
#define RPC_SETFNHWMONTECARLOSV                 338

#define RPC_CREATEDFBSMODEL                     339

#define RPC_SETDFBSMODEL                        340

#define RPC_CREATEDFFXBSPLAINMODEL              341

#define RPC_SETDFFXBSPLAINMODEL                 342

#define RPC_CREATEDFFXBSMODEL                   343

#define RPC_SETDFFXBSMODEL                      344

#define RPC_CREATEDFHWXSIGVAR                   345
 
#define RPC_SETDFHWXSIGVAR                      346

#define RPC_CREATEHW2FSIGVARANALYTIC            347
 
#define RPC_SETHW2FSIGVARANALYTIC               348

#define RPC_CREATEFNHW2FMONTECARLOSV            349
 
#define RPC_SETFNHW2FMONTECARLOSV               350

#define RPC_CREATEDFHWSIGVARTREE                351

#define RPC_SETDFHWSIGVARTREE                   352

#define RPC_CREATELOGDECANA                     353

#define RPC_SETLOGDECANA                        354

#define RPC_CREATELOGDECMC                      355

#define RPC_SETLOGDECMC                         356

#define RPC_CREATELOGDECANA_P                   357

#define RPC_SETLOGDECANA_P                      358

#define RPC_CREATELOGDECMC_P                    359

#define RPC_SETLOGDECMC_P                       360

#define RPC_CREATELOGDECMC_FROM_ANA             361

#define RPC_SETLOGDECMC_FROM_ANA                362

#define RPC_CREATEFRMANA                        363

#define RPC_SETFRMANA                           364

#define RPC_CREATEFRMANA_PORT                   365

#define RPC_SETFRMANA_PORT                      366

#define RPC_CREATEFRMLSMONTECARLO               367

#define RPC_SETFRMLSMONTECARLO                  368

#define RPC_CREATEFRMTREE                       369

#define RPC_SETFRMTREE                          370

#define RPC_CREATEFRMTREE_AUTO                  371

#define RPC_SETFRMTREE_AUTO                     372

#define RPC_FRMSHORTRATEVOLS                    373

#define RPC_CREATEFRMLSMC_AUTO                  374

#define RPC_SETFRMLSMC_AUTO                     375

#define RPC_CREATEFRMLSMC_AUTO_2                376

#define RPC_SETFRMLSMC_AUTO_2                   377

#define RPC_CREATESMILEDMCRNLDC                 378

#define RPC_SETSMILEDMCRNLDC                    379

#define RPC_CREATEFITTEDSMILEDMCRNLDC           380

#define RPC_SETFITTEDSMILEDMCRNLDC              381

#define RPC_CREATESMILEDLDCANA                  382

#define RPC_SETSMILEDLDCANA                     383

#define RPC_CREATESMILEDLDCFROMANA              384

#define RPC_SETSMILEDLDCFROMANA                 385

#define RPC_CREATEGLOBDFBSMODEL                 386

#define RPC_SETGLOBDFBSMODEL                    387

#define RPC_BMCFRM_CREATE                       388

#define RPC_BMCFRM_SET                          389

#define RPC_BMCFRM2CR_CREATE                    390

#define RPC_BMCFRM2CR_SET                       391

#define RPC_CREATEBSSMILEDMODEL                 392

#define RPC_SETBSSMILEDMODEL                    393

#define RPC_CREATEFRMTREE_AUTO_G                394

#define RPC_SETFRMTREE_AUTO_G                   395

#define RPC_CREATEFRMLSMC_AUTO_G                396

#define RPC_SETFRMLSMC_AUTO_G                   397 

#define RPC_CREATEFRMTREE_AUTO_B                398

#define RPC_SETFRMTREE_AUTO_B                   399

#define RPC_CREATE_BSCORRMODEL                  9000

#define RPC_SET_BSCORRMODEL                     9001

#define RPC_GET_SABR_SIGMA                      9002

#define RPC_CREATE_HWSIGVARFROMANA              9003

#define RPC_SET_HWSIGVARFROMANA                 9004

#define RPC_CREATE_HWSIGVARCALIB                9005

#define RPC_SET_HWSIGVARCALIB                   9006 

#define RPC_CREATE_HWSIGVAR_CALIBRATOR          9007

#define RPC_SET_HWSIGVAR_CALIBRATOR             9008

#define RPC_CREATE_CALIBRATE_BSSMILED           9009

#define RPC_SET_CALIBRATE_BSSMILED              9010

#define RPC_CREATE_GET_SABR_SIG_RHO_NU          9011

#define RPC_SET_GET_SABR_SIG_RHO_NU             9012

#define RPC_CREATE_FRMMARKOVTREE                9013

#define RPC_SET_FRMMARKOVTREE                   9014

#define RPC_CREATE_BOOTSCALIB_FRMMODEL          9015

#define RPC_SET_BOOTSCALIB_FRMMODEL             9016

#define RPC_CREATE_CALIBRATION_FRMMODEL         9017

#define RPC_SET_CALIBRATION_FRMMODEL            9018


/* 
   !!!!!!!! Next time the INDEX for a model 
            should begin for instance at : 9000
*/





/*---- OPTIONS ----*/

#define RPC_CREATEOPTION            400

#define RPC_SETOPTION               401

#define RPC_CREATE_EXOPTION         402

#define RPC_SET_EXOPTION            403

#define RPC_VOLIMP                  404

#define RPC_KIMP                    405

#define RPC_GETUNDPRICE             406

#define RPC_BSOPTION                407 

#define RPC_BSDELTA                 408

#define RPC_BSVEGA                  409

#define RPC_BSTHETA                 410

#define RPC_BSGAMMA                 411

#define RPC_CREATEPORTOPTION        412

#define RPC_SETPORTOPTION           413




/*---- BARRIER ----*/
 
#define RPC_CREATE_CONSTBARRIER     430

#define RPC_SET_CONSTBARRIER        431

#define RPC_CREATE_BARRIER          432

#define RPC_SET_BARRIER             433




/*---- SWAPTIONS ----*/

#define RPC_CREATE_SWAPTION             450

#define RPC_SET_SWAPTION                451

#define RPC_CREATE_LIBORSWAPTION        452

#define RPC_SET_LIBORSWAPTION           453

#define RPC_CREATE_GEN_SWAPTION         454

#define RPC_SET_GEN_SWAPTION            455

#define RPC_CREATE_EXOSWAPTION          456

#define RPC_SET_EXOSWAPTION             457

#define RPC_IMPLIED_VOL                 458

#define RPC_CREATE_CFEXOSWAPTION        459

#define RPC_SET_CFEXOSWAPTION           460

#define RPC_IMPLIED_SPREAD_WITH2MODELS  461

#define RPC_CREATE_VARFIX_SWOPT         462

#define RPC_SET_VARFIX_SWOPT            463

#define RPC_CREATE_OPTIONALACCRUALZC    464

#define RPC_SET_OPTIONALACCRUALZC       465

#define RPC_CREATE_FLEXACCSWAPTION      466

#define RPC_SET_FLEXACCSWAPTION         467


/*---- CCY ----*/
 
#define RPC_CREATECCY                   500
 
#define RPC_SETCCY                      501

#define RPC_ISBUSINESSDAY               502

#define RPC_NEXTBUSINESSDAY             503

#define RPC_CREATEISOCCY                504

#define RPC_SETISOCCY                   505

#define RPC_ADJUSTBUSDATE               506

#define RPC_GETDEFIDXFROMCCY            507

#define RPC_GETPAYCALNAME               508

#define RPC_GETCCYNAME                  509

#define RPC_ADDMONTHS                   510

#define RPC_ADDYEARS                    511

#define RPC_GETSPOTDAYSFROMCCY          512




/*---- FOREX ----*/
 
#define RPC_CREATEFOREX                 600
 
#define RPC_SETFOREX                    601
 



/*---- PIBOR FUTURE ----*/
 
#define RPC_CREATEPIBFUT                700

#define RPC_SETPIBFUT                   701



/*---- IRINDEX FUTURE ----*/
 
#define RPC_CREATEIRFUT                 710

#define RPC_SETIRFUT                    711

#define RPC_CREATE_GEN_IRFUT            712

#define RPC_SET_GEN_IRFUT               713
 


/*---- BOND FUTURE ----*/

#define RPC_CREATEBDFUT                 803

#define RPC_SETBDFUT                    804
 
#define RPC_CREATEBDBASKFUT             805
 
#define RPC_SETBDBASKFUT                806

#define RPC_CREATE_NOTIONNAL_GILT_BUND  807

#define RPC_SET_NOTIONNAL_GILT_BUND     808

#define RPC_GETCHEAPEST                 809

#define RPC_GETCONVERSIONFACTOR         810



/*---- PORTFOLIO ----*/

#define RPC_CREATEPF                    901

#define RPC_SETPF                       902

#define RPC_MKTPRICE                    903

#define RPC_PFTHEOPRICE                 904

#define RPC_PFVSKESTIMATE               905

#define RPC_SETPFVSKESTIMATED           906

#define RPC_PFSPLESTIMATE               907

#define RPC_SETPFSPLESTIMATED           908

#define RPC_PFHWXESTIMATE               909
 
#define RPC_SETPFHWXESTIMATED           910
 
#define RPC_PFMODESTIMATE               911
 
#define RPC_SETPFMODESTIMATED           912

#define RPC_PFHWSIGVARESTIMATE          913

#define RPC_SETPFHWSIGVARESTIMATED      914

#define RPC_PFHWSIGVARGLOBMIN           915

#define RPC_SETPFHWSIGVARGLOBMIN        916

#define RPC_PFHWSIGVARPENALFIT          917

#define RPC_SETPFHWSIGVARPENALFIT       918

#define RPC_PFHW2FPENALFIT              919

#define RPC_SETPFHW2FPENALFIT           920

#define RPC_PFHW2FSIGVARESTIMATE        921

#define RPC_SETPFHW2FSIGVARESTIMATED    922

#define RPC_PFHW2FSIGVARGLOBMIN         923

#define RPC_SETPFHW2FSIGVARGLOBMIN      924

#define RPC_PFLOGDECVOLFIT              925

#define RPC_SETPFLOGDECVOLFIT           926

#define RPC_PFINSTLOGDECVOLFIT          927

#define RPC_SETPFINSTLOGDECVOLFIT       928

#define RPC_INSTLOGDECVOLFIT            929

#define RPC_SETINSTLOGDECVOLFIT         930




/*---- STRUCTURE ----*/
 
#define RPC_CREATESTRUCTURE             951
 
#define RPC_SETSTRUCTURE                952




/*----  IASECURITY --- */

#define RPC_CREATE_IASEC                970

#define RPC_SET_IASEC                   971


/*----  REVERSE & REVERSECOUPON --- */

#define RPC_CREATE_REVERSE              980

#define RPC_SET_REVERSE                 981

#define RPC_CREATE_REVERSECOUPON        982

#define RPC_SET_REVERSECOUPON           983
 
#define RPC_CREATE_STRUCTREVERSECOUPON  984
 
#define RPC_SET_STRUCTREVERSECOUPON     985
 
#define RPC_CREATE_STRUCTREVERSECOUPON2 986
 
#define RPC_SET_STRUCTREVERSECOUPON2    987
 
#define RPC_CREATE_REVERSE_CALENDAR     988
 
#define RPC_SET_REVERSE_CALENDAR        989

 
 
 
 
/*---- IR INDEX  ----*/

#define RPC_CREATE_IRINDEX              1000

#define RPC_SET_IRINDEX                 1001

#define RPC_CREATE_LIBOR                1002

#define RPC_SET_LIBOR                   1003

#define RPC_CREATE_IRINDEX_MONEY_MARKET 1004

#define RPC_SET_IRINDEX_MONEY_MARKET    1005

#define RPC_CREATE_CMS                  1006

#define RPC_SET_CMS                     1007

#define RPC_CREATE_FIX                  1008

#define RPC_SET_FIX                     1009



/*----   VOL    -----*/

#define RPC_CREATEVOLCURVELIN           1100

#define RPC_SETVOLCURVELIN              1101

#define RPC_CREATEVOLFLAT               1102

#define RPC_SETVOLFLAT                  1103

#define RPC_GETVOLFROMSUMMIT            1104

#define RPC_SETVOLFROMSUMMIT            1105

#define RPC_COMPUTEVOLATILITY           1106

#define RPC_GETVOLCUBEFROMSUMMIT        1107

#define RPC_SETVOLCUBEFROMSUMMIT        1108

#define RPC_CREATEVOLCUBE               1109

#define RPC_SETVOLCUBE                  1110

#define RPC_GETFXVOLFROMSUMMIT          1111

#define RPC_SETFXVOLFROMSUMMIT          1112

#define RPC_CREATEBUMPVOLATILITY        1113

#define RPC_SETBUMPVOLATILITY           1114



/*---- SWAP LEGS ----*/

#define RPC_CREATE_SWAPLEG              2000

#define RPC_SET_SWAPLEG                 2001

#define RPC_CREATE_FIXEDLEG             2003

#define RPC_SET_FIXEDLEG                2004

#define RPC_CREATE_LIBORLEG             2005

#define RPC_SET_LIBORLEG                2006

#define RPC_CREATE_TMLEG                2007

#define RPC_SET_TMLEG                   2008

#define RPC_CREATE_CMSLEG               2009
 
#define RPC_SET_CMSLEG                  2010

#define RPC_CREATE_CMTLEG               2011

#define RPC_SET_CMTLEG                  2012

#define RPC_CREATESET_FIXEDRATES        2013

#define RPC_SETSET_FIXEDRATES           2014

#define RPC_CREATE_SWAP_WITH_NOTIONNAL  2015

#define RPC_SET_SWAP_WITH_NOTIONNAL     2016

#define RPC_GEN_AMORT                   2017

#define RPC_SET_GEN_AMORT               2018

#define RPC_CUST_FIRST_PERIOD           2019

#define RPC_DISP_SCHED_DATES            2020

#define RPC_DISP_SCHED_VALUES           2021




/*---- SWAPS ----*/

#define RPC_CREATE_SWAP                 3000

#define RPC_SET_SWAP                    3001

#define RPC_CREATE_LIBORSWAP            3002

#define RPC_SET_LIBORSWAP               3003

#define RPC_SWAP_PRICE_TO_RATE          3004

#define RPC_SWAP_RATE_TO_PRICE          3005

#define RPC_GET_SWAP_FROM_SUMMIT        3006

#define RPC_SET_CREATED_SWAP_FROM_SUMMIT    3007

#define RPC_CREATE_GEN_SWAP             3008

#define RPC_SET_GEN_SWAP                3009

#define RPC_IMPLIED_SPREAD              3010

#define RPC_LIBOR_ASSET_SWAP_MARGIN     3011



/*---- CAP & FLOOR ----*/

#define RPC_CREATE_CAPFLOOR           4000

#define RPC_SET_CAPFLOOR              4001

#define RPC_CREATE_LIBORCAPFLOOR      4002

#define RPC_SET_LIBORCAPFLOOR         4003

#define RPC_CREATE_FLEXIBLECAPFLOOR   4004

#define RPC_SET_FLEXIBLECAPFLOOR      4005

#define RPC_CREATE_LIBORFLXCAPFLOOR   4006

#define RPC_SET_LIBORFLXCAPFLOOR      4007

#define RPC_CREATE_GEN_CAPFLOOR       4008

#define RPC_SET_GEN_CAPFLOOR          4009

#define RPC_CAPLETPRICE               4010

#define RPC_CREATE_EXOFLEXCF          4011

#define RPC_SET_EXOFLEXCF             4012

#define RPC_CREATE_RATCHET            4013
 
#define RPC_SET_RATCHET               4014

#define RPC_CREATE_STICKY             4015

#define RPC_SET_STICKY                4016

#define RPC_CREATE_SPREADLEG          4017

#define RPC_SET_SPREADLEG             4018

#define RPC_CREATE_CORRIDORLEG        4019
 
#define RPC_SET_CORRIDORLEG           4020

#define RPC_CREATE_RESTRIKABLELEG     4021

#define RPC_SET_RESTRIKABLELEG        4022

#define RPC_IMPLIED_RANGE             4023

#define RPC_CREATE_MATCAPFLOOR        4024

#define RPC_SET_MATCAPFLOOR           4025 

#define RPC_CREATE_SPREADOPTION       4026

#define RPC_SET_SPREADOPTION          4027

#define RPC_CREATE_DIGITAL            4028

#define RPC_SET_DIGITAL               4029



/*---- RANGE NOTE ----*/

#define RPC_CREATE_RNGNOTE            5000

#define RPC_SET_RNGNOTE               5001

#define RPC_CREATE_LIBORRNGNOTE       5002

#define RPC_SET_LIBORRNGNOTE          5003




/*---- EXERCISE STYLE ----*/

#define RPC_CREATE_EUROPXSTYLE        6000
 
#define RPC_SET_EUROPXSTYLE           6001

#define RPC_CREATE_AMERICANXSTYLE     6002
 
#define RPC_SET_AMERICANXSTYLE        6003

#define RPC_CREATE_BERMUDANXSTYLE     6004
 
#define RPC_SET_BERMUDANXSTYLE        6005

#define RPC_CREATE_CUSTOMXSTYLE       6006
 
#define RPC_SET_CUSTOMXSTYLE          6007




/*---- REFERENCE VALUE ----*/

#define RPC_CREATE_CONSTREFVALUE      6100

#define RPC_SET_CONSTREFVALUE         6101

#define RPC_CREATE_REFVALUE           6102
 
#define RPC_SET_REFVALUE              6103

#define RPC_CREATE_IA3LEVREFVAL       6104

#define RPC_SET_IA3LEVREFVAL          6105




/*----   UTIL  ------ */

#define RPC_GETPRTYBYNAME            7000



/*----   DB    ------ */

#define RPC_DBGETSWAPTION            7500

#define RPC_DBSTORESWAPTION          7501

#define RPC_DBGETBERMUDANSWAPTION    7502

#define RPC_DBSTOREBERMUDANSWAPTION  7503


/*---- CRF ---*/


#define RPC_CREATE_REVERSEFLOAT      8001
#define RPC_SET_REVERSEFLOAT         8002


/*----------------------------------------------------------------------------*/
/*---- end of file  ----*/
