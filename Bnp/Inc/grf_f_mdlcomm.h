
#ifndef GRF_H_MDLCOMM
#define GRF_H_MDLCOMM

#ifdef __NOTDEF

/*
        grf_f_mdlcomm.h
        ===============

        Author: toto 26Nov1999

        Objective: define communications between any model and GRFN
*/

/*	Model communicates with GRFN via the FIRSTMktAtT structure      ,
        containing all the relevant market information for evaluation
        of GRFN at time T	*/

typedef struct {
  /*	For internal use only      , please do not touch	*/
  GrfnEvent *evt;
  SrtSample smp;

  /*	Information accessible to the user      , read-only	*/
  /*	(completed on call to function INIT)	*/

  /*	Date	*/
  double tme; /*	In yeapasswrrs from today	*/
  long dte;   /*	As a date	*/

  /*	Info about dfs	*/

  /*	num_df[i] = number of dfs required from curve number i	*/
  int num_df[MAXUNDERLYING];
  /*	df_mat[i][j] = maturity (in years from today) of the df number j
          from curve number i	*/
  double *df_mat[MAXUNDERLYING];
  /*	Same      , only dates	*/
  double *df_mat_dte[MAXUNDERLYING];

  /*	Information to be provided by the user before evaluation      ,
          once per node - space is allocated on call to function INIT	*/

  /*	Domestic numeraire	*/
  double dom_num;

  /*	Equity spots and associated numeraires	*/
  double spt[MAXUNDERLYING];
  double spt_num[MAXUNDERLYING];

  /*	Fx and associated numeraires	*/
  double fx[MAXUNDERLYING];
  double fx_num[MAXUNDERLYING];

  /*	Other state variables	*/
  double statevar[MAXUNDERLYING][MAXSTATEVAR];

  /*	DFs	*/
  /*	df[i][j] = df number j from curve number i */
  double *df[MAXUNDERLYING];
} FIRSTMktAtT;

typedef struct {
  /*	For internal use only      , please do not touch	*/
  GrfnDeal *gd;
  SrtUndInfo *und_info;

  /*	Underlyings	*/
  int num_und;
  SrtUndPtr und_ptr[MAXUNDERLYING];

  /*	Number of events in table	*/
  int num_evts;
  /*	Array of the events	*/
  FIRSTMktAtT *evts;
} FIRSTAllMkts, *GRFNCOMMSTRUCT;

#endif
#endif