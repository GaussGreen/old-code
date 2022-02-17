
#ifndef GRF_H_MDLCOMM
#define GRF_H_MDLCOMM

/*
        grf_f_mdlcomm.h
        ===============

        Author: Toto 26Nov1999

        Objective: define communications between any model and GRFN
*/

/*
        Model communicates with GRFN via the local FIRSTMktAtT structure      ,
        containing all the relevant market information for evaluation
        of GRFN at time T and the global FIRSTAllMkts structure      ,
        that contains an array of FIRSTMktAtT
*/

#define ONE_MINUTE 1.90259E-06

/*
        Structures
        ==========
*/

typedef struct {

  /*	Old SORT representation of market information	*/
  GrfnEvent *evt;
  SrtSample smp;

} FIRSTMktAtT;

typedef struct {

  /*	Old SORT representation	*/
  GrfnDeal *gd;
  SrtUndInfo *und_info;

  /*	Underlying data	*/
  int num_und;
  Und_Data *und_data;
  SrtUndPtr und_ptr[MAXUNDERLYING];

  /*	Number of events in table	*/
  int num_evt;
  /*	Number of columns	*/
  int num_cols;

  /*	Event dates	*/
  long *dts;
  /*	Event times (in years)	*/
  double *tms;
  /*	Array of the events	*/
  FIRSTMktAtT *evt;
  /*	1 if events are AM      , 0 otherwise	*/
  int *am;

} FIRSTAllMkts, *GRFNCOMMSTRUCT;

/*
        Init function
        =============

        Inputs:	GRFN tableau
                        Domestic underlying name
                        GRFN param

        Output:	A pointer on a FIRSTAllMkts structure      ,
                        fully allocated      ,
                        with requests on market parameters (df maturities)

        Initialise and process GRFN tableau and value historical events
*/

Err FIRSTInitMktStruct(
    /*	GRFN Tableau	*/
    int init_num_evt_dts, long *init_evt_dts, long tab_rw, long tab_cl,
    char ***tab_str, int **mask, long aux_width, long *aux_len, double **aux,
    /*	Domestic underlying name	*/
    char *dom_nme,
    /*	GRFN param	*/
    SrtGrfnParam *prm,
    /*	GRFN choice: backward (-1) or forward (+1)	*/
    int *back_or_for,
    /*	Output	*/
    GRFNCOMMSTRUCT comm);

/*
        Functions used to get the global deal information
        =================================================
*/

/*	1.	Get underlyings to be used to value the deal	*/

Err FIRSTGetUndFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Number of underlyings needed to value the deal	*/
    int *num_und,
    /*	Pointer is allocated inside      ,
            and must be freed by call to FIRSTFreeUndFromDeal	*/
    SrtUndPtr **und_ptr);

/*	2.	Free the previous function result	*/

Err FIRSTFreeUndFromDeal(int num_und, SrtUndPtr **und_ptr);

/*	3.	Get number of columns	*/

Err FIRSTGetNumColFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Information	*/
    int *num_col);

/*	4.	Get pv of past (to be added to last column)	*/

Err FIRSTGetPvOfPastFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	PV to be added to last column before simulation (MC) or after (Tree)
     */
    double *pv_of_past);

/*	5.	Get the maximum number of dfs required	*/

Err FIRSTGetMaxNumDfFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Max num df	*/
    int *max_num_df);

/*	6.	Get event dates	*/

Err FIRSTGetEvtDatesFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Number of event dates	*/
    int *num_evt,

    /*	Pointers are allocated inside      ,
            and must be freed by call to FIRSTFreeEvtDatesFromDeal	*/

    /*	Event dates and times	*/
    long **evt_dts, double **evt_tms);

/*	7.	Free the previous function result	*/

Err FIRSTFreeEvtDatesFromDeal(int num_evt, long **evt_dts, double **evt_tms);

/*	8.	Get events and df maturities required to calculate them	*/

Err FIRSTGetEventInfoFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	As output from FIRSTGetEvtDatesFromDeal	*/
    int num_evt,
    /*	As output from FIRSTGetUndFromDeal	*/
    int num_und,
    /*	As output from FIRSTGetMaxNumDfFromDeal	*/
    int max_num_df,

    /*	All pointers are allocated inside      ,
            and must be freed by call to FIRSTFreeEventInfoFromDeal	*/

    /*	Events themselves	*/
    FIRSTMktAtT **evt,

    /*	Wether American	*/
    int **am,

    /*	Information relative to df required for event evaluation	*/

    /*	num_df_mat[i][j] = number of df required for event i from underlying j
     */
    int ***num_df_mat,
    /*	df_mat_dts[i][j][k] and df_mat_tms[i][j][k]
            = maturity of the df number k required from underlying j at event i
     */
    long ****df_mat_dts, double ****df_mat_tms);

/*	9.	Free the previous function result	*/

Err FIRSTFreeEventInfoFromDeal(int num_evt, int num_und, int max_num_df,
                               FIRSTMktAtT **evt, int **am, int ***num_df_mat,
                               long ****df_mat_dts, double ****df_mat_tms);

/*
        Functions used to get the local deal information
        ================================================
*/

/*	1.	Get event and df maturities required at one given time	*/

Err FIRSTGetLocalEventInfoFromDeal(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Current time	*/
    double time,
    /*	Return NULL if no event at this time      , or pointer on the event
     */
    FIRSTMktAtT **evt,

    /*	POINTERS num_df_mat      , df_mat_dts and df_mat_tms
            ARE NOT ALLOCATED INSIDE
            they must be allocated prior to function call	*/

    /*	num_df_mat[i] = number of df required from underlying i	*/
    int *num_df_mat,
    /*	df_mat_dts[i][j] and df_mat_tms[i][j]
            = maturity of the df number j required from underlying i	*/
    long **df_mat_dts, double **df_mat_tms);

/*	2.	Set DF values prior to valuation	*/

Err FIRSTSetDFValue(
    /*	Pointer on the event to be valued      ,
            must be a valid non-empty event	*/
    FIRSTMktAtT *evt,
    /*	Index of the underlying	*/
    int und_idx,
    /*	Index of the df	*/
    int df_idx,
    /*	Value of df	*/
    double df_val);

/*	3.	Set state variable values prior to valuation	*/

Err FIRSTSetSVValue(
    /*	Pointer on the event to be valued      ,
            must be a valid non-empty event	*/
    FIRSTMktAtT *evt,
    /*	Index of the underlying	*/
    int und_idx,
    /*	Type of the statevar (SPOT      , R      , PHI      , ...)	*/
    int sv_type,
    /*	Value of statevar	*/
    double sv_val);

/*
        Functions used to eval event
        ============================

        DFs and stevars are supposed to have been set accordingly by caller
        before calling this function
*/

Err FIRSTEvalEvent(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Pointer on the event to be valued      ,
            must be a valid non-empty event	*/
    FIRSTMktAtT *evt,
    /*	Number of columns	*/
    int num_col,
    /*	Information needed by Forward/Backward models only
            Please ask Eric Fournie what these are about...
            If you don't want to use these      , please set them to default
     */
    int type_eval, /*	Default is 2	*/
    double *fwd,   /*	Default is NULL	*/
    double *cur,   /*	Default is NULL	*/
    /*	Vector of PVs so far of columns      , will be updated by function
            tree only	*/
    double *col_pvs,
    /*	Cash flow
            MC only	*/
    double *cash_flow);

// IMPLEMENTATION OF PVR
Err FIRSTEvalEventCredit(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm,
    /*	Pointer on the event to be valued      ,
            must be a valid non-empty event	*/
    FIRSTMktAtT *evt,
    /*	Number of columns	*/
    int num_col,
    /*	Information needed by Forward/Backward models only
            Please ask Eric Fournie what these are about...
            If you don't want to use these      , please set them to default
     */
    int type_eval, /*	Default is 2	*/
    double *fwd,   /*	Default is NULL	*/
    double *cur,   /*	Default is NULL	*/
    /*	Vector of PVs so far of columns      , will be updated by function
            tree only	*/
    double *col_pvNr, double *col_pvR, double *col_pvSource, int *sourceStatus,
    /*	Cash flow
            MC only	*/
    double *cash_flow);

#if 0
Err FIRSTEvalEvent2(
/*	As output from FIRSTInitMktStruct	*/
GRFNCOMMSTRUCT			comm      ,
/*	Pointer on the event to be valued      , 
	must be a valid non-empty event	*/
FIRSTMktAtT*			evt      ,
/*	Number of columns	*/
int						num_col      ,
/*	Information needed by Forward/Backward models only
	Please ask Eric Fournie what these are about...
	If you don't want to use these      , please set them to default	*/
int						type_eval      ,	/*	Default is 2	*/
double*					fwd      ,		/*	Default is NULL	*/
double*					cur      ,		/*	Default is NULL	*/
/*	Vector of PVs so far of columns      , will be updated by function
	tree only	*/
double*					col_pvs      ,
/*	Cash flow
	MC only	*/
double*					cash_flow);

#endif

/*
        Functions used to free everything after valuation (or on error)
        ===============================================================
*/

Err FIRSTFreeMktStruct(
    /*	As output from FIRSTInitMktStruct	*/
    GRFNCOMMSTRUCT comm);

/*
        A simple example
        ================
*/

/*	Specific parameters of model	*/

typedef struct {
  void *nothing;
} * DetIrUndSpec;

/*	1.	Function to free underlying	*/

Err SrtFreeDetermUnd(void *undPtr);

/*	2.	Function to initialise deterministic model	*/

Err SrtInitDetermUnd(char *undName, char *ycName);

/*	3.	Function to value GRFN deal in deterministic model	*/

Err FirstValueDealInDeterministicIrModel(
    /*	GRFN Tableau	*/
    int num_evt_dts, long *evt_dts, long tab_rw, long tab_cl, char ***tab_str,
    int **mask, long aux_width, long *aux_len, double **aux,
    /*	Domestic underlying name	*/
    char *dom_nme,
    /*	GRFN param	*/
    SrtGrfnParam *prm,
    /*	Premium	*/
    double *price);

#endif