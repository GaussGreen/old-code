/****************************************************************************/
/*      Declaration file for structure of deal data.                        */
/****************************************************************************/
/*      fxkoseries_t.h                                                      */
/****************************************************************************/

/*
$Header$
*/



typedef struct                        /* Deal input data structure          */
{                                     /*                                    */
    char                              /*                                    */
        LoS,                          /* Long or short the option (L or S)  */
        CoP,                          /* Call or put option ('C' or 'P')    */
        KnockIoO,                     /* Knock out or in                    */
        IoO,                          /* Inside or outside range            */
        ExerFreq,                     /* E=european, N=american, M, Q, S, A */
        KoFreq,                       /*                                    */
        SmoothFlag;

    int                               /*                                    */
        NbExer,                       /* Number of exercise dates           */
        NbKoDates;                    /* Number of KO dates                 */
   
    long                              /*                                    */
        Exer[MAXNBDATE],              /* Exercise dates                     */
        KoDates[MAXNBDATE];           /* KO dates                           */
   
    double                            /*                                    */
        Notional,                     /* Notional in foreign currency       */
        Strike[MAXNBDATE],            /* Strikes                            */
        LoBarrier[MAXNBDATE],
        HiBarrier[MAXNBDATE],
        Rebate[MAXNBDATE];


} FXKOSERIES_DATA;
