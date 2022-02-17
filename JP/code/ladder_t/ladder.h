/****************************************************************************/
/*      Declaration file for structure of deal data.                        */
/****************************************************************************/
/*      LADDER.h                                                            */
/****************************************************************************/


typedef struct                        /* Deal input data structure        */
{                  
    char    
                SoZ,                  /**< Swap or zero coupon              */
                LoS,                  /**< JPM long or short the option     */
                Style,                /**< Option style                     */
                RibSmoothing,         /**< Smoothing of range profile       */
                RibIdxFreq[2],        /**< Rib index frequency              */
                RibIdxDCC[2],         /**< Rib index day count              */
                Smoothing,            /**< Smoothing flag                   */
                PayFreqSt,            /**< Payment frequency of sticky leg  */
                PayFreqFl,            /**< Payment frequency of floating leg*/
                DayCountSt,           /**< Day count for sticky payment     */
                DayCountFl,           /**< Day count for floating payment   */
                CompSt,               /**< Simple or compound sticky cpn    */
                CompFl,               /**< Simple or compound floating cpn  */
                ArrearsSt,            /**< Arrears (Y or N) of sticky leg   */
                ArrearsFl,            /**< Arrears (Y or N) of floating leg */
                StubConv,             /**< Stub convention (F or B)         */
                IdxBaseFl,            /**< Day count of floating index      */
                IdxFreqFl,            /**< Frequency of floating index      */
                IdxNameFl[200],       /**< Name of floating index           */
                IdxOnFl,              /**< TRUE if used False if not        */
                IdxBaseSt[2],         /**< Day count of sticky indices      */
                IdxFreqSt[2],         /**< Frequency of sticky indices      */
                IdxOnSt[2],           /**< TRUE if used False if not        */
                IdxNameSt[2][200],    /**< Name of sticky indices           */
                CalcStats,            /**<'Y' - yes, 'N' - no stat calculation */
                OptStatStyle,         /**< Exercise stats style             */
                AoM;                  /**< Additive or multiplicative       */

    int 
                NbAmort,              /**< Number of amortization dates     */
                NbExer,               /**< Number of exercise dates         */
                NbStepUpFl,           /**< Number of step up dates          */
                NbStepUpSt,           /**< Number of sticky dates           */
                NbFixing,             /**< Number of prev fixing dates      */
                FixingGivenFl,        /**< TRUE if user has spec float fix  */
                FixingGivenSt,        /**< TRUE if user has spec sticky fix */
                CplxIsRib,            /**< TRUE if complex leg is Rib       */
                NbRibObsDates,        /**< Number of rib observation dates  */
                NbPastRibObsDates,    /**< Number of past rib obs dates     */
                RibIdxMat[2],         /**< Rib index maturity in months     */
                RibIdxIoD[2],         /**< Rib index curve                  */
                RibIdxOn[2],          /**< TRUE if rib obs used False if not*/
                NbRibObsInPer[MAXNBDATE],/**< Nb Rib obs in ladder period   */
                MaxNbRib,             /**< Max nb Rib obs in ladder period  */
                IdxMatFl,             /**< Maturity of floating index       */
                FirstResetI,          /**< First reset date idx for sticky  */
                EndResetI,            /**< Last reset date idx for sticky   */
                IdxMatSt[2],          /**< Maturity of sticky index         */
                IdxIoDFl,             /**< Curve for floating index         */
                IdxIoDSt[2],          /**< Curve for sticky index           */
                NbStates,             /**< Number of state variables        */
                NbStDev;              /**< 'width' of state var estimiated  */
    long 
                AmortDate[MAXNBDATE], /**< Amortization dates               */
                FlDate[MAXNBDATE],    /**< Step up date                     */
                StDate[MAXNBDATE],    /**< St date                          */
                ExerDate[MAXNBDATE],  /**< Exercise dates                   */
                FixingDate[MAXNBDATE],/**< Past fixings of sticky idx dates */
                OptStatDates[MAXNBDATE], /**< Exercise stats dates          */
                *RibObsDate,          /**< Rib obs dates                    */
                *RibObsEffDate,       /**< Rib obs eff dates                */
                FirstResetDateSt,     /**< First reset date for sticky      */
                FirstResetDateFl,     /**< First reset date for funding     */
                AccStDate,            /**< Swap accrual start date          */
                OptNbStats,           /**< Number of Exercise stats         */
                MatDate;              /**< Swap final maturity date         */
    double                                        
                Notional,             /**< Absolute value of notional       */
                NotionalSign,         /**< Absolute value of notional       */
                OrgNotPerc,           /**< Original notional percentage     */
                Amort[MAXNBDATE],     /**< Amortization array               */
                *RibHiBarrier,        /**< Range barrier level              */
                *RibLoBarrier,        /**< Range barrier level              */
                *RibInRangeWeight,    /**< In-range obsevent weight         */
                *RibOutRangeWeight,   /**< Out-of-range obsevent weight     */
                RibIdxWeight[2],      /**< Rib index weights                */
                RibPastObsPerc,          /**< Rib past obs in range %       */
                RibPastObsWeight,        /**< Past Rib obs / curr period    */
                RibPastObsInRangeWeight, /**< Past Rib obs, in range weight */
                RibPastObsOutRangeWeight,/**< Past Rib obs, out range weight*/
                StickyCoef[MAXNBDATE],/**< Sticky coefficient               */
                CapSt[MAXNBDATE],     /**< Cap on sticky cpn                */
                FloorSt[MAXNBDATE],   /**< Floor on sticky cpn              */
                UpRateSt[MAXNBDATE],  /**< Up step in sticky cpn            */
                MidRateSt[MAXNBDATE], /**< Mid step in sticky cpn           */
                DownRateSt[MAXNBDATE],/**< Down step in sticky cpn          */
                IdxWeightSt[2][MAXNBDATE],/**< Weight of sticky index       */
                IdxObsWeightSt[2],    /** Observation weight of sticky index*/
                SprdSt[MAXNBDATE],    /**< Spread of sticky coupon          */
                CapFl[MAXNBDATE],     /**< Cap on funding cpn               */
                FloorFl[MAXNBDATE],   /**< Floor on funding cpn             */
                UpRateFl[MAXNBDATE],  /**< Up rate in funding cpn           */
                BarrierLo[MAXNBDATE], /**< Low barrier                      */
                BarrierHi[MAXNBDATE], /**< High barrier                     */
                Leverage[MAXNBDATE],  /*   Leverage                         */
                IdxFloor[MAXNBDATE],  /*   Floor on the idx in ladder leg   */
                IdxCap[MAXNBDATE],    /*   Cap on the idx in the ladder leg */
                Strike[MAXNBDATE],    /**< Strikes                          */
                FixingFl,             /**< 1st fixing of floating index     */
                FixingSt[2][MAXNBDATE],/**< Past fixings of sticky indices  */
                FixingRibPerc[MAXNBDATE],/**< Past Rib % in range           */
                IdxWeightFl,          /**< Weight of flaoting index         */
                FirstLevel,           /**< lkbk for the very 1st reset      */
                InitOuts,             /**< Cumulative notional up value date*/
                InitState,            /**< Initial state variable           */
                LastMinState,         /**< Last state variable              */
                LastMaxState,         /**< Last state variable              */
                LastFwdState;         /**< Last state variable              */

    double      **State;              /**< State variable slices            */

} LADDER_DATA;
