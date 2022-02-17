#ifndef NAG_TYPES
#define NAG_TYPES  

/* <nag_types.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Contains the NAG defined types.
 *
 * Shah Datardina, NAG Ltd., Oxford, U.K., 1996.
 *
 * Mark 4 re-issue, 1996.
 * Mark 5 revised. IER-2138 (Feb 1998).
 * New types for Mark 5 & changes relating to e04mfc/e04nfc 
 * revisions. 1997
 */

#define MatrixTranspose_start 1
#define MatrixTriangle_start 4
#define MatrixUnitTriangular_start 6
#define OperationSide_start 8
#define PivotType_start 10
#define SequenceDirection_start 14
#define NormType_start 16
#define MatrixType_start 20
#define Nag_VectorOp_start 27
#define Nag_TransformDirection_start 29
#define Nag_QuadWeight_start 31
#define Nag_BoundInterval_start 35
#define Nag_TrigTransform_start 38
#define Nag_GaussFormulae_start 40
#define Nag_MCMethod_start 44
#define Nag_ErrorControl_start 46
#define Nag_MeshSet_start 49
#define Nag_RK_task_start 51
#define Nag_RK_method_start 53
#define Nag_ErrorAssess_start 56
#define Nag_SolDeriv_start 58
#define Nag_2d_Scat_Method_start 61
#define Nag_DerivType_start 63
#define Nag_Start_start 65
#define Nag_PrintType_start 72
#define Nag_GradChk_start 84
#define Nag_DPrintType_start 94
#define Nag_CheckType_start 100
#define Nag_DWantType_start 103
#define Nag_DerivInfo_start 107
#define Nag_DerivSet_start 113
#define Nag_FunType_start 118
#define Nag_LinFun_start 120
#define Nag_InitType_start 123
#define Nag_BoundType_start 128
#define Nag_ProblemType_start 134
#define Nag_EndState_start 153
#define Nag_DegenJob_start 186
#define Nag_SetStateMode_start 189
#define Nag_SparseQP_LU_Type_start 191
#define Nag_SparseQP_LU_SolveType_start 194
#define Nag_SparseQP_BasisFactType_start 197
#define Nag_CrashType_start 199
#define Nag_ScaleType_start 203
#define Nag_WhereElements_start 207
#define Nag_InitRotation_start 209
#define Nag_Select_Eigenvalues_start 211
#define Nag_SolveSystem_start 213
#define Nag_SparseNsym_Method_start 219
#define Nag_SparseNsym_Fact_start 222
#define Nag_SparseNsym_Piv_start 224
#define Nag_SparseNsym_Zeros_start 228
#define Nag_SparseNsym_Dups_start 231
#define Nag_SparseSym_Method_start 234
#define Nag_SparseSym_Fact_start 236
#define Nag_SparseNsym_PrecType_start 238
#define Nag_SparseSym_Weight_start 241
#define Nag_SparseSym_Bisection_start 243
#define Nag_SparseSym_Norm_start 245
#define Nag_SparseSym_Term_start 248
#define Nag_SparseSym_CheckData_start 250
#define Nag_SparseSym_Piv_start 252
#define Nag_SparseSym_PrecType_start 255
#define Nag_SparseSym_Dups_start 259
#define Nag_SparseSym_Zeros_start 261
#define Nag_TailProbability_start 263
#define Nag_Scores_start 268
#define Nag_Ties_start 274
#define Nag_IncludeWeight_start 279
#define Nag_IncludeMean_start 281
#define Nag_UpdateObserv_start 283
#define Nag_SumSquare_start 285
#define Nag_Initialize_start 287
#define Nag_Link_start 289
#define Nag_RegType_start 297
#define Nag_PsiFun_start 301
#define Nag_SigmaEst_start 307
#define Nag_CovMatrixEst_start 311
#define Nag_SigmaSimulEst_start 314
#define Nag_PrinCompMat_start 316
#define Nag_PrinCompScores_start 320
#define Nag_Weightstype_start 324
#define Nag_RotationLoading_start 327
#define Nag_TransNorm_start 329
#define Nag_RotationScale_start 335
#define Nag_FacMat_start 337
#define Nag_FacScoreMethod_start 340
#define Nag_FacRotation_start 342
#define Nag_GroupCovars_start 344
#define Nag_MahalDist_start 346
#define Nag_DiscrimMethod_start 348
#define Nag_PriorProbability_start 350
#define Nag_MatUpdate_start 353
#define Nag_DistanceType_start 355
#define Nag_VarScaleType_start 358
#define Nag_ClusterMethod_start 362
#define Nag_DendOrient_start 368
#define Nag_Eigenvalues_start 372
#define Nag_ScaleCriterion_start 374
#define Nag_Blocks_start 376
#define Nag_DiscreteDistrib_start 379
#define Nag_PopVar_start 381
#define Nag_Smooth_Type_start 383
#define Nag_FreqTime_start 385
#define Nag_Likelihood_start 387
#define Nag_LagWindow_start 391
#define NagMeanOrTrend_start 395
#define Nag_LoggedSpectra_start 398
#define Nag_state_start 400
#define Nag_ab_input_start 402
#define Nag_ObserverForm_start 404
#define Nag_ControllerForm_start 406
#define Nag_OutputType_start 408
#define Nag_BB_Fail_start 412
#define Nag_MIP_ProbType_start 415
#define Nag_Node_Selection_start 421
#define Nag_Var_Selection_start 427
#define Nag_Branch_Direction_start 431
#define Nag_NodeStatus_start 435
#define Nag_SortOrder_start 442
#define Nag_SearchMatch_start 444
#define Nag_HashError_start 446


#ifdef NAG_INT_BOOL
typedef int Boolean;
#else
#ifdef NAG_CHAR_BOOL
typedef char Boolean;
#else
typedef signed char Boolean;
#endif
#endif

typedef long Integer;

#ifdef NAG_VOID_STAR
typedef void NAG_HUGE * Pointer;          /* Generic pointer */
#else
typedef char NAG_HUGE * Pointer;
#endif

typedef struct {                /* NAG Complex Data Type */
  double re,im;
} Complex;

/*
 *  Enum types used by the Library
 */

/* Following are used by f06 functions*/
typedef enum { NoTranspose=MatrixTranspose_start, Transpose, ConjugateTranspose } MatrixTranspose;
typedef enum { UpperTriangle=MatrixTriangle_start, LowerTriangle } MatrixTriangle;
typedef enum { UnitTriangular=MatrixUnitTriangular_start, NotUnitTriangular } MatrixUnitTriangular;
typedef enum { LeftSide=OperationSide_start, RightSide } OperationSide;
typedef enum { BottomPivot=PivotType_start, TopPivot, VariablePivot, FixedPivot } PivotType;
typedef enum { ForwardSequence=SequenceDirection_start, BackwardSequence } SequenceDirection;
typedef enum { OneNorm=NormType_start, FrobeniusNorm, MaxAbsValue, InfinityNorm} NormType;
typedef enum { General=MatrixType_start, UpperTriangular, LowerTriangular, SymmetricUpper,
	       SymmetricLower, HermitianUpper, HermitianLower } MatrixType;

/* Following used in c02 chapter */
typedef struct
{
  Boolean ovflow, unflow;
} Ac02af;

typedef struct
{
  Integer expdep;
  double finity;
  Integer lrgexp;
  double sqrtfy, sqrtty, tiny, deps;
  Integer emin, emax;
} Bc02af;

typedef struct
{
  Boolean ovflow, unflow;
} Ac02ag;

typedef struct
{
  double finity, sqrtfy, sqrtty, tiny, deps;
  Integer eminm1, emaxm1, expdep, lrgexp;
} Bc02ag;

typedef struct
{
  double dpnewl, dpnewu, deps, temp, fact;
  Integer dbase, mnexp, mxexp, newl, newu;
} Cc02ag;


/* Following used in c06 chapter */
typedef enum {Nag_Convolution=Nag_VectorOp_start, Nag_Correlation} Nag_VectorOp;/*used in c06ekc*/
typedef enum {Nag_ForwardTransform=Nag_TransformDirection_start, Nag_BackwardTransform
} Nag_TransformDirection; 

/* Following used in d01 chapter */
typedef enum {
  Nag_Alg=Nag_QuadWeight_start, Nag_Alg_loga, Nag_Alg_logb, Nag_Alg_loga_logb
  } Nag_QuadWeight;

typedef enum
{
  Nag_UpperSemiInfinite=Nag_BoundInterval_start, Nag_LowerSemiInfinite, Nag_Infinite
  } Nag_BoundInterval;

typedef enum
{
  Nag_Cosine=Nag_TrigTransform_start, Nag_Sine
  } Nag_TrigTransform;

typedef enum
{
  Nag_Legendre=Nag_GaussFormulae_start, Nag_Rational, Nag_Laguerre, Nag_Hermite
  } Nag_GaussFormulae;

typedef enum
{
  Nag_OneIteration=Nag_MCMethod_start, Nag_ManyIterations
  } Nag_MCMethod;

/* Following used by d02 functions */
typedef enum {Nag_Relative=Nag_ErrorControl_start, Nag_Absolute, Nag_Mixed} Nag_ErrorControl;
/*typedef enum {Nag_InitMeshNotSet=-1, Nag_UserInitMesh,
		Nag_DefInitMesh} Nag_MeshSet;
*/
typedef enum {Nag_UserInitMesh=Nag_MeshSet_start,
		Nag_DefInitMesh} Nag_MeshSet;
typedef enum {Nag_RK_range=Nag_RK_task_start, Nag_RK_onestep} Nag_RK_task;
typedef enum {Nag_RK_2_3=Nag_RK_method_start, Nag_RK_4_5, Nag_RK_7_8} Nag_RK_method;
typedef enum {Nag_ErrorAssess_off=Nag_ErrorAssess_start, Nag_ErrorAssess_on} Nag_ErrorAssess;
/* typedef enum {Nag_InitSolDerivNotSet=-1, Nag_Sol, Nag_Der, Nag_SolDer}
Nag_SolDeriv;
*/
typedef enum {Nag_Sol=Nag_SolDeriv_start, Nag_Der, Nag_SolDer} Nag_SolDeriv;

/* Following used by e01 functions*/    
typedef enum {Nag_RC=Nag_2d_Scat_Method_start, Nag_Shep} Nag_2d_Scat_Method;


/* Following used by e02 functions*/    
typedef enum { Nag_LeftDerivs=Nag_DerivType_start, Nag_RightDerivs} Nag_DerivType;

/* Following used by d01 d02 e02 and e04 functions */
typedef enum { Nag_StartNotSet=Nag_Start_start  , Nag_Cold, Nag_Warm, Nag_Hot,
	       Nag_NewStart, Nag_ReStart, Nag_Continue } Nag_Start;


/* Following used by e04 and g13 functions (& some by h02) */
typedef enum {
Nag_PrintNotSet=Nag_PrintType_start  , Nag_NoPrint, Nag_Soln_Root, Nag_Soln, Nag_Iter, Nag_Iter_Long,
Nag_Soln_Root_Iter, Nag_Soln_Iter, Nag_Soln_Iter_Long, Nag_Soln_Iter_Const, 
Nag_Soln_Iter_Diag, Nag_Soln_Iter_Full
} Nag_PrintType;

typedef enum {
Nag_ChkNotSet=Nag_GradChk_start  , Nag_NoCheck, Nag_SimpleCheck, Nag_CheckObj, Nag_CheckCon,
Nag_CheckObjCon, Nag_XSimpleCheck, Nag_XCheckObj, Nag_XCheckCon,
Nag_XCheckObjCon
} Nag_GradChk;

typedef enum {
  Nag_D_NotSet=Nag_DPrintType_start  , Nag_D_NoPrint, Nag_D_Print, Nag_D_Sum, Nag_D_Full,
  Nag_D_Dbg
} Nag_DPrintType;

typedef enum {
Nag_ObjCheck=Nag_CheckType_start, Nag_ConCheck, Nag_DiffInt
} Nag_CheckType;

typedef enum {
Nag_DWant_NotSet=Nag_DWantType_start  , Nag_Grad_HessDiag, Nag_HessFull, Nag_Grad_HessFull
} Nag_DWantType;

typedef enum {
Nag_DInfo_NotSet=Nag_DerivInfo_start  , 
Nag_Deriv_OK, 
Nag_Fun_Constant, 
Nag_Fun_LinearOdd,
Nag_2ndDeriv_Large, 
Nag_1stDeriv_Small
} Nag_DerivInfo;

typedef enum {
Nag_DerivNotSet=Nag_DerivSet_start  , Nag_SomeG_SomeJ, Nag_AllG_SomeJ, Nag_SomeG_AllJ,
Nag_AllG_AllJ
} Nag_DerivSet;

typedef enum {
Nag_Deriv=Nag_FunType_start, Nag_NoDeriv
} Nag_FunType;

typedef enum {
Nag_LinFunNotSet=Nag_LinFun_start  , Nag_Lin_Deriv, Nag_Lin_NoDeriv
} Nag_LinFun;

typedef enum {
Nag_InitNotSet=Nag_InitType_start  , Nag_Init_None, Nag_Init_F_G_H, Nag_Init_All, Nag_Init_H_S
} Nag_InitType;

typedef enum {
Nag_BoundNotSet=Nag_BoundType_start  , Nag_Bounds, Nag_BoundsZero, Nag_BoundsEqual,
Nag_NoBounds, Nag_NoBounds_One_Call
} Nag_BoundType;

typedef enum {
Nag_ProbTypeNotSet=Nag_ProblemType_start  , Nag_FP, Nag_LP, Nag_QP1, Nag_QP2, Nag_QP3, Nag_QP4,
Nag_LS1, Nag_LS2, Nag_LS3, Nag_LS4, Nag_SparseFP, Nag_SparseLP, Nag_SparseQP,
Nag_SparseFPE, Nag_SparseFPL, Nag_SparseFPS, Nag_SparseQPP, Nag_SparseQPS	
} Nag_ProblemType;

typedef enum {
Nag_EndStateNotSet=Nag_EndState_start  , Nag_Feasible, Nag_Optimal, Nag_Deadpoint, Nag_Weakmin,
Nag_Unbounded, Nag_Infeasible, Nag_Too_Many_Iter, Nag_Hess_Too_Big, Nag_Cycling, 
Nag_Not_Converged, Nag_Not_Kuhn_Tucker, Nag_Deriv_Error, 
Nag_Hess_Indefinite, Nag_Basis_Ill_Cond, Nag_Basis_Singular, Nag_Out_Of_Workspace,

/* Successful outcomes */
Nag_MIP_Best_ISol, Nag_MIP_Stop_First_ISol, 

/* Succesfully executed but no solution */
Nag_MIP_No_ISol, 

/* Failure at root */
Nag_MIP_Root_Unbounded, Nag_MIP_Root_Infeasible, Nag_MIP_Root_Max_Itn, 
Nag_MIP_Root_Big_Hess, 

/* Integer solution found but could be suboptimal */
Nag_MIP_Max_Itn_ISol, Nag_MIP_Max_Nodes_ISol,
Nag_MIP_Big_Hess_ISol, Nag_MIP_Max_Depth_ISol,

/* No integer solution found but didn't complete succesfully */
Nag_MIP_Big_Hess_No_ISol, Nag_MIP_Max_Itn_No_ISol, Nag_MIP_Max_Nodes_No_ISol,
Nag_MIP_Max_Depth_No_ISol, 

Nag_MIP_User_Stop
} Nag_EndState;


typedef enum {
Nag_DegenInit=Nag_DegenJob_start, Nag_DegenOptimal, Nag_DegenEndCycle
} Nag_DegenJob;

typedef enum {
Nag_StatesInternal=Nag_SetStateMode_start, Nag_StatesExternal
} Nag_SetStateMode;

typedef enum { 
Nag_LU_TypeB=Nag_SparseQP_LU_Type_start, Nag_LU_TypeBS, Nag_LU_TypeBT 
} Nag_SparseQP_LU_Type;

typedef enum { 
Nag_LU_SolveL=Nag_SparseQP_LU_SolveType_start, Nag_LU_SolveB, Nag_LU_SolveBt
} Nag_SparseQP_LU_SolveType;

typedef enum { 
Nag_BasisFactB=Nag_SparseQP_BasisFactType_start, Nag_BasisFactBS 
} Nag_SparseQP_BasisFactType;

typedef enum { 
Nag_CrashNotSet=Nag_CrashType_start, Nag_NoCrash, Nag_CrashOnce, Nag_CrashTwice 
} Nag_CrashType;

typedef enum { 
Nag_ScaleNotSet=Nag_ScaleType_start, Nag_NoScale, Nag_RowColScale, Nag_ExtraScale 
} Nag_ScaleType;

/* Following used by f01qdc, f01qec */
typedef enum { Nag_ElementsIn=Nag_WhereElements_start, Nag_ElementsSeparate } Nag_WhereElements;

/* Following used by f02 functions */
typedef enum { Nag_Supplied=Nag_InitRotation_start, Nag_NotSupplied } Nag_InitRotation;
typedef enum {Nag_Select_Modulus=Nag_Select_Eigenvalues_start, Nag_Select_RealPart} Nag_Select_Eigenvalues;

/* Following used by f04mcc */
typedef enum { Nag_LDLTX=Nag_SolveSystem_start, Nag_LDX, Nag_DLTX, Nag_LLTX, Nag_LX, Nag_LTX
} Nag_SolveSystem;


typedef enum { Nag_SparseNsym_RGMRES=Nag_SparseNsym_Method_start, Nag_SparseNsym_CGS,
               Nag_SparseNsym_BiCGSTAB} Nag_SparseNsym_Method;

typedef enum { Nag_SparseNsym_ModFact=Nag_SparseNsym_Fact_start,
               Nag_SparseNsym_UnModFact } Nag_SparseNsym_Fact;

typedef enum { Nag_SparseNsym_NoPiv=Nag_SparseNsym_Piv_start, Nag_SparseNsym_UserPiv,
               Nag_SparseNsym_PartialPiv, Nag_SparseNsym_CompletePiv } Nag_SparseNsym_Piv;


typedef enum { Nag_SparseNsym_RemoveZeros=Nag_SparseNsym_Zeros_start, Nag_SparseNsym_KeepZeros, 
               Nag_SparseNsym_FailZeros } Nag_SparseNsym_Zeros;


typedef enum { Nag_SparseNsym_RemoveDups=Nag_SparseNsym_Dups_start, Nag_SparseNsym_SumDups, 
               Nag_SparseNsym_FailDups} Nag_SparseNsym_Dups;




typedef enum { Nag_SparseSym_CG=Nag_SparseSym_Method_start, Nag_SparseSym_Lanczos} Nag_SparseSym_Method;
typedef enum { Nag_SparseSym_ModFact=Nag_SparseSym_Fact_start, Nag_SparseSym_UnModFact } Nag_SparseSym_Fact;


typedef enum { Nag_SparseNsym_NoPrec=Nag_SparseNsym_PrecType_start, Nag_SparseNsym_SSORPrec,
               Nag_SparseNsym_JacPrec } Nag_SparseNsym_PrecType;



typedef enum { Nag_SparseSym_Weighted=Nag_SparseSym_Weight_start, Nag_SparseSym_UnWeighted} Nag_SparseSym_Weight;
typedef enum { Nag_SparseSym_Bisect=Nag_SparseSym_Bisection_start, Nag_SparseSym_NoBisect} Nag_SparseSym_Bisection;
typedef enum { Nag_SparseSym_OneNorm=Nag_SparseSym_Norm_start, Nag_SparseSym_TwoNorm,
               Nag_SparseSym_InfNorm} Nag_SparseSym_Norm;

typedef enum { Nag_SparseSym_TermLanczos=Nag_SparseSym_Term_start, Nag_SparseSym_TermLanczosCG } Nag_SparseSym_Term;
typedef enum { Nag_SparseSym_NoCheck=Nag_SparseSym_CheckData_start, Nag_SparseSym_Check } Nag_SparseSym_CheckData;


typedef enum { Nag_SparseSym_NoPiv=Nag_SparseSym_Piv_start, Nag_SparseSym_MarkPiv,
                                         Nag_SparseSym_UserPiv } Nag_SparseSym_Piv;


typedef enum { Nag_SparseSym_Prec=Nag_SparseSym_PrecType_start, Nag_SparseSym_NoPrec, Nag_SparseSym_SSORPrec,
               Nag_SparseSym_JacPrec } Nag_SparseSym_PrecType;

typedef enum { Nag_SparseSym_SumDups=Nag_SparseSym_Dups_start, 
               Nag_SparseSym_RemoveDups} Nag_SparseSym_Dups;

typedef enum { Nag_SparseSym_KeepZeros=Nag_SparseSym_Zeros_start, Nag_SparseSym_RemoveZeros} Nag_SparseSym_Zeros;


typedef struct Nag_Sparse_Comm_tag
{

  Integer idatax[20];
  double rdatax[20];
  Integer idata[20];
  double rdata[20];
  double *work;
  Integer *iwork;  
  Integer lwork;
  Integer liwork;
  
  Integer f11gac_lwreq;

  /* saved variables in unit f11bby */
  Integer f11bby_iter;
  Integer f11bby_jump;
  Integer f11bby_j;
  Integer f11bby_kase;  

  /* saved variables in unit f11gbz */

  Integer f11gbz_ib, f11gbz_id, f11gbz_ie2, f11gbz_infoch, f11gbz_ip, f11gbz_ir, f11gbz_iterm,
          f11gbz_itn, f11gbz_its, f11gbz_iw, f11gbz_iwgt, f11gbz_ix, f11gbz_maxitn, 
          f11gbz_maxits, f11gbz_monit, f11gbz_n, f11gbz_norm, f11gbz_sigcmp, f11gbz_sigcmx,
          f11gbz_stage;

  double f11gbz_anorm, f11gbz_bnorm, f11gbz_sigerr, f11gbz_sigerc, f11gbz_sigmax, f11gbz_sigmac,
         f11gbz_stplhs, f11gbz_stprhs, f11gbz_talpha, f11gbz_tbeta, f11gbz_tol, f11gbz_xnorm;

  Boolean f11gbz_floop, f11gbz_next;
  Nag_SparseSym_PrecType f11gbz_precon;


  /* saved variables in unit f11gby */

 double f11gby_alpha, f11gby_anorm, f11gby_beta1, f11gby_beta2, f11gby_betain, f11gby_bnorm, f11gby_bnorm2,
        f11gby_gammab, f11gby_pi, f11gby_rho1, f11gby_rho2,
        f11gby_rnorm2, f11gby_sigerr, f11gby_sigerc, f11gby_sigmax, f11gby_sigmac, f11gby_stplhs,
        f11gby_stprhs, f11gby_tol, f11gby_xnorm, f11gby_xnorm0, f11gby_zeta, f11gby_zetab;

 Integer f11gby_ib, f11gby_id, f11gby_ie2, f11gby_infoch, f11gby_iterm, f11gby_itn, f11gby_its,
         f11gby_iv1, f11gby_iv2, f11gby_iw1, f11gby_iw2, f11gby_iwgt, f11gby_ix,
         f11gby_maxitn, f11gby_maxits, f11gby_monit, f11gby_n, f11gby_norm, f11gby_resid,
         f11gby_sigcmp, f11gby_sigcmx, f11gby_stage;


 Boolean f11gby_floop, f11gby_next, f11gby_usexc;
 Nag_SparseSym_PrecType f11gby_precon;

  /* saved variable in unit f11gbx */
  
  Boolean f11gbx_fpass;

  /* saved variables in unit f11gbv */
  
  double f11gbv_rho, f11gbv_xi;

  /* saved variables in unit f11gbs */
  
  Boolean f11gbs_dov1;

  /* saved variables in unit f11gbr */
  
  double f11gbr_deltab, f11gbr_epslon;


  /* saved variables in unit f11baz */
  
  Integer f11baz_idatax[20];
  double f11baz_rdatax[20];

  /* saved variables in unit f11gbw */
  
  double f11gbw_eps, f11gbw_epslon, f11gbw_mu1, f11gbw_mu1l, f11gbw_mun, f11gbw_munh, 
         f11gbw_sigmx0, f11gbw_sigmx1, f11gbw_sigtol, f11gbw_ta2, f11gbw_tb2, f11gbw_tnorm;

  Integer f11gbw_maxits;

  Boolean f11gbw_sigcmp;

  /* saved variables in unit f11gbf */

   Integer f11gbc_ifaill, f11gbc_ifailm, f11gbc_info, f11gbc_irevcx, f11gbc_kill,
           f11gbc_idata[20], f11gbc_lworkl;

   Boolean f11gbc_done, f11gbc_first;

   double f11gbc_rdata[20];
 
  /* saved variables in unit f11gbq */

   double f11gbq_factol, f11gbq_stprhx, f11gbq_tolc, f11gbq_toleps, f11gbq_tolf, f11gbq_toll,
          f11gbq_xcbnrm, f11gbq_xlbnrm, f11gbq_xnorm;
   Integer f11gbq_loop;
   Boolean f11gbq_fterm;
   
  /* saved variables in unit f11gbu */
   
   double f11gbu_epslon, f11gbu_stprhx, f11gbu_toleps, f11gbu_tolf, f11gbu_xnorm0;
   Integer f11gbu_loop, f11gbu_nrestl, f11gbu_nrestc;
   Boolean f11gbu_fterm;



  /* saved variables in unit f11bbz */

  Integer f11bbz_iwgt, f11bbz_norm, f11bbz_iterm, f11bbz_monit,  f11bbz_m, f11bbz_n,
    f11bbz_stage, f11bbz_ib, f11bbz_ic,  f11bbz_ih,
   f11bbz_iq, f11bbz_ir, f11bbz_is, f11bbz_it, f11bbz_ix, f11bbz_infoch,  f11bbz_imonit, 
   f11bbz_maxitn,  f11bbz_itn;

  double f11bbz_anorm, f11bbz_bnorm, f11bbz_rnrm20, f11bbz_xnorm, f11bbz_rnorm2, f11bbz_xnorm2,
  f11bbz_sigmax, f11bbz_stplhs, f11bbz_stprhs, f11bbz_tol;

  Boolean f11bbz_anrcmp, f11bbz_precon, f11bbz_next, f11bbz_fcall;

 /* saved variables in unit f11bbu */

  Boolean f11bbu_first, f11bbu_usedot;

 /* saved variable in unit f11bbv */

   double f11bbv_tolf;

/* saved variables in unit f11bbu */

  double f11bbr_alpha, f11bbr_rho;
  Boolean f11bbr_do1;

/* saved variables in unit f11bbw */


  double  f11bbw_stpr;
  Integer f11bbw_m1;
  double  f11bbw_tnorm;
  Integer f11bbw_im;
  Boolean f11bbw_sigcmp;
  Integer f11bbw_im1;
  double  f11bbw_tinorm;
  Boolean f11bbw_restrt;
  double  f11bbw_tol1;

/* saved variables in unit f11bbt */

  Integer f11bbt_iwgt, f11bbt_norm;
  Boolean f11bbt_next;
  Integer f11bbt_n;
  Boolean f11bbt_fcall;
  Integer f11bbt_irbar, f11bbt_stage;
  double  f11bbt_anorm, f11bbt_bnorm;
  Integer f11bbt_iterm, f11bbt_monit;
  double  f11bbt_xnorm;
  Integer f11bbt_ib, f11bbt_ip, f11bbt_ir, f11bbt_it, f11bbt_ix, f11bbt_iy, f11bbt_infoch;
  Boolean f11bbt_anrcmp, f11bbt_precon;
  double  f11bbt_sigmax;
  Integer f11bbt_maxitn;
  double  f11bbt_stplhs, f11bbt_stprhs;
  Boolean f11bbt_restrt;
  Integer f11bbt_itn;
  double  f11bbt_tol;

/* saved variables in unit f11bbq */

  double  f11bbq_tolf;
  Integer f11bbq_loop;
  Boolean f11bbq_pass1;
  Boolean f11bbq_frest;
  double  f11bbq_xnmax, f11bbq_toleps, f11bbq_stplhx, f11bbq_stprhx, f11bbq_zz1, f11bbq_zz2;

/* saved variables in unit f11bbp */

  Boolean f11bbp_doxn;
  Boolean f11bbp_next;
  Boolean f11bbp_fcall;
  Boolean f11bbp_anrcmp, f11bbp_finish, f11bbp_precon;
  Boolean f11bbp_restrt;

  Integer f11bbp_imonit, f11bbp_maxitn;
  Integer f11bbp_igam, f11bbp_itau;
  Integer f11bbp_iwgt, f11bbp_norm;
  Integer f11bbp_m, f11bbp_n;
  Integer f11bbp_igamp, f11bbp_irbar, f11bbp_stage;
  double  f11bbp_xnorm;
  double  f11bbp_anorm, f11bbp_bnorm;

  Integer f11bbp_iterm, f11bbp_monit;
  Integer f11bbp_ib, f11bbp_iq, f11bbp_ir, f11bbp_ix, f11bbp_infoch;
  Integer f11bbp_idx, f11bbp_itn;
  Integer f11bbp_ipx;

  double f11bbp_sigmax;
  double f11bbp_stplhs, f11bbp_stprhs;
  double f11bbp_tol;


/* saved variables in unit f11bbn */

  double  f11bbn_rhox, f11bbn_alpha, f11bbn_tol1, f11bbn_tol2, f11bbn_rn1;
  Integer f11bbn_i1, f11bbn_j1, f11bbn_n1, f11bbn_im, f11bbn_jn1, f11bbn_imn;
  Boolean f11bbn_doq;

/* saved variables in unit f11bbm */

  double  f11bbm_tolf;
  Integer f11bbm_jitn, f11bbm_kitn, f11bbm_loop;
  Boolean f11bbm_pass1;
  Boolean f11bbm_frest;
  Integer f11bbm_irest;
  double  f11bbm_xnmax, f11bbm_toleps, f11bbm_stplhx, f11bbm_stprhx, f11bbm_zz1, f11bbm_zz2;


/* saved variables in unit f11bbf */

  Integer f11bbf_irevcx;
  Boolean f11bbf_done;
  Boolean f11bbf_first;
  Integer f11bbf_info, f11bbf_kill;
  Integer f11bbf_idata[20];
  double  f11bbf_rdata[20];
  Integer f11bbf_ifaill, f11bbf_ifailm, f11bbf_lworkl;


} Nag_Sparse_Comm;


/* Following used by g01 chapter */
typedef enum {Nag_LowerTail=Nag_TailProbability_start, Nag_UpperTail, Nag_TwoTailSignif,
		Nag_TwoTailConfid, Nag_TwoTail}
Nag_TailProbability;

typedef enum {Nag_RankScores=Nag_Scores_start, Nag_NormalScores, Nag_BlomScores,
		Nag_TukeyScores, Nag_WaerdenScores, Nag_SavageScores}
Nag_Scores;

typedef enum { Nag_AverageTies=Nag_Ties_start, Nag_LowestTies, Nag_HighestTies,
Nag_RandomTies, Nag_IgnoreTies} Nag_Ties;

/* Following used by g02 chapter */
typedef enum {Nag_WeightedEstimate=Nag_IncludeWeight_start, Nag_UnweightedEstimate} Nag_IncludeWeight;
typedef enum {Nag_MeanInclude=Nag_IncludeMean_start, Nag_MeanZero} Nag_IncludeMean;
typedef enum {Nag_ObservAdd=Nag_UpdateObserv_start, Nag_ObservDel}   Nag_UpdateObserv;
typedef enum {Nag_AboutMean=Nag_SumSquare_start, Nag_AboutZero} Nag_SumSquare;
typedef enum {Nag_FirstCall=Nag_Initialize_start, Nag_Update} Nag_Initialize;
typedef enum {Nag_Expo=Nag_Link_start, Nag_Iden, Nag_Log, Nag_Sqrt, Nag_Reci, 
	      Nag_Logistic, Nag_Probit, Nag_Compl } Nag_Link;
typedef enum {Nag_RegNotSet=Nag_RegType_start  , Nag_HuberReg, Nag_MallowsReg, Nag_SchweppeReg}
              Nag_RegType;
typedef enum {Nag_PsiNotSet=Nag_PsiFun_start  , Nag_Lsq, Nag_HuberFun, Nag_HampelFun, Nag_AndrewFun, 
              Nag_TukeyFun} Nag_PsiFun;
typedef enum {Nag_SigmaNotSet=Nag_SigmaEst_start  , Nag_SigmaRes, Nag_SigmaConst, Nag_SigmaChi}
              Nag_SigmaEst;
typedef enum {Nag_CovNotSet=Nag_CovMatrixEst_start  , Nag_CovMatAve, Nag_CovMatObs} Nag_CovMatrixEst;
typedef enum {Nag_SigmaSimul=Nag_SigmaSimulEst_start, Nag_SigmaBypas} Nag_SigmaSimulEst;

/* Following used by g03 chapter */
typedef enum {Nag_MatCorrelation=Nag_PrinCompMat_start, Nag_MatStandardised, Nag_MatSumSq,
		Nag_MatVarCovar} Nag_PrinCompMat;
typedef enum {Nag_ScoresStand=Nag_PrinCompScores_start, Nag_ScoresNotStand, Nag_ScoresUnitVar,
		Nag_ScoresEigenval} Nag_PrinCompScores;
typedef enum {Nag_NoWeights=Nag_Weightstype_start, Nag_Weightsfreq, Nag_Weightsvar} Nag_Weightstype;
typedef enum {Nag_RoLoadStand=Nag_RotationLoading_start, Nag_RoLoadNotStand} Nag_RotationLoading;
typedef enum {Nag_NoTransNorm=Nag_TransNorm_start, Nag_Orig, Nag_OrigCentroid, 
		Nag_Norm, Nag_OrigNorm, Nag_OrigNormCentroid}
Nag_TransNorm;
typedef enum {Nag_LsqScale=Nag_RotationScale_start, Nag_NotLsqScale} Nag_RotationScale;
typedef enum {Nag_DataCorr=Nag_FacMat_start, Nag_DataCovar, Nag_MatCorr_Covar} Nag_FacMat;
typedef enum {Nag_FacScoreRegsn=Nag_FacScoreMethod_start, Nag_FacScoreBart} Nag_FacScoreMethod;
typedef enum {Nag_FacRotate=Nag_FacRotation_start, Nag_FacNoRotate} Nag_FacRotation;
typedef enum {Nag_EqualCovar=Nag_GroupCovars_start, Nag_NotEqualCovar} Nag_GroupCovars;
typedef enum {Nag_SamplePoints=Nag_MahalDist_start, Nag_GroupMeans} Nag_MahalDist;
typedef enum {Nag_DiscrimEstimate=Nag_DiscrimMethod_start, Nag_DiscrimPredict} Nag_DiscrimMethod;
typedef enum {Nag_EqualPrior=Nag_PriorProbability_start, Nag_GroupSizePrior, Nag_UserPrior} Nag_PriorProbability;
typedef enum {Nag_MatUp=Nag_MatUpdate_start, Nag_NoMatUp} Nag_MatUpdate;
typedef enum {Nag_DistAbs=Nag_DistanceType_start, Nag_DistEuclid, Nag_DistSquared} Nag_DistanceType;
typedef enum {Nag_VarScaleStd=Nag_VarScaleType_start, Nag_VarScaleRange, Nag_VarScaleUser, Nag_NoVarScale} Nag_VarScaleType;


typedef enum {Nag_SingleLink=Nag_ClusterMethod_start, Nag_CompleteLink, Nag_GroupAverage,
	              Nag_Centroid, Nag_Median, Nag_MinVariance} Nag_ClusterMethod;

typedef enum {Nag_DendNorth=Nag_DendOrient_start, Nag_DendSouth,
                    Nag_DendEast, Nag_DendWest} Nag_DendOrient;

typedef enum {Nag_AllEigVals=Nag_Eigenvalues_start, Nag_LargeEigVals} Nag_Eigenvalues;

typedef enum {Nag_Stress=Nag_ScaleCriterion_start, Nag_SStress} Nag_ScaleCriterion;

/* Following used by g04 chapter */
typedef enum {Nag_NoBlocks=Nag_Blocks_start, Nag_SerialBlocks, Nag_ParallelBlocks} Nag_Blocks;









/* Following used by g05 chapter */
typedef enum {Nag_PDF=Nag_DiscreteDistrib_start, Nag_CDF} Nag_DiscreteDistrib;

/* Following used by g07 chapter */
typedef enum {Nag_PopVarEqual=Nag_PopVar_start, Nag_PopVarNotEqual} Nag_PopVar;

/* Following used by g10 chapter */
typedef enum {Nag_4253H=Nag_Smooth_Type_start, Nag_3RSSH} Nag_Smooth_Type;

/* Following used by g12 chapter */
typedef enum {Nag_Freq=Nag_FreqTime_start, Nag_NoFreq} Nag_FreqTime;

/* Following used by g13 chapter */
typedef enum {Nag_CriteriaNotSet=Nag_Likelihood_start  , Nag_LeastSquares,
Nag_Exact, Nag_Marginal
} Nag_Likelihood;

typedef enum {Nag_Rectangular=Nag_LagWindow_start, Nag_Bartlett,
Nag_Tukey, Nag_Parzen
} Nag_LagWindow;

typedef enum {Nag_NoCorrection=NagMeanOrTrend_start, Nag_Mean, Nag_Trend}
NagMeanOrTrend;

typedef enum {Nag_Unlogged=Nag_LoggedSpectra_start, Nag_Logged} Nag_LoggedSpectra;

/* Following used by g13 Kalman filter routines */
typedef enum {Nag_next_state=Nag_state_start, Nag_curr_state} Nag_state;
typedef enum {Nag_ab_prod=Nag_ab_input_start, Nag_ab_sep} Nag_ab_input;
typedef enum {Nag_UH_Observer=Nag_ObserverForm_start, Nag_LH_Observer} Nag_ObserverForm;
typedef enum {Nag_UH_Controller=Nag_ControllerForm_start, Nag_LH_Controller} Nag_ControllerForm;

/* h02 chapter */
typedef enum 
{
  Nag_OutputNotSet=Nag_OutputType_start  , Nag_NoOutput, Nag_MPS_Summary, Nag_MPS_List
} Nag_OutputType;


typedef enum
{
  Nag_BB_OK=Nag_BB_Fail_start, Nag_BB_Alloc_Fail, Nag_BB_Internal_Error
} Nag_BB_Fail;

typedef enum
{
  Nag_MIP_TypeNotSet=Nag_MIP_ProbType_start  , Nag_MILP, Nag_MIQP1,
  Nag_MIQP2, Nag_MIQP3, Nag_MIQP4  /* & sparse cases? */
} Nag_MIP_ProbType;

typedef enum
{
  Nag_NodSel_NotSet=Nag_Node_Selection_start, Nag_Deep_Search, Nag_Broad_Search,
  Nag_MinObj_Search, Nag_DeepBroad_Search, Nag_DeepMinObj_Search
} Nag_Node_Selection;

typedef enum
{
  Nag_VarSel_NotSet=Nag_Var_Selection_start, Nag_First_Int, Nag_Nearest_Half,
  Nag_Use_Priority
} Nag_Var_Selection;

typedef enum
{
  Nag_BrDir_NotSet=Nag_Branch_Direction_start, Nag_Branch_Left, Nag_Branch_Right,
  Nag_Branch_InitX
} Nag_Branch_Direction;

typedef enum
{
  Nag_NS_NotSet=Nag_NodeStatus_start  , Nag_NS_NotSolved, Nag_NS_NotBranched, 
  Nag_NS_Integer, Nag_NS_Bounded, Nag_NS_Infeasible, Nag_NS_Terminated  
} Nag_NodeStatus;


/* m01 chapter */
typedef enum {Nag_Ascending=Nag_SortOrder_start, Nag_Descending} Nag_SortOrder;
typedef enum {Nag_First=Nag_SearchMatch_start, Nag_Last} Nag_SearchMatch;

/* x04 chapter */
typedef enum { Nag_HashOK=Nag_HashError_start, Nag_HashTableTooBig} Nag_HashError;
  

/* Nag Error Structure */
#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_ERRHAN)(char NAG_HUGE *, int, char NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_ERRHAN)();
#endif

typedef struct {
  int code;                     /* Out: Error Code */
  Boolean print;                        /* In: print? yes/no */
  char message[NAG_ERROR_BUF_LEN];      /* InOut: Error message */
  NAG_ERRHAN handler;           /* In: Error handling function */
  Integer errnum;               /* May hold useful value for some errors */
} NagError;

/* Structure containing void pointer for user to communicate
 * to a user defined function, used by d02's.
 */
typedef struct
{
    Pointer p;
} Nag_User;

typedef struct comm_struct {
Integer flag;   /* user error flag set to 0 before each call to lsqfun, test
		   for != 0 on exit */
Boolean first;  /* TRUE on first call to objective function, FALSE after */
Integer nf;     /* output: number of function evaluations or number of times
		 residuals have been calculated. Calls to objective function */
Boolean it_prt; /* Nag_Iteration print */
Boolean it_maj_prt; /* Nag major iteration print (e04ucc) */
Boolean sol_sqp_prt;/* Solution print (e04ucc)*/
Boolean sol_prt;/* Solution print or, for e04ucc, major itn. sol. print */
Boolean rootnode_sol_prt; /* Solution at rootnode of Branch & Bound */
Boolean node_prt; /* Summary of node solution in B&B */
Boolean rootnode_prt;  /* Summary of root node solution in B&B */
Boolean g_prt;  /* Derivative check info. from e04ucy available */
Boolean new_lm; /* New Lagrange multipliers calculated */


/* User workspace pointers */
Pointer p;      /* generic pointer to user workspace */
Integer NAG_HUGE *iuser; /* user workspace */
double NAG_HUGE *user;   /* user workspace */ 

/* -- Members used internally -- */

/* Nag workspace pointers */
double NAG_HUGE *nag_w;  /* Nag workspace */
Integer NAG_HUGE *nag_iw;/* Nag integer workspace */
Pointer nag_p;  /* Nag generic workspace */

   /* Following are specifically intended to provide workspace
      for a print function (e.g. to pass array of values, constructed 
      within the function, to user-defined print function).  The ptrs 
      should be used *only* to provide wkspace local to a print fun 
      & not to provide any kind of persistent data storage (nag_w, 
      nag_iw do not always follow this convention; hence the need 
      for the new members.) */

double *nag_print_w;
Integer *nag_print_iw;

Integer nrealloc;  /* Counter for reallocations in sparse optimizers */
} Nag_Comm;

/* structure for specifying output stream to
 * message printing function x04baz,
 * not visible to user */
typedef struct
{
    Boolean set;    /* Has stream been set before? */
    Boolean st_out; /* TRUE if stdout stream being used */
    Boolean st_err; /* TRUE if stderr stream being used */
    Boolean st_in;  /* TRUE if stdin stream being used */
    char name[NAG_FILE_LEN]; /* filename to which message to be output */
    char NAG_HUGE *chapter;  /* points to string identifying chapter in use */
    int error;     /* non-zero denotes an error has occured
		    * during use */
} Nag_FileSt;

/* User visible structure for carrying
 * message printing information
 */
typedef struct
{
    Integer init_mesg1;
    Integer init_mesg2;
    char NAG_HUGE *res_file; /* filename to which results to be output */
    Boolean print;  /* default printing if TRUE */
    int code;       /* message code */
    char NAG_HUGE *name;     /* user specified name of the nag routine
		     * which the user is calling */
    char NAG_HUGE *opt_name; /* user specified name of option set being read */
#ifdef NAG_PROTO
   /* In: User defined message printing function */
    void (*print_mesg)(char NAG_HUGE *str, int code, char NAG_HUGE *name);
#else
    void (*print_mesg)();
#endif
} Nag_Mesg;
	    
/* d01 structure types */
typedef struct
{
  Integer num_subint;
  Integer fun_count;
  double NAG_HUGE *sub_int_beg_pts;
  double NAG_HUGE *sub_int_end_pts;
  double NAG_HUGE *sub_int_result;
  double NAG_HUGE *sub_int_error;
} Nag_QuadProgress;

typedef struct
{
  Integer intervals;
  Integer fun_count;
  Integer subints_per_interval;
  double NAG_HUGE *interval_error;
  double NAG_HUGE *interval_result;
  Integer NAG_HUGE *interval_flag; 
} Nag_QuadSubProgress;


/* Structure types for d02 functions */

typedef struct
{
    double alpha[12], beta[12], psi[12], v[12], w[12], sig[13], g[13], gi[11];
    double xold, hold, told, xsave, tstar, twou;
    Integer init, ibegin, itol, iinteg, itstop, inflop, iquit, iv[10], ns,
	   kord, kold, ksteps, kle4, kprev, ivc, kgi;
    Boolean start, phase1, nornd, stiff, intout;
} Nag_ad02;



typedef struct
{
    double zero, u, fouru, sru, u34, u78, srbig, delsgn, troots, tleft, svtnew;
    Integer krootp, inrotp;
    Boolean gstop, pgstop, root, roots, needg, discop, newgeq, search, pserch;
} Nag_bd02;


typedef struct
{
  Integer nfun, nimp;
  Boolean repeat;
  Integer l1, l2, ng, n1, nc;
} Nag_ad02ra; /* Introduced at Mark 5 to make d02rac, d02gac, d02gbc 
	       * thread-safe.
	       */

typedef struct
{
    double ts, gts, fd, ftd, fm, fmx, prev, min, max, avg, amb, ambs, 
	    sbma, sf, dezu34;
    Integer ic, kmc, ising, kount, ktest;
    Boolean locext;
} Nag_dd02;

typedef struct {
    double aez, agtrl, at, atm, atmtr, brd, brds, ddt, deltg, deltgm, 
	    deltr, dg, dtg, dtgk, dtgp, dtqr, fac, g21, g32, gat, gbr, gbrc, 
	    gbrmin, gbrmns, gct, glocmx, gmsign, goldsv, got, gres, gright, 
	    gt21, gt32, gtint, gtqr, gtrl, gtt, gttmr, gvalue, gvdif, gvmax, 
	    gvmin, omu78, opu78, pu, qc0, qc1, qc2, qd, qdq, qr1, qr1n, qr1p, 
	    quadmn, rd, rdg, rem, remd, rez, rootno, sdg, signg, slope1, 
	    slope2, slopep, sotgv, ss, sstu78, t21, t31, t32, taddon, tcb, 
	    tcf, tclose, tdif, tdmn, tdmx, tfar, tgbrmn, tgvdif, tgvmt, tint, 
	    tintp, tk, tkold, tlbk, tnew, tp, tpoint, tprev, tqr, tqr1, tqr2, 
	    tqrf, tqrmp, tqrp, tqs, tqt, trd, trj, trl, trn, trnow, troot, 
	    tsave, tslop, tsru, tstep, tt, tu78, tut, tz, tzero, tzerok, tzj, 
	    tzk;
    double grno, gt, sttu78, tadd[2], tr[3];
} Nag_ed02;

typedef struct {
    Boolean backr, mroot, multr, peak, qrreal, skpmmt;
} Nag_fd02;

typedef struct {
    Integer i, i1, i2, i21, i3, iadd, ic1, ic2, icase, icbr, icbrk, idtg, iga,
	     ijk, ik, ikmr, iktqr, indx, isi, isin, isr, itgvt, itry, j, j1, 
	    j2, j23, j3, jp1, jr, jsr, k, kge, kp1, krootn, krooto, kroots, 
	    kskpmt, ktry, l, ltqbig, ltr, mmr, mmrin, mmrk, nsr, ntin, ntinm;
    int ipath;
} Nag_gd02;

typedef struct
{
    double a, absdx, da, delx, dfdub, dfdxb, dx, fbnd, power, relper, ynorm,
	   ypnorm;
    Integer icase, k, lk;
} Nag_jd02;

typedef struct
{
    double two[13], gstr[13], absh, big, erk, erkm1, p5eps, round, u;
    Integer ifail, km1, km2, knew, kp1, kp2;
} Nag_kd02;

typedef struct Nag_Adams
{
  Integer init1;  /* initialisation checks */
  Integer init2;
  Integer neqf;  /* number of o.d.e's */
  Integer neqg;  /* number of event functions defined for root finding */
  Nag_Start state; /* new system, continuation current system or restart */
  Boolean alter_g; /* TRUE if event functions redefined */
  Boolean sophist; /* type of root search technique */
  Boolean vectol;  /* vector or scalar error control */
  double NAG_HUGE *a_tol;    /* absolute tolerance(s) */
  double NAG_HUGE *rtol;    /* relative tolerance(s) */
  Boolean crit;    /* if integration is not to go beyond tcrit */
  double tcrit;    /* value beyond which integration must not go */
  Boolean one_step; /* return after one step or after interval */
  double hmax;      /* if not 0.0 hmax is bound on absolute step size */
  Integer max_step; /* maximum number of attempted steps */
  Boolean root;    /* TRUE if root found */
  Integer badcmp; /* also in fail.errnum, index of zero error weight */
  double hnext;   /* next step */
  double hlast;   /* last step */
  double tolfac;  /* tolerance scale factor */
  double tcurr;   /* current value of t */
  Integer nsuccess;  /* number of successful steps since start */
  Integer nfail;    /* number of failed steps since start */
  Integer ord_last;  /* order of method used last successfully used */
  Integer ord_next;  /* order of method to be attempted next */
  Integer NAG_HUGE *events;  /* array pointer, info. about event functions */
  double NAG_HUGE *resids;  /* array pointer, value of kth event function at root */
  double NAG_HUGE *yp;  /* approx. derivative of y[] from integration */
  double NAG_HUGE *yy;
  double NAG_HUGE *p;
  double NAG_HUGE *phi;
  
  Integer qwc_state; /* whether error occured in call to d02qwc */
  Integer qfc_state; /* whether integrator called before d02qzc */
  Integer NAG_HUGE *giwork, *siwork;     /* Workspace pointers  */
  double NAG_HUGE *w, *gwork, *swork;
#ifdef NAG_PROTO
  void (NAG_CALL *free)(struct Nag_Adams NAG_HUGE *opt);
#else
  void (*free)();
#endif
  Integer index;  /* index of event equation for which root detected */
  Integer type;   /* type of root */
  Nag_ad02 ad02qf;
  Nag_bd02 bd02qf;
  Boolean crash;
  Integer izflag;
  Nag_dd02 dd02qf;
  Nag_ed02 ed02qf;
  Nag_fd02 fd02qf;
  Nag_gd02 gd02qf;
  Integer kroo;
  Nag_jd02 jd02qf;
  Nag_kd02 kd02qf;
} Nag_ODE_Adams;



typedef struct
{
    double tstrt, tnd, dir, hstrt, tolr;
    Integer neqn;
} Nag_ad02pd;

typedef struct
{
    double t, h, told, hold;
    Integer nfcn, svnfcn, okstp, flstp;
    Boolean first, last;
} Nag_bd02pd;

typedef struct
{
    Integer prthrs, prerst, prwt, pryold, prscr, pry, pryp, 
    prstgs, printp, lnintp;
} Nag_cd02pd;

typedef struct
{
    double a[169], b[13], c[13], bhat[13], r[66], e[7];
    Integer ptr[13], nstage, methd, mintp;
    Boolean intp;
} Nag_dd02pd;

typedef struct
{
    double toosml, cost, safety, expon, stbrad, tanang, rs, rs1, rs2, rs3, rs4;
    Integer order, lststg, maxtry, nsec;
    Boolean fsal;
} Nag_ed02pd;

typedef struct
{
    double maxerr, locmax;
    Integer gnfcn, przstg, przy, przyp, przers, przerr, przynu;
    Boolean erason, erasfl;
} Nag_fd02pd;

typedef struct
{
    double mcheps, dwarf, rndoff, sqrrmc, cubrmc, tiny;
    Integer outch;
} Nag_gd02pd;

typedef struct
{
    Boolean utask;
} Nag_hd02pd;

typedef struct
{
    char rec[800];
} Nag_jd02pd;

typedef struct Nag_RK
{
#ifdef NAG_PROTO
  void (*free) (struct Nag_RK NAG_HUGE *);
#else
  void (*free) ();
#endif
  double NAG_HUGE *work;
  Nag_ad02pd ad02pd;
  Nag_bd02pd bd02pd;
  Nag_cd02pd cd02pd;
  Nag_dd02pd dd02pd;
  Nag_ed02pd ed02pd;
  Nag_fd02pd fd02pd;
  Nag_gd02pd gd02pd;
  Nag_hd02pd hd02pd;
  Nag_jd02pd jd02pd;

Integer init1;          /* flags indicating structure has been initialised, */
Integer init2;          /* set to IDUMMY and INIT2DUMMY = -23456 on initialisation */

/* d02pcc saved variables */
   double utend;
   double tlast;

/* d02pdc saved variables */
  Boolean chkeff;
  Integer ntend;
  double havg;
  double errold;
  Boolean phase2;
  Integer ynew;
  Integer ypold;
  Integer jflstp;

/* d02pxc saved variables */
  Integer nwntsv;
  Integer startp;
  Boolean inintp;

/* d02pdm saved variables */
  Integer svsta[8];

/* d02pyc arguments */
  Integer totfcn;
  Integer  stpcst;
  double waste;
  Integer stpsok;
  double hnext;
/* End d02pyc arguments */

/* d02pxc work array pointer and saved size of array */
  double NAG_HUGE *wrkint;
  Integer save_nwant;

} Nag_ODE_RK;

/* Structures used by d02ejc */

typedef struct {
    Integer itrace, idev;
} Nag_ad02nm;

typedef struct {
    double hold;
    Integer nq, nqu, nst, nre, nje, niter, ninter, kcur;
} Nag_bd02nm;

typedef struct {
    double hmxstt;
} Nag_cd02nm;

typedef struct {
    double h, el0, tn, hu, hmin, hmxi, dsqu, els, wb, elc, ts;
    Integer init, lewt, lacor, lsavr, ldae, mxstep, mxhnil, nhnil, nslast;
} Nag_dd02nm;

typedef struct {
    Integer ifunc, istep, ier, iewset, iretur, idacnt, idaold, initl, iode, 
	    icrash, ifn, idachk;
} Nag_ed02nm;

typedef struct {
    double dunflo, uround;
    Integer iovflo;
} Nag_fd02nm;

typedef struct {
    double damp, rjnorm, crate;
    Integer maxit;
} Nag_gd02nm;

typedef struct {
    Integer inorm;
} Nag_hd02nm;

typedef struct {
    double big, qtymin, rootn;
} Nag_jd02nm;

typedef struct {
    Boolean ihit;
} Nag_kd02nm;

typedef struct {
    double a_toli, big, ewti, h0, rh, rtoli, size, tcrit, tem, hnext, 
	    tolsf, tp, hmax;
    Integer i, ifj, iflag, ichmax, imxer, j, kgo, lenrw, mxhnl0, mxstp0, n, 
	    isave, nold;
    Boolean crit1, crit2, nt1stp;
} Nag_ld02nm;

typedef struct {
    char singlr[7];
} Nag_md02nm;

typedef struct {
    Integer irow, icol;
} Nag_nd02nm;

typedef struct {
    double con, di, fac, hl0, r, r0, srur, yi, yj, yjj;
    Integer i, ier, ii, j, j1, jj, n;
} Nag_rd02nm;

typedef struct {
    Integer iwp;
} Nag_sd02nm;

typedef struct {
    double rdm10[7], bdel, rdm11, bel1h, rdm12, bdcon, rdm13[34];
    Integer idm11[17];
} Nag_yd02nm;

typedef struct {
    char odcode[7];
} Nag_zd02nm;

typedef struct {
    double const1, const2, const3, const4, const5, const6;
    Boolean start, speclh;
} Nag_nd02nn;

typedef struct {
    Integer nzeros;
    Boolean nonzer, somzer, allzer;
} Nag_pd02nn;

typedef struct {
    double hprop;
    Boolean change;
} Nag_rd02nn;

typedef struct {
    double hold;
} Nag_sd02nn;

typedef struct {
    double errcwe;
} Nag_td02nn;

typedef struct {
    double ddn, dsm, dup, eljh, terk, terkm1, terkp1, exdn, exsm, exup, r,
	     rh, rhdn, rhsm, rhup, told, terkm2, ddn1;
    Integer i, jc, iredo, iret, j, jb, kgo, ncf, iz, jstart, kflag;
    Boolean raise;
} Nag_ud02nn;

typedef struct {
    double dcon, del, delp, d1, d2, el1h;
    Integer maxcor, isave, m, iz, i;
} Nag_yd02nn;

typedef struct {
    double rownd, conit, el[13], elco[156]      /* was [13][12] */, rmax, 
	    tesco[36]   /* was [3][12] */, ccmax, rc;
    Integer ialth, ipup, lmax, meo, nslp, icf, jcur, l, meth, miter, maxord, 
	    msbp, mxncf, n;
} Nag_ad02xk;

typedef struct {
    double hhused;
    Boolean petzld, fniter;
} Nag_bd02xk;

typedef struct {
    char jceval[7], outopt[7], gopt[7];
} Nag_cd02ej;

typedef struct Nag_BDF
{
  Nag_ad02nm ad02nm;
  Nag_bd02nm bd02nm;
  Nag_cd02nm cd02nm;
  Nag_dd02nm dd02nm;
  Nag_ed02nm ed02nm;
  Nag_fd02nm fd02nm;
  Nag_gd02nm gd02nm;
  Nag_hd02nm hd02nm;
  Nag_jd02nm jd02nm;
  Nag_kd02nm kd02nm;
  Nag_ld02nm ld02nm;
  Nag_md02nm md02nm;
  Nag_nd02nm nd02nm;
  Nag_rd02nm rd02nm;
  Nag_sd02nm sd02nm;
  Nag_yd02nm yd02nm;
  Nag_zd02nm zd02nm;

  Nag_nd02nn nd02nn;
  Nag_pd02nn pd02nn;
  Nag_rd02nn rd02nn;
  Nag_sd02nn sd02nn;
  Nag_td02nn td02nn;
  Nag_ud02nn ud02nn;
  Nag_yd02nn yd02nn;

  Nag_ad02xk ad02xk;
  Nag_bd02xk bd02xk;

  Nag_cd02ej cd02ej;
  Integer icall; /* Thread safety, icall was declared static in d02nmu */
} Nag_ODE_BDF;

/* End of structures used by d02ejc */


/* Structures used by e01 functions */
typedef struct{
  Integer init1;  /* initialisation checks */
  Integer init2;
  Integer exit_status; /* set to magic number if non-fatal error during
			  set up */
  Nag_2d_Scat_Method method;
  Integer m;
  double  rnw;
  Boolean mem_set;
  double  NAG_HUGE *x;
  double  NAG_HUGE *y;
  double  NAG_HUGE *f;
  double  NAG_HUGE *fnodes;
  Integer NAG_HUGE *triang;
  double  NAG_HUGE *grads;
} Nag_Scat_Struct;

typedef struct{
  double  rnw;
  double  rnq;
  Integer nw;
  Integer nq;
  Integer minnq;
} Nag_E01_Opt;


/* Structures used by e02 functions */
typedef struct{
  Integer n;
  double NAG_HUGE *lamda;
  double NAG_HUGE *c;
  Integer init1;
  Integer init2;
} Nag_Spline;

typedef struct{
  Integer nx;
  double NAG_HUGE *lamda;
  Integer ny;
  double NAG_HUGE *mu; 
  double NAG_HUGE *c;
  Integer init1;
  Integer init2;
} Nag_2dSpline;





/* Typedefs for pointers to functions used by e04 functions  */
#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_E04ABC_FUN)(double, double NAG_HUGE *, Nag_Comm NAG_HUGE *);
typedef void (NAG_CALL * NAG_E04BBC_FUN)(double, double NAG_HUGE *, double NAG_HUGE *, Nag_Comm NAG_HUGE *);

typedef void (NAG_CALL * NAG_E04FCC_FUN)(Integer, Integer, 
	  double NAG_HUGE *, double NAG_HUGE *, Nag_Comm NAG_HUGE *);   

typedef void (NAG_CALL * NAG_E04GBC_FUN)(Integer, Integer,
	  double NAG_HUGE *, double NAG_HUGE *, double NAG_HUGE *, Integer, Nag_Comm NAG_HUGE *);   
typedef void (NAG_CALL * NAG_E04LBC_FUN)(Integer, double NAG_HUGE *, double NAG_HUGE *, 
					 double NAG_HUGE *, Nag_Comm NAG_HUGE *);
typedef void (NAG_CALL * NAG_E04LBC_HESS)(Integer, double NAG_HUGE *, double NAG_HUGE *, 
					  double NAG_HUGE *, Nag_Comm NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_E04ABC_FUN)();
typedef void (NAG_CALL * NAG_E04BBC_FUN)();
typedef void (NAG_CALL * NAG_E04FCC_FUN)();
typedef void (NAG_CALL * NAG_E04GBC_FUN)();
typedef void (NAG_CALL * NAG_E04LBC_FUN)();
typedef void (NAG_CALL * NAG_E04LBC_HESS)();
#endif
/* Structures used by e04 functions */

/* Nag_Fun contains function pointers for use by least squares routines,
   these pointers are assigned the objective function and the function
   which calls the objective function. This assignment is made in the top-level
   function e04fcc or e04gbc. Objective function may either supply
   derivatives (lsf_deriv) or not (lsf_noderiv). The type of function
   assigned is denoted by value of member type.
*/
typedef struct {
  Nag_FunType type;

#ifdef NAG_PROTO
  void (*call_lsf_noderiv)(Integer m, Integer n, NAG_E04FCC_FUN lsf,
			   double NAG_HUGE x[], double NAG_HUGE fvec[], double NAG_HUGE fjac[],
			   Integer tdj, Nag_Comm NAG_HUGE *comm);

  NAG_E04FCC_FUN lsf_noderiv;

  void (*call_lsf_deriv)(Integer m, Integer n, NAG_E04GBC_FUN lsf,
			 double NAG_HUGE x[], double NAG_HUGE fvec[], double NAG_HUGE fjac[],
			 Integer tdj,Nag_Comm NAG_HUGE *comm);

  NAG_E04GBC_FUN lsf_deriv;
#else
  void (*call_lsf_noderiv)();
  void (*lsf_noderiv)();
  void (*call_lsf_deriv)();
  void (*lsf_deriv)();
#endif

} Nag_Fun;

typedef struct
{
    Integer loclc[20]; /* workspace partition indices */
    Integer litotl;    /* total integer workspace in use */
    Integer lwtotl;    /* total real workspace in use */
    Integer lennam, ldt, ncolt, ldq;
} Nag_ae04mf;

typedef struct {
  Nag_ProblemType prob;
  Integer m;
} Nag_e04nfu;

typedef struct {
  double asize, dtmax, dtmin;
} Nag_de04nb;

typedef struct
{
    double tolx0, tolinc;
    Integer idegen, kdegen, ndegen, itnfix, nfix[2];
} Nag_ce04mf;

typedef struct
{
    double alfa, trulam;
    Integer isdel, jdel, jadd;
    Boolean header, prnt;
} Nag_de04mf;

typedef struct
{
    Integer idbg;
    Boolean lcdbg;
    Integer ilcdbg[5];
    Boolean cmdbg;
    Integer icmdbg[5];
} Nag_ee04mf;

typedef struct   /* e04xaz variables saved from a previous call */
{
    double cdsave, fdsave, hsave, oldh, rho, sdsave, ce1big, ce2big, te2big;
} Nag_XazSt;

typedef struct   /* e04ucj, e04uck variables saved from a previous call */
{
    Boolean braktd, crampd, extrap, moved, wset;
    Boolean vset; /* e04ucj only */
    Integer nsamea, nsameb;
    double a, b, xtry, factor, xw,fw, tolmax;
    double gw; /* e04uck only */
    double fa, xv,fv; /* e04ucj only */
} Nag_UcjkSt;

typedef struct   /* e04uc0, e04uc1 variables saved from a previous call, consequence of 
                    not updating e04dgc at Mark 4 */
{
    Boolean braktd, crampd, extrap, moved, wset;
    Boolean vset; /* e04ucj only */
    Integer nfsrch, nsamea, nsameb;
    double a, b, alfuzz, factor, xw,fw, tolmax;
    double rtmin, gw; /* e04uck only */
    double fa, xv,fv; /* e04ucj only */
} Nag_Uc01St;

typedef struct
{
   Integer lfdset; /* 0 = finite diff interval not supplied, e04xay calculates
		    * and puts in hforwd[] and hcntrl[].
		    * 1 = One interval for forwd. (fdint) &
		    *     one for central (cdint).
		    * 2 = forwd. & cent. intervals supplied bu user, held in
		    *     hforwd[] and hcntrl[]. Checked by e04xay */
   Integer ncdiff; /* number of missing gradient elements from obj. function */
   Integer nfdiff; /* number of missing Jacobian elements from con. function */
   Integer lvldif; /* Used in e04ucv, e04ucz, e04upz routines for e04ucc,
		    * e04upc.
		    * e04ucv not called by anything according to calldall.
		    * 0 means no finite differences wanted
		    * 1 means forward differences wanted
		    * 2 means central differences wanted  */
/* *******************************************************************
 *   double NAG_HUGE *hforward;
 *   double NAG_HUGE *hcentral;
 *   double fd_int;
 *   double cd_int;
 * ******* should these related variables be included in the future? */ 
} Nag_Deriv_Inf;

typedef struct
{
    Integer inpdbg[5];
    Boolean print; /* was debug */
} Nag_DebugSt;

typedef struct
{
    Nag_GradChk type;       /* was lvrfyc */
    Integer obj_start;  /* was jverfy[0] */
    Integer obj_stop;   /* was jverfy[1] */
    Integer con_start;  /* was jverfy[2] */
    Integer con_stop;   /* was jverfy[3] */
    int g_error;
} Nag_Grad_Chk_St;

typedef struct
{
    Integer n_elements; /* number of elements assigned by user in objfun() */
    double dir_deriv; /* directional derivative */
    double fd_approx; /* forward difference approximation */
    Boolean correct;
    double max_error;
    Integer max_subfunction; /* Subfunction index for max deviation between user & fd derivs */
    Integer max_constraint; /* Constraint index for max deviation between user & fd derivs */
} Nag_SimSt;

typedef struct
{
    Boolean assigned;
    double hopt;
    double gdiff;
    Integer iter;
    Boolean correct;
    char NAG_HUGE *comment;
    Boolean zero;
    double error;
} Nag_CompSt;

typedef struct
{
    double hfd;
    double hcd;
    double errmax;
} Nag_FDIntSt;

typedef struct
{
    Nag_Grad_Chk_St g_chk;
    Integer nfdiff;       /* number of missing non-constant g elements */
    Integer ncdiff;       /* number of missing non-constant J elements */
    Integer nf_const;     /* number of missing constant g elements */
    Integer nc_const;     /* number of missing constant J elements */
    Nag_DerivSet deriv_level;
    Nag_DerivSet old_deriv_level;
    Nag_SimSt   f_sim;        /* simple check of grad of f */
    Nag_SimSt   c_sim;        /* simple check of Jacobian of c */
    Nag_CompSt  NAG_HUGE *f_comp;       /* component check of grad of f */
    Nag_CompSt  NAG_HUGE *c_comp;       /* component check of Jacobian of c */
    Nag_FDIntSt NAG_HUGE *intervals;    /* fd intervals and error */
    Nag_DPrintType print_deriv; /* Print level for derivative checking info */
} Nag_GPrintSt;

typedef struct
{
    Boolean refactor;
    Integer jmax;
    double errmax;
    Boolean moved;
    Integer nmoved;
    Boolean feasible;
    Boolean rowerr;
    Boolean first;
    Boolean header;
    Boolean prnt;
} Nag_QP_Print;


typedef struct
{
    Integer m;  /* number of residuals */
    Integer n;  /* number of variables */
    Integer numinf; /* Number of infeasibilities */
    double NAG_HUGE *x;  /* current point x[] */
    double  f;  /* objective function value */ 
    double objf; /* Objective function value */
    double NAG_HUGE *g;  /* gradient of objfun at current point or estimate,
		   may be out of date for variables on bounds */

    double NAG_HUGE *diagt; /* diagonal elements of t */
    double NAG_HUGE *fvec; /* value of residuals */
    double NAG_HUGE *fjac; /* matrix of first derivatives at current point */
    double NAG_HUGE *conjac; /* matrix of first derivatives at non-linear constraint functions */
    Integer tdj;  /* trailing dimension of fjac[] */
    Integer tdconjac;  /* trailing dimension of conjac[] */
    double step;  /* step length of line search */
    double xk_norm;   /* distance moved in last iteration */
    double gpj_norm;  /* norm of projected gradient (g of free variables) */
    double cond;   /* estimate of condition number of proj. Hessian */
    double cond_h;   /* lower bound on the condition number of
			proj. Hessian h*/
    double cond_hz;   /* lower bound on the condition number of
			proj. Hessian hz */
    double cond_t;   /* lower bound on the condition number of
			the matrix of predicted active constraints */
    double cond_r;    /* Lower bound on condition number of Rz */

    Boolean iter_conv;     /* Status of iterates converged convergence test */
    Boolean norm_gz_small; /* Status of size of projected gradient convergence test */
    Boolean violtn_small;  /* Status of size of constraint violations convergence test */

    Boolean update_modified;  /* TRUE if q-Newton update modified to ensure pos. def. */
    Boolean c_diff;           /* TRUE if central diff approxns made to unspecified grads */
    Boolean qp_not_feasible;  /* TRUE if the QP subproblem was infeasible */
    Boolean step_limit_exceeded;  /* TRUE if the linesearch step limit was exceeded */
    Boolean refactor;         /* TRUE if the approximate Hessian has been refactorised */

    char NAG_HUGE *mgrmsg; /* Contains the codes M,I,L,C,R */

    double norm_nlnviol; /* Norm of non-linear constraint violations */

    Integer grade; /* grade of Jacobian matrix */
    double NAG_HUGE *s;     /* singular values of Jacobian */
    Integer NAG_HUGE *state; /* inf. on which variable are free or on bounds */
    Integer iter;   /* Number of iterations - major or minor depending on context */
    Integer major_iter; /* Number of major iterations at end of e04ucc */
    Integer minor_iter; /* Minor iteration count at end of e04ucc major iteration */
    Integer nfun; /* Objective function evaluation for line search */
    double merit; /* augmented Lagrangian merit function */
    double violtn; /* Euclidean norm of the residuals of violated
	       or predicted active set */
    
    double penalty; /* Euclidean norm of the vector of penalty
		       parameters */
    
    Integer itn; /* number of major iters/10000 (e04ncj) */
    Integer nf;    /* number of function evaluations not including those
		      for derivative testing or for estimating fd intervals */
    Integer nctotal; /* n + lin constraints + nolin constraints */


    double NAG_HUGE *lncons_val; /* e04nbx lin constraints' value in the current working set */

    Boolean major_iter_soln; /* True if soln is for a major iteration */
    Boolean local_search; /* A local search has just been performed */

    Integer nclin;   /* Number of linear constraints */
    Integer ncnlin;   /* Number of non-linear constraints */
    Integer nactiv;  /* Number of active constraints (dimension of diagt) */
    Integer jdel;    /* Index of constraint deleted */
    Integer jadd;    /* Index of constraint added */
    Integer ninf;    /* Number of infeasiblities */
    double sum_infeasibilities; /* Sum of infeasibilities */
    Integer bnd;     /* Number of bound constraints in working set */
    Integer lin;     /* Number of linear constraints in working set */
    Integer nln;     /* Number of non-linear constraints in working set */
    double NAG_HUGE *t;       /* Upper triangular matrix T, */
    Integer tdt;     /* Trailing dimension for T */
    double NAG_HUGE *r;       /* Upper triangular matrix R, */
    double NAG_HUGE *diagr;   /* Diagonal of upper triangular matrix R, */
    Integer tdr;     /* Trailing dimension for r */
    Integer nart;    /* Number of artificial constraints in working set */
    Integer nartm;    /* Multipliers for the artificial constraints in
			 working set (i.e., nz-nrz-1)*/
    Integer nrz;     /* Number of columns of matrices Zr and R */
    Integer nz;     /* Number of columns of matrices Z */
    double norm_gz;  /* Norm of reduced gradient Gz */
    double norm_gf;  /* Norm of gradient w.r.t free variables */
    Integer nopt;     /* Number of non-optimal multipliers */
    double min_lm;   /* Value of multiplier associated with deleted constraint. */
    double condt;    /* Lower bound on condition number of T */
    double condr;    /* Lower bound on condition number of Rz */
    double rzz;      /* Last diagonal element of D */
    Boolean rset;    /* TRUE if r, condr and rzz are set */
    double NAG_HUGE *ax;      /* Current values of linear constraints Ax */
    double NAG_HUGE *cx;      /* Current values of non-linear constraints C */
    Integer nrank; /* e04ucc Required to allow diag elements of r to be printed */
    Integer nfixed; /* e04ucc Required to allow bound constraints to be printed */

    /* The following members relate to Lagrange multipliers */

    Integer NAG_HUGE *kx; /* Indices of bound constraints with associated multipiers.
		  * kx[i] is the index of the variable which is active on 
		  * a bound.  bclambda[i] is the associated value; 
		  * i = 0, 1, ..., bnd -1
		  */

    Integer NAG_HUGE *kactive; /* Inverted kactive[], indices correspond
		       * directly to lambda[]. E.g value of kactive[i]
		       * is the index of the constraint with
		       * multiplier lambda[i], i = 0, 1, ..., lin-1
		       */
    double NAG_HUGE *lambda;   /* In e04ncj (QP iteration printing), lambda is
		       * the array of multipliers associated with the
		       * active linear constraints, i = 0,1,...,lin-1.
		       *
		       * For e04nbx (solution printing) lambda contains
		       * the lagrange multipliers for bounds constraints,
		       * linear constraints, and nonlinear constraints
		       * (if any), in that order.  Indexing is by
		       * variable and constraint index (so zeros are
		       * given for LMs corresp. to non-active constraints.
		       */
    double NAG_HUGE *bclambda; /* Multipliers for the bound constraints 0, ..., bnd-1 */
    double NAG_HUGE *gq;       /* Multipliers of artificial constraints 0, ..., nart-1*/
   
   /* These are same as Nag_QP_Print */
    Integer jmax;
    double errmax;
    Boolean moved;
    Integer nmoved;
    Boolean feasible;
    Boolean rowerr;
    Boolean first;

    double NAG_HUGE *bl;
    double NAG_HUGE *bu;
    Nag_EndState endstate;
    Nag_GPrintSt NAG_HUGE *gprint; /* Substructure containing results of
			   * e04ucy derivative check */
/* Extra members for the simplex case */
    double fmin;
    double fmax;
    double NAG_HUGE *simplex; /*(n+1)*n array containig the simplex vertices */
/* Extra members for e04ucc */
    char NAG_HUGE *prbtyp;
/* Extra members for e04unc */
    Integer tdfjac;  /* Trailing dimension of fjac */

/* Members for e04nkc */

    /* -- Solution -- */
    Boolean col;    /* TRUE if the solution data refers to a column (variable)
		     * rather than a row (constraint).
		     */
    Integer index;  /* Index of row or column (1..m, 1..n resp.) */
    char *name;     /* Row/col name */
    char *key;      /* "D" = degen bas/sbas  "I" = infeas b/sb
		     * "A" = altern opt (degen nonb dual)
		     * "N" = nonopt sb/nonb (infeas dual)
		     */
    char *sstate;  /* State of variable/constraint */
    double val;     /* Value of variable/constraint activity */
    double blo;     /* Lower bound of var/constr */
    double bup;     /* Upper "       "           */
    double lmult;   /* Lagrange mult of bnd/constraint */
    /* double resid;   *//* Residual/slack activity of variable /constraint */
    double objg;      /* Objective gradient */ 

    /* -- Iteration -- */
    Boolean just_feas;  /* True if problem has just become feasible */
    Boolean qp;         /* True if qp problem                       */
    Integer pprice;     /* Partial price indicator                  */
    double rgval;       /* Red grad/cost of var in pricing opern.   */
    Integer sb_add;     /* Variable selected to enter super bas set */
    Integer sb_leave;   /* Variable leaving superbasic set          */
    Integer b_leave;    /* Variable leaving basis -> nonbasic       */
    Integer bswap_leave;/* Variable leaving basis in bas<->sbas swap*/
    double pivot;       /* Value of pivot                           */
    Integer nnz_l;      /* Number of nnz in lower basis factor      */
    Integer nnz_u;      /* Number of nnz in upper basis factor      */
    Integer ncp;        /* Number of compressions                   */
    double norm_rg;     /* Norm of reduced gradient                 */
    Integer nsb;        /* Number of superbasics                    */

    /* Members for H02 (MIP) */
    Boolean root;
    Integer nnodes;
    Integer node_num;
    Integer parent_node;
    char **crnames;
    Nag_NodeStatus node_status;
    Integer branch_index;
    double x_before;
    double x_after;
    double x_lo;
    double x_up;
    Integer depth;

    /* e04lbc */
    Boolean posdef;
    
} Nag_Search_State;

typedef struct {        /* Used by e04xyc, e04xyz when reading options from file */
    Boolean options;
    Boolean begin;
    Boolean end;
    Boolean eof;
} Nag_Opt_Found;



/* Typedef for pointers to user defined e04 print function */
#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_E04_PRINTFUN)(const Nag_Search_State NAG_HUGE *,
			    Nag_Comm NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_E04_PRINTFUN)();
#endif

typedef struct {
char NAG_HUGE *f_name;           /* name of the function for which the structure
			   has been initialised  */
char NAG_HUGE *f_longname;       /* Nag long name of function for which structure
			 * has been initialised */
Integer init1;          /* flags indicating structure has been initialised, */
Integer init2;          /* set to IDUMMY and INIT2DUMMY = -23456 on initialisation */
Nag_LinFun minlin;      /* e04gbc, e04kbc: enum specifies linear minimisation routine.*/
Integer grade;          /* e04fcc, e04gbc: output, grade of Jacobian */
Integer nf;             /* output: number of function evaluations or number of times
			   residuals have been calculated. Calls to objective function */
Integer iter;           /* e04dgc, e04fcc, e04gbc, e04jbc, e04ucc: number of
			   iterations performed by routine */
Integer max_iter;	/* iteration limit. */
Integer minor_max_iter;	/* minor iteration limit. */
Boolean min_infeas;     /* e04mfc, e04nfc: minimise sum of infeasibilities
			   when TRUE and problem infeasible. Default = FALSE */
double f_est;           /* e04jbc, e04kbc, e04dgc: estimate of function value at
			   minimum */
double step_max;        /* e04fcc, e04gbc, e04jbc, e04kbc: Estimate of distance to
			   solution on entry. default = 100000.0 */
double step_limit; 	/* e04ucc: maximum change at first step of line search.
                           Default = 2.0 */
double max_line_step;   /* e04dgc: maximum step length, default = 10**10 */
double f_prec;          /* e04dgc, e04ucc: function precision, default = eps**0.9 */
double linesearch_tol;  /* e04dgc, e04ucc: linesearch tolerance, default = 0.9
			   e04fcc, e04gbc, e04jbc, e04kbc: accuracy of linear
			   minimisation, eta. Default = 0.5 but if n = 1
			   default = 0.0. */
double optim_tol;       /* e04dgc, e04ucc: optimality tolerance.
			   Default = r_eps**0.8.
			   e04fcc, e04gbc, e04jb, e04kbcf: accuracy
			   of solution, xtol. 
			   Default = sqrt(eps) or 0.0. */
double NAG_HUGE *conf;		/* e04ucc */
Integer nag_conf;
double NAG_HUGE *conjac;		/* e04ucc */
Integer nag_conjac;
/*double NAG_HUGE *r;		 e04ucc */
/* Integer nag_r; */
double NAG_HUGE *h;		/* e04ucc */
Integer nag_h;

Boolean h_unit_init;    /* e04unc */
Integer h_reset_freq;   /* e04unc */

double NAG_HUGE *s;              /* e04fcc, e04gbc: singular values of Jacobian,
			   minimum size n.*/
Integer nag_s;
double NAG_HUGE *v;              /* e04fcc, e04gbc: matrix V of SVD, min size n*n. */ 
Integer nag_v;
Integer tdv;            /* Trailing dimension for array v */
double NAG_HUGE *delta;          /* e04jbc: Difference intervals deriv. estimation.
			   minimum size n.*/
Integer nag_delta;
double NAG_HUGE *hesl;           /* e04jbc, e04kbc: factor L minimum size max(n*(n-1)/2, 1) */
Integer nag_hesl;
double NAG_HUGE *hesd;           /* e04jbc, e04kbc: factor D minimum size n */
Integer nag_hesd;
Integer NAG_HUGE *state;         /* e04jbc, e04kbc, e04nfc, e04mfc, e04ucc: */
Integer nag_state;      /*  which var. or con. on bounds, size n + nclin + ncnln */
Nag_InitType init_state;        /* e04jbc, e04kbc: enum, which arrays initialised */
Boolean obj_deriv;      /* e04ucc: If true, all objective gradient elements supplied */
Boolean con_deriv;      /* e04ucc: If true, all constraint gradients supplied */
double f_diff_int;      /* e04ucc: forward difference interval */
double c_diff_int;      /* e04ucc: central difference interval */
Integer obj_check_start;/* e04dgc, e04ucc: element to start checking obj. gradients,
		   only if verify_grad Nag_CheckObj or Nag_CheckObjCon. Default = 1 */
Integer obj_check_stop; /* e04dgc, e04ucc: element to stop checking obj. gradients,
		   only if verify_grad Nag_CheckObj or Nag_CheckObjCon. Default = n */
Integer con_check_start;/* e04dgc, e04ucc: element to start checking con. gradients,
		   only if verify_grad Nag_CheckCon or Nag_CheckObjCon. Default = 1 */
Integer con_check_stop; /* e04dgc, e04ucc: element to stop checking con. gradients,
		   only if verify_grad Nag_CheckCon or Nag_CheckObjCon. Default = n */
Nag_GradChk verify_grad;  /* e04dgc, e04ucc: type of check on user gradients,
			   default Nag_SimpleCheck */
Boolean hessian;        /* e04ucc: return the hessian. Default = FALSE */
Boolean local_search;   /* e04jbc, e04kbc: perform local search, default TRUE */
Boolean deriv_check;    /* e04gbc: Check derivatives using e04yac if TRUE,
				   default = TRUE. */
Boolean list;           /* e04dgc, e04ucc: list options, default = TRUE */
Boolean print_gcheck;	/* e04dgc: If TRUE print grad. check results. Default = TRUE */
Nag_DPrintType print_deriv; /* e04dgc: Print grad. check results
				Nag__D_Full */
Nag_PrintType print_level;      /* e04dgc, e04ucc, e04fcc, e04jbc, e04kbc:
				type of solution printout from
				major iterations */
Nag_PrintType minor_print_level; /* print level for QP solver called
				  * from e04ucc or e04upc. */
Integer print_iter;     /* e04gbc, e04fcc: print_iter will be replaced by
					   print_level.
					   print every print_iter iterations
					   print_iter = 0, last iter only,
					   print_iter = -1 no printout. */
char outfile[NAG_FILE_LEN];     /* file to which output should be given */
#ifdef NAG_PROTO
  NAG_E04_PRINTFUN print_fun;  /* In: user defined monitoring function */ 
#else
  void (*print_fun)();  /* In: user defined monitoring function */ 
#endif

double NAG_HUGE *ax;          /* output array pointer, values of linear constraints Ax */
Integer nag_ax;
double NAG_HUGE *lambda;     /* output array pointer, values of Lagrange multipliers */
Integer nag_lambda;
Integer fcheck;      /* input, check frequency */
double crash_tol;    /* input, crash tolerance */
Integer reset_ftol;  /* input, expand frequency, reset feasibility tolerance
		      * every reset_ftol iterations */
Integer fmax_iter;   /* input, feasibility phase iteration limit */
double ftol;         /* input, feasibility tolerance */
double lin_feas_tol; /* linear feasibility tolerance */
double nonlin_feas_tol; /* nonlinear feasibility tolerance */
Integer hrows;       /* input, number of rows of Hessian */
double inf_bound;    /* input, infinite bound size */
double inf_step;     /* input, infinite step size */
Integer max_df;         /* input, maximum degrees of freedom */
Nag_ProblemType prob;    /* input, enum specifies problem type */
double rank_tol;     /* input, condition number of R */
Boolean form_hessian;  /* TRUE if wish to recover Hessian factor e04ncc/ucc */
Nag_Start start;         /* input, Cold or Warm start */
Integer debug_iter;       /* Start debug output at this iteration */
Integer minor_debug_iter; /* Start QP debug output at this iteration */
Boolean debug;            /* Debug output required */
Boolean used;           /* TRUE if options structure reused */
Boolean unitq;       /* save, used for warm start in e04nfc */
Integer nfree;       /* save, used for warm start in e04nfc */
Integer nactiv;      /* save, used for warm start in e04nfc */

/* These are used in the e04ucc stringent */
Integer nprob;
Integer n;
Integer nclin;
Integer ncnlin;
double bndlow;
double bndupp;

/* New e04nkc members */
Boolean minimize;      /* Dirn. of optn. (min or max)     */
Integer factor_freq;   /* Basis refactorisation frequency */
Integer partial_price; /* Relevant to LP problems         */
Integer max_sb;        /* Maximum number of superbasics   */
double lu_factor_tol;  /* LU factorisation tolerance      */
double lu_sing_tol;    /* LU singularity tolerance        */
double lu_update_tol;  /* LU update tolerance             */
double pivot_tol;      /* Pivot selection criterion       */
double scale_tol;      /* Scale tolerance                 */
Nag_CrashType crash;   /* Crash option                    */
Nag_ScaleType scale;   /* Scale option                    */

char prob_name[MPS_NAME_LEN+1];
char obj_name[MPS_NAME_LEN+1];
char bnd_name[MPS_NAME_LEN+1];
char rhs_name[MPS_NAME_LEN+1];
char range_name[MPS_NAME_LEN+1];

char **crnames;          /* Array of column and row names   */
Integer nag_crnames;

Integer nsb;           /* Output - number of superbasics  */

/* e04nkc undocumented options */
Integer max_restart;   /* Maximum no. of restarts on exhausting wkspace   */
Integer max_basis_len; /* User assigned workspace length                  */
Integer max_compress;  /* Maximum number of compressions before a restart */
Integer max_basis_nfactor;  /* No. basis refactors before assumed singlr. */

/* e04mzc options (MPS Reading) */
double col_lo_default;
double col_up_default;
double infinity;

Integer ncol_approx;
Integer nrow_approx;

double est_density;

Nag_OutputType output_level;

/* e04xac options */
Nag_DWantType deriv_want;

Boolean use_hfwd_init;
double f_prec_used; /* Output */

} Nag_E04_Opt;


/* Structure to handle variables which are SAVEd in the Fortran
 * version of sparse QP (avoiding use of static locals)
 */

typedef struct {
  Integer itnfix; /* Used in e04nlq */
  Integer nfix[2];   /* Used in e04nlq */
  double tolx0;   /* Used in e04nlq */
  Integer nbfac;  /* Used in e04nmj */
  double umin;    /* Used in e04nmj */
} Nag_SparseQP_Save_Vars;


/* New e04 temporary types */
typedef struct {
  Integer locls[20];
} Nag_ae04nc;

typedef struct {
  Integer locnp[35];
} Nag_ae04uc;

typedef struct {
  Integer locnl[20];
} Nag_ae04up;

typedef struct {
  Integer tdt, ncolt, tdq;
} Nag_be04nb;

typedef struct {
  Integer ilsdbg[5];
  Boolean lsdbg;
} Nag_ce04nc;

typedef struct {
  double rhomax, rhonrm, rhodmp, scale;
  Boolean incrun;
} Nag_de04uc;

typedef struct {
  double rcndbd, rfrobn, drmax, drmin;
} Nag_ee04nb;

typedef struct {
  Integer icmdbg[5];
  Boolean cmdbg;
} Nag_fe04nb;

typedef struct {
  Integer nactiv, nfree, nz;
  Boolean unitq;
} Nag_fe04nc;

typedef struct {
  Integer inpdbg[5];
  Boolean npdbg;
} Nag_fe04uc;

/* Structcures used by g02 functions */

typedef struct
{
  Integer max_iter;
  double  tol, eps;
  Integer print_val;
  Integer init1, init2;
} Nag_GlmAcc;

typedef struct
{
  Integer rank;
  double  df;
  double  rss, dev;
  double  NAG_HUGE *se_beta;
  double  NAG_HUGE *cov;
} Nag_GlmRes;

typedef struct
{
  Integer max_iter;
  Integer iter_mon;
  double  cpsi;
  double  cucv;
  double  dchi;
  double  hpsi[3];
  double  tol;
  Integer init1;
  Integer init2;
} Nag_RegData;

typedef struct
{
  double  NAG_HUGE *rs;
  double  NAG_HUGE *wgt;
  double  NAG_HUGE *work;
} Nag_ResWrk;

typedef struct
{
  Integer max_iter;
  Integer indm;
  Integer iter_mon;
  double  bl;
  double  bd;
  double  tol;
  Integer init1;
  Integer init2;
} Nag_CorrData;



/* Structcures used by g13 functions */
	
typedef struct
{
  Integer itc;
  double  rss; 
  double objf; 
  double NAG_HUGE *para; 
  double NAG_HUGE *sd; 
  double df;
  Integer npara; 
  Integer npe; 
  Integer NAG_HUGE *mtyp; 
  Integer NAG_HUGE *mser;
} Nag_UserPrintFun;

/* Typedef for pointers to user defined e04 print function */
#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_G13_PRINTFUN)(const Nag_UserPrintFun NAG_HUGE *,
			    Nag_Comm NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_G13_PRINTFUN)();
#endif
typedef struct {
  Boolean used;         /* TRUE if options structure reused */
  char NAG_HUGE *f_name;
  Integer init1;
  Integer init2;
  Boolean cfixed;          /* integer => Boolean */
  Nag_Likelihood criteria;     /* integer => enum */
  Integer max_iter;
  double alpha;
  double beta;
  double delta;
  double gamma;
  Integer iter;
  double NAG_HUGE *cm;
  Integer tdcm;
  Integer nag_cm;
  double NAG_HUGE *res;
  Integer lenres;
  Integer nag_res;
  double NAG_HUGE *zt;
  Integer tdzt;
  Integer nag_zt;
  double NAG_HUGE *noise;
  Integer nag_noise;
  Integer isttf;
  double NAG_HUGE *sttf;
  Integer nag_sttf;
  Integer nsttf;
  Nag_PrintType print_level;    /* g13bec
				type of solution printout from
				major iterations */
  Boolean list; /* List Option structure */
#ifdef NAG_PROTO
  NAG_G13_PRINTFUN print_fun; /* In: User print function */
#else
  void (*print_fun)();          /* In: User print function */
#endif
  char outfile[NAG_FILE_LEN];
} Nag_G13_Opt;


typedef struct {
  Integer p;
  Integer d;
  Integer q;
  Integer bigp;
  Integer bigd;
  Integer bigq;
  Integer s;
} Nag_ArimaOrder;


typedef struct {
  Integer init3;
  Integer NAG_HUGE *b;
  Integer nag_b;
  Integer NAG_HUGE *q;
  Integer nag_q;
  Integer NAG_HUGE *p;
  Integer nag_p;
  Integer NAG_HUGE *r;
  Integer nag_r;
} Nag_TransfOrder;


/* Structures used by h02 chapter */

/* Typedef for pointers to user defined h02 print function */
#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_H02_PRINTFUN)(const Nag_Search_State NAG_HUGE *,
			    Nag_Comm NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_H02_PRINTFUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_H02_QPFUN)(Integer, Integer,
					double NAG_HUGE *, Integer, double NAG_HUGE *,
					double NAG_HUGE *, Nag_Comm NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_H02_QPFUN)();
#endif



  /* Options for MPS reading */
typedef struct
{
  Integer init1;
  Integer init2;

  char prob_name[MPS_NAME_LEN+1];
  char bnd_name[MPS_NAME_LEN+1];
  char obj_name[MPS_NAME_LEN+1];
  char range_name[MPS_NAME_LEN+1];
  char rhs_name[MPS_NAME_LEN+1];

  double col_lo_default;
  double col_up_default;
  double infinity;

  Integer ncol_approx;
  Integer nrow_approx;

  Boolean sparse;
  double est_density;

  Boolean requ_ivar;

  Boolean list;
  Nag_OutputType output_level;
  char outfile[NAG_FILE_LEN];
}
Nag_MPS_Opt;

typedef struct
{
  /* Name/type */
  char *f_name;
  char *f_longname;
  Nag_MIP_ProbType prob;
  
  /* 'Hidden stuff' */
  Integer init1;  /* Initialisation checks */
  Integer init2;
  Boolean used;
  
  /* User defined printing etc. */
  Boolean list;
  Nag_PrintType print_level;
  NAG_H02_PRINTFUN print_fun;
  char outfile[NAG_FILE_LEN];

  /* Algorithm parameters */
  /* -- LP */
  Integer max_iter;   /* Subproblem iteration limit */
  double inf_bound;  
  double feas_tol;    /* Feasibility tolerance for constraints */

  /* -- QP */
  Integer hrows;
  Integer max_df;
  double rank_tol;

  /* -- IP */
  Boolean first_soln;
  Integer max_depth;  /* Max depth of BB tree */
  Integer max_nodes;  /* Max number of nodes to be solved */
  double int_obj_bound; 
                      /* Upper bound on objective at integer solution */
  double soln_tol; /* Tolerance using to decide to fathom a node by bounds  */
  double int_tol;  /* Max deviation of var from integral value for int soln.*/
  Nag_Node_Selection nodsel;        /* Node selection strategy     */
  Nag_Var_Selection varsel;    /* Variable selection strategy */
  Nag_Branch_Direction branch_dir; 
  double *priority; /* Pointer to array of priorities for branching */

  /* Solution details */
  double *lower;
  double *upper;
  double *lambda;
  Integer *state;

  Integer nag_lower;
  Integer nag_upper;
  Integer nag_state;
  Integer nag_lambda;

  /* MPS & row/col names */
  char prob_name[MPS_NAME_LEN+1];
  char obj_name[MPS_NAME_LEN+1];
  char bnd_name[MPS_NAME_LEN+1];
  char rhs_name[MPS_NAME_LEN+1];
  char range_name[MPS_NAME_LEN+1];
  
  char **crnames;          /* Array of column and row names   */
  Integer nag_crnames;

  /* MPS Reading options (h02buc) */
  double col_lo_default;
  double col_up_default;
  double infinity;

  Integer ncol_approx;
  Integer nrow_approx;

  Nag_OutputType output_level;

  Integer n_ivar; /* Output: number of integer variables in MPS file */

  /* Undocumented debug options */
  Boolean dbg;
  Boolean dbg_rt;
  Boolean dbg_rqip;
  Nag_EndState dbg_es;

} Nag_H02_Opt;

typedef struct 
{
  Integer init1;     /* Set to magic number by initialisation function */
  Integer init2;
  Boolean used;

  Nag_MIP_ProbType type;
  
  Integer n_ivar;     /* Number of integer variables                      */
  Integer *ivar;   /* Array of integer variable indices                */

  Integer n;           /* Number of variables   */
  Integer m;           /* Number of constraints */

  double *bl;    /* Lower bounds on variables & constraints */
  double *bu;    /* Upper bounds on variables & constraints */

  char **names;     /* Names of variables and constraints */
   
  double *c;        /* Linear objective term    */
  double *h;        /* Hessian term             */
  Integer tdh;      

  NAG_H02_QPFUN qphess;

  double *a;        /* Linear constraints */
  Integer tda;      

  /* Names of bounds etc. */
  char prob_name[MPS_NAME_LEN+1];
  char bnd_name[MPS_NAME_LEN+1];
  char obj_name[MPS_NAME_LEN+1];
  char range_name[MPS_NAME_LEN+1];
  char rhs_name[MPS_NAME_LEN+1];

  /* Additional data for sparse problem */
  Integer nnz;
  Integer iobj;
  Integer *ha;
  Integer *ka;

} Nag_MIP_Problem;

typedef struct 
{
  double objval;
  double *x;
  double *lambda;
  Integer *state;
  Nag_EndState endstate;
} Nag_MIP_Solution;

typedef struct
{
  char *prbtyp;
  Boolean vertex;
  Boolean cset;
  Integer *iw;
  double *w;
  Nag_ae04mf *ae04mf;
  Nag_ee04mf *ee04mf;
  Nag_E04_Opt *e04_opt;
  Nag_Comm *e04_comm;
  Nag_FileSt *e04_stream;
} Nag_e04nfg_Args;

typedef struct _s_Nag_HashWrap Nag_HashWrap;
struct _s_Nag_HashWrap
{
  Pointer data;       /* The user's data */
  char *name; 
  Nag_HashWrap *next;
};

typedef struct
{
  Nag_HashWrap *pool;
  Nag_HashWrap *next_free;
  unsigned long pool_size;
  unsigned long n_used;
  
  Nag_HashWrap **grow;
  unsigned long total_pool_size;
  Integer ngrow;
  Integer maxgrow;
} Nag_HashPool;

typedef struct
{
  Nag_HashWrap *hw_array;  /* Allocate array size table_size */
  Nag_HashPool hpool;
  unsigned long table_size;       /* Ideally, a prime number */
} Nag_HashTable;


#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_C05ADC_FUN)(double);
#else
typedef double (NAG_CALL * NAG_C05ADC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_C05NBC_FUN)(Integer, double NAG_HUGE *,
					  double NAG_HUGE *, Integer NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_C05NBC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_C05PBC_FUN)(Integer, double NAG_HUGE *, double NAG_HUGE *,
					  double NAG_HUGE *, Integer, Integer NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_C05PBC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_C05ZBC_FUN)(Integer, double NAG_HUGE *, double NAG_HUGE *,
					  double NAG_HUGE *, Integer, Integer NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_C05ZBC_FUN)();
#endif

/* Multi-threading versions of the c05s above */
#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_C05SDC_FUN)(double, Nag_User NAG_HUGE *);
#else
typedef double (NAG_CALL * NAG_C05SDC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_C05TBC_FUN)(Integer, double NAG_HUGE *,
					  double NAG_HUGE *, Integer NAG_HUGE *, Nag_User NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_C05TBC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_C05UBC_FUN)(Integer, double NAG_HUGE *, double NAG_HUGE *,
	  double NAG_HUGE *, Integer, Integer NAG_HUGE *, Nag_User NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_C05UBC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_C05ZCC_FUN)(Integer, double NAG_HUGE *, double NAG_HUGE *,
	  double NAG_HUGE *, Integer, Integer NAG_HUGE *,Nag_User NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_C05ZCC_FUN)();
#endif





/* Typedefs for pointers to functions used by e04's */
#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_E04CCC_FUN)(Integer n, double NAG_HUGE *, 
					double NAG_HUGE *, Nag_Comm NAG_HUGE *);
typedef void (NAG_CALL * NAG_E04UCC_FUN)(Integer, double NAG_HUGE *, double NAG_HUGE *,
					double NAG_HUGE *, Nag_Comm NAG_HUGE *);   
typedef void (NAG_CALL * NAG_E04UNC_OBJFUN)(Integer, Integer, double NAG_HUGE *, double NAG_HUGE *,
					   double NAG_HUGE *, Integer, Nag_Comm NAG_HUGE *);

typedef void (NAG_CALL * NAG_E04UCC_CONFUN)(Integer, Integer, Integer NAG_HUGE *,
					   double NAG_HUGE *, double NAG_HUGE *, double NAG_HUGE *,
					   Nag_Comm NAG_HUGE *);
typedef void (NAG_CALL * NAG_E04UNC_CONFUN)(Integer, Integer, Integer NAG_HUGE *,
					   double NAG_HUGE *, double NAG_HUGE *, double NAG_HUGE *,
					   Nag_Comm NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_E04CCC_FUN)();   
typedef void (NAG_CALL * NAG_E04UCC_FUN)();   
typedef void (NAG_CALL * NAG_E04UCC_CONFUN)();
typedef void (NAG_CALL * NAG_E04UNC_CONFUN)();
#endif

typedef NAG_E04UCC_FUN NAG_E04DGC_FUN;
typedef NAG_E04GBC_FUN NAG_E04YAC_FUN;


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_E04JBC_FUN)(Integer, double NAG_HUGE *, double NAG_HUGE *,
					double NAG_HUGE *, Nag_Comm NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_E04JBC_FUN)();
#endif   
typedef NAG_E04JBC_FUN NAG_E04HBC_FUN;
typedef NAG_E04JBC_FUN NAG_E04HCC_FUN;
typedef NAG_E04JBC_FUN NAG_E04KBC_FUN;


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_E04NFC_FUN)(Integer, Integer,
					double NAG_HUGE *, Integer, double NAG_HUGE *,
					double NAG_HUGE *, Nag_Comm NAG_HUGE *);
#else
typedef void (NAG_CALL * NAG_E04NFC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_E04NKC_HESSFUN)
	(Integer, double[], double[], Nag_Comm *);
#else
typedef void (NAG_CALL * NAG_E04NKC_HESSFUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_E04NKC_HESSFUN_WRAP)
	(NAG_E04NKC_HESSFUN, Integer, Integer, Integer, double[],
	Integer[], Integer[], double[], double[], Integer,
	Integer[], Integer, double[], Integer, Integer[],
	double[], Nag_Comm *);
#else
typedef void (NAG_CALL * NAG_E04NKC_HESSFUN_WRAP)();
#endif

#ifdef NAG_PROTO
typedef Integer (NAG_CALL * NAG_M01_FUN)(const Pointer, const Pointer);
#else
typedef Integer (NAG_CALL * NAG_M01_FUN)();
#endif


#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_D01_FUN)(double);
#else
typedef double (NAG_CALL * NAG_D01_FUN)();
#endif

typedef NAG_D01_FUN NAG_D01AJC_FUN;
typedef NAG_D01_FUN NAG_D01AKC_FUN;
typedef NAG_D01_FUN NAG_D01ALC_FUN;
typedef NAG_D01_FUN NAG_D01AMC_FUN;
typedef NAG_D01_FUN NAG_D01ANC_FUN;
typedef NAG_D01_FUN NAG_D01APC_FUN;
typedef NAG_D01_FUN NAG_D01AQC_FUN;
typedef NAG_D01_FUN NAG_D01ASC_FUN;
typedef NAG_D01_FUN NAG_D01BAC_FUN;

#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_D01GBC_FUN)(Integer, double NAG_HUGE *);
#else
typedef double (NAG_CALL * NAG_D01GBC_FUN)();
#endif

#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_D01FCC_FUN)(Integer, double NAG_HUGE *);
#else
typedef double (NAG_CALL * NAG_D01FCC_FUN)();
#endif

/* Multi-threading versions of the d01s above */

#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_D01_TS_FUN)(double, Nag_User NAG_HUGE *comm);
#else
typedef double (NAG_CALL * NAG_D01_TS_FUN)();
#endif


typedef NAG_D01_TS_FUN NAG_D01SJC_FUN;
typedef NAG_D01_TS_FUN NAG_D01SKC_FUN;
typedef NAG_D01_TS_FUN NAG_D01SLC_FUN;
typedef NAG_D01_TS_FUN NAG_D01SMC_FUN;
typedef NAG_D01_TS_FUN NAG_D01SNC_FUN;
typedef NAG_D01_TS_FUN NAG_D01SPC_FUN;
typedef NAG_D01_TS_FUN NAG_D01SQC_FUN;
typedef NAG_D01_TS_FUN NAG_D01SSC_FUN;
typedef NAG_D01_TS_FUN NAG_D01TAC_FUN;

#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_D01XBC_FUN)(Integer, double NAG_HUGE *, Nag_User NAG_HUGE *comm);
#else
typedef double (NAG_CALL * NAG_D01XBC_FUN)();
#endif

#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_D01WCC_FUN)(Integer, double NAG_HUGE *, Nag_User NAG_HUGE *comm);
#else
typedef double (NAG_CALL * NAG_D01WCC_FUN)();
#endif





#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02CJC_FUN)(Integer neq, double x, double NAG_HUGE y[],
					double NAG_HUGE f[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02CJC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02CJC_OUTFUN)(Integer neq, double NAG_HUGE *xsol,
					   double NAG_HUGE y[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02CJC_OUTFUN)();
#endif

#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_D02CJC_GFUN)(Integer neq, double x, double NAG_HUGE y[],
					   Nag_User NAG_HUGE *comm);
#else
typedef double (NAG_CALL * NAG_D02CJC_GFUN)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02EJC_FUN)(Integer neq, double x, double NAG_HUGE y[],
					double NAG_HUGE *f, Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02EJC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02EJC_PFUN)(Integer neq, double x, double NAG_HUGE y[], 
					 double NAG_HUGE pw[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02EJC_PFUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02EJC_OUTFUN)(Integer neq, double NAG_HUGE *xsol,
					   double NAG_HUGE y[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02EJC_OUTFUN)();
#endif

#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_D02EJC_GFUN)(Integer neq, double x, double NAG_HUGE y[],
					   Nag_User NAG_HUGE *comm);
#else
typedef double (NAG_CALL * NAG_D02EJC_GFUN)();
#endif



#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAZ_fcn)(Integer neq, double x, double NAG_HUGE y[],
		       double NAG_HUGE f[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAZ_fcn)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAZ_g)(Integer neq, double eps, double NAG_HUGE ya[], double NAG_HUGE yb[],
		       double NAG_HUGE bc[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAZ_g)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAZ_fcnep)(Integer neq, double x, double eps, double NAG_HUGE y[],
			   double NAG_HUGE f[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAZ_fcnep)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAZ_fcna)(Integer neq, double x, double NAG_HUGE a[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAZ_fcna)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAZ_fcnb)(Integer neq, double x, double NAG_HUGE a[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAZ_fcnb)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAZ_jacobe)(Integer neq, double x, double eps, double NAG_HUGE y[],
			    double NAG_HUGE f[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAZ_jacobe)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAZ_jacobg)();
#else
typedef void (NAG_CALL * NAG_D02RAZ_jacobg)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAZ_jaceps)(Integer neq, double x, double eps, double NAG_HUGE y[],
			    double NAG_HUGE f[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAZ_jaceps)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAZ_jacgep) (Integer neq, double eps, double NAG_HUGE ya[], double NAG_HUGE yb[],
			    double NAG_HUGE bcep[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAZ_jacgep)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02PCC_FUN)(Integer neq, double t, double NAG_HUGE y[],
		       double NAG_HUGE yp[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02PCC_FUN)();
#endif

typedef NAG_D02PCC_FUN NAG_D02PDC_FUN;
typedef NAG_D02RAZ_fcn NAG_D02GAC_FUN;

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02GBC_FUN)(Integer neq, double x, double NAG_HUGE f[],
					 Nag_User NAG_HUGE *com);
#else
typedef void (NAG_CALL * NAG_D02GBC_FUN)();
#endif

#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02GBC_GFUN)(Integer neq, double x, double NAG_HUGE g[],
					 Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02GBC_GFUN)();
#endif




#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02QFC_FUN)(Integer neqf, double x, double NAG_HUGE y[],
					double NAG_HUGE f[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02QFC_FUN)();
#endif


#ifdef NAG_PROTO
typedef double (NAG_CALL * NAG_D02QFC_GFUN)(Integer neqf, double x, double NAG_HUGE y[],
					   double NAG_HUGE yp[], Integer k,
					   Nag_User NAG_HUGE *comm);
#else
typedef double (NAG_CALL * NAG_D02QFC_GFUN)();
#endif



#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAC_FUN)(Integer neq, double x, double eps,
					double NAG_HUGE y[], double NAG_HUGE f[],
					Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAC_FUN)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAC_GFUN)(Integer neq, double eps, double NAG_HUGE ya[],
					 double NAG_HUGE yb[], double NAG_HUGE bc[],
					 Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAC_GFUN)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAC_JFUN) (Integer neq, double x, double eps,
					    double NAG_HUGE y[], double NAG_HUGE f[],
					    Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAC_JFUN)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAC_JGFUN)(Integer neq, double eps,
					    double NAG_HUGE ya[], double NAG_HUGE yb[],
					    double NAG_HUGE aj[], double NAG_HUGE bj[],
					    Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAC_JGFUN)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAC_JEPSFUN)(Integer neq, double x, double eps,
					    double NAG_HUGE y[], double NAG_HUGE f[],
					    Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAC_JEPSFUN)();
#endif


#ifdef NAG_PROTO
typedef void (NAG_CALL * NAG_D02RAC_JGEPSFUN)(Integer neq, double eps,
					    double NAG_HUGE ya[], double NAG_HUGE yb[],
					    double NAG_HUGE bcep[], Nag_User NAG_HUGE *comm);
#else
typedef void (NAG_CALL * NAG_D02RAC_JGEPSFUN)();
#endif

#endif  /* not NAG_TYPES */
