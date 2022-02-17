/****************************************************************************/
/*      Set up core model-dependent functions                               */
/****************************************************************************/
/*      fix3model.c                                                         */
/****************************************************************************/



#include "fix123head.h"


/* Global model interface -- this should be wrapped in a FIX3_MODEL_INTERFACE
 * data structure, but this would require a lot of changes in products */

/* Volatility bootstrapping */
int     (*Fix3_BFactor)          (double*, long ,long, char, char, MKTVOL_DATA* , T_CURVE const*) = NULL;
int     (*Fix3_SpotVol)          (MKTVOL_DATA*, T_CURVE const*) = NULL;
int     (*Fix3_Interp_SpotVol)   (FIX3_TREE_DATA*, MKTVOL_DATA*) = NULL;
int     (*Fix3_IndexVol)         (double *, long, long, char, char,MKTVOL_DATA*,T_CURVE const*) = NULL;
int     (*Fix3_Filtered_SpotVol) (MKTVOL_DATA *, T_CURVE const*) = NULL;
int     (*Fix3_GenIndexVol)      (double *, long, long, long, long, char, char, MKTVOL_DATA *, T_CURVE const* ) = NULL;
int     (*Fix3_IndexLimits)      (double* ,double*, long, long, long, int, char, char,  MKTVOL_DATA*, T_CURVE const*) = NULL;
int     (*Fix3_Cet_Main)         (int,T_CURVE const*,MKTVOL_DATA *,FIX3_TREE_DATA const*) = NULL;
/* Tree construction */
int     (*Fix3_Tree_Limits)      (int *,int *,int *,int *,int **,int **,int ***,int ***,int, FIX3_TREE_DATA *tree_data, MKTVOL_DATA *mktvol_data) = NULL;
int     (*Fix3_Drift_1D)         (MKTVOL_DATA *,FIX3_TREE_DATA *) = NULL;
int     (*Fix3_Drift_2D)         (MKTVOL_DATA *,FIX3_TREE_DATA *) = NULL;
int     (*Fix3_Drift_3D)         (MKTVOL_DATA *,FIX3_TREE_DATA *) = NULL;
int     (*Fix3_Lattice)          (FIX3_DEV_DATA *,int,int,MKTVOL_DATA *,FIX3_TREE_DATA *) = NULL;
int     (*Fix3_Dev)              (double *,int,int,int,FIX3_DEV_DATA const*,FIX3_TREE_DATA const*) = NULL;
/* Model I/O */
int     (*Fix3_Param_Check)      (int,MKTVOL_DATA *,FIX3_TREE_DATA *) = NULL;
int     (*Fix3_Cet_Model_Manager)(MKTVOL_DATA const*,FIX3_TREE_DATA const*,MKTVOL_DATA *,FIX3_TREE_DATA *) = NULL;  
int     (*Fix3_PackModelData)    (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *) = NULL;
int     (*Fix3_UnPackModelData)  (MKTVOL_DATA *,FIX3_TREE_DATA *,long *,double *,long *,double *) = NULL;


int     Fix3_Model_Interface_Init (int ModelChoice)
{
    switch (ModelChoice)
    {
    case FIX3_ORIGINAL:

        Fix3_BFactor           = Fix3_BFactor_Classic;
        Fix3_SpotVol           = Fix3_SpotVol_Classic;
        Fix3_Interp_SpotVol    = Fix3_Interp_SpotVol_Classic;
        Fix3_IndexVol          = Fix3_IndexVol_Classic;
        Fix3_Filtered_SpotVol  = Fix3_Filtered_SpotVol_Classic;
        Fix3_GenIndexVol       = Fix3_GenIndexVol_Classic;
        Fix3_IndexLimits       = Fix3_IndexLimits_Classic;
        Fix3_Cet_Main          = Fix3_Cet_Classic;

        Fix3_Tree_Limits       = Fix3_Tree_Limits_Classic;
        Fix3_Drift_1D          = Fix3_Drift_1D_Classic;
        Fix3_Drift_2D          = Fix3_Drift_2D_Classic;
        Fix3_Drift_3D          = Fix3_Drift_3D_Classic;
        Fix3_Lattice           = Fix3_Lattice_Classic;
        Fix3_Dev               = Fix3_Dev_Classic;

        Fix3_Param_Check       = Fix3_Param_Check_Original;
        Fix3_Cet_Model_Manager = Fix3_Cet_Model_Manager_Classic;
        Fix3_PackModelData     = Fix3_PackModelData_Classic;
        Fix3_UnPackModelData   = Fix3_UnPackModelData_Classic;
        
        return (SUCCESS);

    case FIX3_CLASSIC:

        Fix3_BFactor           = Fix3_BFactor_Classic;
        Fix3_SpotVol           = Fix3_SpotVol_Classic;
        Fix3_Interp_SpotVol    = Fix3_Interp_SpotVol_Classic;
        Fix3_IndexVol          = Fix3_IndexVol_Classic;
        Fix3_Filtered_SpotVol  = Fix3_Filtered_SpotVol_Classic;
        Fix3_GenIndexVol       = Fix3_GenIndexVol_Classic;
        Fix3_IndexLimits       = Fix3_IndexLimits_Classic;
        Fix3_Cet_Main          = Fix3_Cet_Classic;

        Fix3_Tree_Limits       = Fix3_Tree_Limits_Classic;
        Fix3_Drift_1D          = Fix3_Drift_1D_Classic;
        Fix3_Drift_2D          = Fix3_Drift_2D_Classic;
        Fix3_Drift_3D          = Fix3_Drift_3D_Classic;
        Fix3_Lattice           = Fix3_Lattice_Classic;
        Fix3_Dev               = Fix3_Dev_Classic;

        Fix3_Param_Check       = Fix3_Param_Check_Classic;
        Fix3_Cet_Model_Manager = Fix3_Cet_Model_Manager_Classic;
        Fix3_PackModelData     = Fix3_PackModelData_Classic;
        Fix3_UnPackModelData   = Fix3_UnPackModelData_Classic;
        
        return (SUCCESS);

    case FIX3_TIMEDEP:

        Fix3_BFactor           = Fix3_BFactor_TimeDep;
        Fix3_SpotVol           = Fix3_SpotVol_TimeDep;
        Fix3_Interp_SpotVol    = Fix3_Interp_SpotVol_TimeDep; 
        Fix3_IndexVol          = Fix3_IndexVol_TimeDep;
        Fix3_Filtered_SpotVol  = Fix3_Filtered_SpotVol_TimeDep;
        Fix3_GenIndexVol       = NULL;
        Fix3_IndexLimits       = Fix3_IndexLimits_Classic;
        Fix3_Cet_Main          = Fix3_Cet_Classic;

        Fix3_Tree_Limits       = Fix3_Tree_Limits_TimeDep;
        Fix3_Drift_1D          = Fix3_Drift_1D_TimeDep;
        Fix3_Drift_2D          = Fix3_Drift_2D_TimeDep;
        Fix3_Drift_3D          = Fix3_Drift_3D_TimeDep;
        Fix3_Lattice           = Fix3_Lattice_TimeDep;
        Fix3_Dev               = Fix3_Dev_Classic;

        Fix3_Param_Check       = Fix3_Param_Check_TimeDep;
        Fix3_Cet_Model_Manager = Fix3_Cet_Model_Manager_TimeDep;
        Fix3_PackModelData     = Fix3_PackModelData_TimeDep;
        Fix3_UnPackModelData   = Fix3_UnPackModelData_TimeDep;


        return (SUCCESS);

    case FIX3_SMD:

        Fix3_BFactor           = Fix3_BFactor_Classic;
        Fix3_SpotVol           = Fix3_SpotVol_Smd;
        Fix3_Interp_SpotVol    = Fix3_Interp_SpotVol_Smd;
        Fix3_IndexVol          = Fix3_IndexVol_Smd;
        Fix3_Filtered_SpotVol  = Fix3_SpotVol_Smd; /* Smd does not filter */
        Fix3_GenIndexVol       = NULL;
        Fix3_IndexLimits       = Fix3_IndexLimits_Classic;
        Fix3_Cet_Main          = Fix3_Cet_Classic;

        Fix3_Tree_Limits       = Fix3_Tree_Limits_Smd;
        Fix3_Drift_1D          = NULL;
        Fix3_Drift_2D          = Fix3_Drift_2D_Smd;
        Fix3_Drift_3D          = NULL;
        Fix3_Lattice           = Fix3_Lattice_Smd;
        Fix3_Dev               = Fix3_Dev_Classic;

        Fix3_Param_Check       = Fix3_Param_Check_Smd;
        Fix3_Cet_Model_Manager = Fix3_Cet_Model_Manager_Smd;
        Fix3_PackModelData     = Fix3_PackModelData_Smd;
        Fix3_UnPackModelData   = Fix3_UnPackModelData_Smd;

        return (SUCCESS);

    case FIX3_TMX:

        Fix3_BFactor           = Fix3_BFactor_Tmx;
        Fix3_SpotVol           = Fix3_SpotVol_Tmx;
        Fix3_Interp_SpotVol    = Fix3_Interp_SpotVol_Tmx;
        Fix3_IndexVol          = Fix3_IndexVol_Tmx;
        Fix3_Filtered_SpotVol  = Fix3_SpotVol_Tmx; /* Tmx does not filter */
        Fix3_GenIndexVol       = NULL;
        Fix3_IndexLimits       = Fix3_IndexLimits_Classic;
        Fix3_Cet_Main          = Fix3_Cet_Tmx;

        Fix3_Tree_Limits       = Fix3_Tree_Limits_Classic;
        Fix3_Drift_1D          = Fix3_Drift_1D_Tmx;
        Fix3_Drift_2D          = NULL;
        Fix3_Drift_3D          = NULL;
        Fix3_Lattice           = Fix3_Lattice_Tmx;
        Fix3_Dev               = Fix3_Dev_Tmx;

        Fix3_Param_Check       = Fix3_Param_Check_Tmx;
        Fix3_Cet_Model_Manager = NULL;
        Fix3_PackModelData     = Fix3_PackModelData_Tmx;
        Fix3_UnPackModelData   = Fix3_UnPackModelData_Tmx;

        return (SUCCESS);

     case FIX3_E2Q:

        Fix3_BFactor           = Fix3_BFactor_E2Q;
        Fix3_SpotVol           = Fix3_SpotVol_E2Q;
        Fix3_Interp_SpotVol    = Fix3_Interp_SpotVol_E2Q;
        Fix3_IndexVol          = Fix3_IndexVol_E2Q;
        Fix3_Filtered_SpotVol  = Fix3_Filtered_SpotVol_E2Q;
        Fix3_GenIndexVol       = Fix3_GenIndexVol_E2Q;
        Fix3_IndexLimits       = Fix3_IndexLimits_Classic;
        Fix3_Cet_Main          = Fix3_Cet_Classic;

        Fix3_Tree_Limits       = Fix3_Tree_Limits_Classic;
        Fix3_Drift_1D          = Fix3_Drift_1D_E2Q;
        Fix3_Drift_2D          = Fix3_Drift_2D_E2Q;
        Fix3_Drift_3D          = Fix3_Drift_3D_E2Q;
        Fix3_Lattice           = Fix3_Lattice_E2Q;
        Fix3_Dev               = Fix3_Dev_Classic;

        Fix3_Param_Check       = Fix3_Param_Check_Classic;
        Fix3_Cet_Model_Manager = Fix3_Cet_Model_Manager_E2Q;
        Fix3_PackModelData     = Fix3_PackModelData_E2Q;
        Fix3_UnPackModelData   = Fix3_UnPackModelData_E2Q;

        return (SUCCESS);

    default:

        return (FAILURE);
    }

}
