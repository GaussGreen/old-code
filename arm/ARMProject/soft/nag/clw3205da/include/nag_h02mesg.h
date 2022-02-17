/* <nag_h02mesg.h>
 *
 * Copyright 1997 Numerical Algorithms Group.
 *
 * Mark 5, 1997.
 *
 */

#ifndef NAG_H02MESG
#define NAG_H02MESG

#define NM_NO_MESG 0
#define NM_FEAS_TOL 1
#define NM_FIRST_SOLN 2
#define NM_MAX_DEPTH 3
#define NM_MAX_NODES 4
#define NM_MAX_NODES_ALL 5
#define NM_INT_TOL 6
#define NM_INT_OBJ_BOUND 7
#define NM_SOLN_TOL 8
#define NM_NODSEL 9
#define NM_VARSEL 10
#define NM_BRANCH_DIR 11
#define NM_MIP_PRINT_LEVEL 12
#define NM_H02_NUMVAR 13
#define NM_NUM_INTVAR 14
#define NM_LOWER_MA 15
#define NM_UPPER_MA 16
#define NM_H02_PROBNAME 17
#define NM_BNDNAME 18
#define NM_OBJNAME 19
#define NM_RANGENAME 20
#define NM_RHSNAME 21
#define NM_QSECTNAME 22
#define NM_H02_COL_LO_DEFAULT 23
#define NM_H02_COL_UP_DEFAULT 24
#define NM_H02_INFINITY 25
#define NM_H02_NCOL_APPROX 26
#define NM_H02_NROW_APPROX 27
#define NM_H02_MPS_SPARSE 28
#define NM_H02_MPS_OUTLEVEL 29
#define NM_MPS_ECHO_LINE 30
#define NM_MPS_FINISHED_READ 31
#define NM_MPS_FINISHED_READ_INT 32
#define NM_MPS_NAMES_SELECT 33
#define NM_MPS_DATA_ASSIGNED 34
#define NM_MIP_VAR_HEAD 35
#define NM_MIP_LCON_HEAD 36
#define NM_MIP_NOBND 37
#define NM_MIP_LOBND 38
#define NM_MIP_UPBND 39
#define NM_MIP_BOTHBND 40
#define NM_MIP_ROOT_SOLN_RES 41
#define NM_MIP_OPT_ROOT_SOLN 42
#define NM_MIP_ROOT_OBJ 43
#define NM_MIP_ROOT_TERMINATED 44
#define NM_MIP_BBEXIT_NODES 45
#define NM_MIP_BBEXIT_ONE_NODE 46
#define NM_MIP_OPT_IP 47
#define NM_MIP_FIRST_IP 48
#define NM_MIP_SUBOPT_IP 49
#define NM_MIP_NODE_HEAD 50
#define NM_MIP_NODE_SUMMARY 51
#define NM_MIP_ROOTNODE_SUMMARY 52
#define NM_MIP_NO_IP_SOLN 53
#define NM_MIP_INT_SOL_FOUND 54
#define NM_H02_LAMBDA 55
#define NM_H02_CRNAMES 56
#define NM_H02_NEWLINE 57
#define NM_H02_INTRO 58
#define NM_H02_SET_TEXT 59
#define NM_H02_SET_DOUBLE 60
#define NM_H02_SET_INT 61
#define NM_H02_RESET_DEF 62
#define NM_H02_RESET 63
#define NM_H02_TEXT_VALID_RANGE 64
#define NM_H02_SET_DOUBLE_MIN 65
#define NM_H02_BND_NAME 66
#define NM_H02_END_LINE 67
#define NM_H02_FINAL_OBJ 68
#define NM_H02_FINAL_STATUS 69
#define NM_H02_FINAL_SUMINF 70
#define NM_H02_HROWS 71
#define NM_H02_INF_BOUND 72
#define NM_H02_MACH_PREC 73
#define NM_H02_MAX_DF 74
#define NM_H02_MAX_ITER_1 75
#define NM_H02_MEM_ALLOC 76
#define NM_H02_MIN_SUMINF 77
#define NM_H02_NUM_LIN_CON 78
#define NM_H02_OBJ_NAME 79
#define NM_H02_OUTFILE 80
#define NM_H02_PARAM_TITLE 81
#define NM_H02_PRINT_LEVEL 82
#define NM_H02_PROB 83
#define NM_H02_PROB_NAME 84
#define NM_H02_RANGE_NAME 85
#define NM_H02_RANK_TOL 86
#define NM_H02_RHS_NAME 87
#define NM_H02_SOLN_RES 88
#define NM_H02_STATE 89

#ifdef NAG_MESG
char *nag_h02mesg[] =
{
  ": Dummy message for Chapter h02",
  "NM_FEAS_TOL: feas_tol................%9.2e",
  "NM_FIRST_SOLN: first_soln..............    %s",
  "NM_MAX_DEPTH:     max_depth...............%9ld\n",
  "NM_MAX_NODES: max_nodes...............%9ld",
  "NM_MAX_NODES_ALL: max_nodes...............ALL_NODES",
  "NM_INT_TOL:     int_tol.................%9.2e\n",
  "NM_INT_OBJ_BOUND: int_obj_bound...........%9.2e",
  "NM_SOLN_TOL:     soln_tol................%9.2e\n",
  "NM_NODSEL: nodsel......%s",
  "NM_VARSEL:     varsel...........%s\n",
  "NM_BRANCH_DIR: branch_dir.......%s",
  "NM_MIP_PRINT_LEVEL: print_level...%s",
  "NM_H02_NUMVAR: Number of variables...........%3ld\n",
  "NM_NUM_INTVAR: Number of integer variables...%3ld\n\n",
  "NM_LOWER_MA: lower...................     %s\n",
  "NM_UPPER_MA: upper...................     %s\n",
  "NM_H02_PROBNAME: PROBLEM name............ %s",
  "NM_BNDNAME: BOUNDS set.............. %s",
  "NM_OBJNAME:     objective............... %s\n",
  "NM_RANGENAME: RANGE set............... %s",
  "NM_RHSNAME:     RHS set................. %s\n",
  "NM_QSECTNAME: QSECTION name........... %s",
  "NM_H02_COL_LO_DEFAULT: col_lo_default..........%9.2e",
  "NM_H02_COL_UP_DEFAULT:     col_up_default..........%9.2e\n",
  "NM_H02_INFINITY: infinity................%9.2e\n",
  "NM_H02_NCOL_APPROX: ncol_approx.............%9ld",
  "NM_H02_NROW_APPROX:     nrow_approx.............%9ld\n",
  "NM_H02_MPS_SPARSE: sparse..................    %s    ",
  "NM_H02_MPS_OUTLEVEL: output_level......%s",
  "NM_MPS_ECHO_LINE: %s",
  "NM_MPS_FINISHED_READ: \n\nMPS file successfully read.\n\
Number of lines read: %9ld\nNumber of columns:    %9ld\n\
Number of rows:       %9ld  (including objective row)\n",  
  "NM_MPS_FINISHED_READ_INT: \n\nThe MPS file has successfully been read.\n\
Number of lines read: %9ld\nNumber of columns:    %9ld  (of which %1ld are integers)\n\
Number of rows:       %9ld  (including objective row)\n",  
  "NM_MPS_NAMES_SELECT: \nMPS Names Selected:\n\
Problem   %8s\nObjective %8s  RHS       %8s\nRANGES    %8s  BOUNDS    %8s\n",
  "NM_MPS_DATA_ASSIGNED: \nMPS data successfully assigned to problem data.\n",
  "NM_MIP_VAR_HEAD: \nVarbl   State    Value       Lower Bound  Upper Bound    \
Lagr Mult    Residual\n\n",
  "NM_MIP_LCON_HEAD: \nConstr  State    Value       Lower Bound  Upper Bound    \
Lagr Mult    Residual\n\n",
  "NM_MIP_NOBND: %-8.8s  %.2s %13.5e      None         None      %11.3e %11.3e\n",
  "NM_MIP_LOBND: %-8.8s  %.2s %13.5e  %12.4e      None     %11.3e %11.3e\n",  
  "NM_MIP_UPBND: %-8.8s  %.2s %13.5e      None     %12.4e  %11.3e %11.3e\n",           
  "NM_MIP_BOTHBND: %-8.8s  %.2s %13.5e  %12.4e %12.4e  %11.3e %11.3e\n",               
  "NM_MIP_ROOT_SOLN_RES: \nSolution of root problem:\n",
  "NM_MIP_OPT_ROOT_SOLN: \nOptimal solution of root %s problem found.\n",
  "NM_MIP_ROOT_OBJ: \nObjective value of root %s problem = %15.7e\n\n",
  "NM_MIP_ROOT_TERMINATED: \nRoot problem terminated.\n",
  "NM_MIP_BBEXIT_NODES: \nExit from branch and bound tree search after \
%1ld nodes.\n",
  "NM_MIP_BBEXIT_ONE_NODE: \nExit from branch and bound tree search after \
1 node.\n",
  "NM_MIP_OPT_IP: \nOptimal IP solution found.\n",
  "NM_MIP_FIRST_IP: \nFirst IP solution found - search terminated as requested.\n",
  "NM_MIP_SUBOPT_IP: \nIP solution found - may be suboptimal.\n",
  "NM_MIP_NODE_HEAD: \n  Node Parent     Obj     Varbl     Value     Lower     \
Upper      Value Depth  \n\
   No   Node     Value   Chosen    Before     Bound     Bound      After\n",
  "NM_MIP_NODE_SUMMARY: %5ld %5ld %s%5ld %9.2e %s%s %9.2e %5ld\n",
  "NM_MIP_ROOTNODE_SUMMARY: %5ld       %s\n",
  "NM_MIP_NO_IP_SOLN: \nNo IP solution was found.\n",
  "NM_MIP_INT_SOL_FOUND:    *** Integer Solution ***\n\n",
  "NM_H02_LAMBDA: lambda..................     %s\n",
  "NM_H02_CRNAMES:     crnames..............%s\n",
  "NM_H02_NEWLINE: \n",

  /* h02xyc introduction */
  "NM_H02_INTRO: \nOptional parameter setting for %s.\n\
--------------------------------------\n\n\
Option file: %s\n\n",

  /* Option range checking output */
  "NM_H02_SET_TEXT: %s set to %s\n",
  "NM_H02_SET_DOUBLE: %s set to %7.2e\n",
  "NM_H02_SET_INT: %s set to %1ld\n",
  "NM_H02_RESET_DEF: This change will cause the default value of %s\n\
to be reset to %7.2e.\n",
  "NM_H02_RESET: The value of %s has been reset from %7.2e to %7.2e\n\
because the valid range for %s is: %s.\n",
  "NM_H02_TEXT_VALID_RANGE: The valid range for %s is: %s.\n",
  "NM_H02_SET_DOUBLE_MIN: %s set to its minimum value %7.2e\n",
  "NM_H02_BND_NAME: bnd_name................ %8.8s\n",  
  "NM_H02_END_LINE: \n",
  "NM_H02_FINAL_OBJ: \nFinal %s objective value = %15.7e\n\n",
  "NM_H02_FINAL_STATUS: --++FRLLULEQTF",
  "NM_H02_FINAL_SUMINF: \nFinal sum of infeasibilities = %15.7e\n",
  "NM_H02_HROWS: hrows...................%9ld",
  "NM_H02_INF_BOUND: inf_bound...............%9.2e    ",
  "NM_H02_MACH_PREC:     machine precision.......%9.2e\n",
  "NM_H02_MAX_DF: max_df..................%9ld\n",
  "NM_H02_MAX_ITER_1: max_iter................%9ld",
  "NM_H02_MEM_ALLOC: \nMemory allocation:\n",
  "NM_H02_MIN_SUMINF: \nMinimum sum of infeasibilities = %15.7e\n",
  "NM_H02_NUM_LIN_CON: Linear constraints............%3ld    ",
  "NM_H02_OBJ_NAME: obj_name................ %8.8s    ",
  "NM_H02_OUTFILE: outfile.................%9s\n",
  "NM_H02_PARAM_TITLE: \nParameters to %s\n--------------------\n\n",
  "NM_H02_PRINT_LEVEL: print_level...%s",
  "NM_H02_PROB: prob....................%s    ",
  "NM_H02_PROB_NAME: prob_name............... %8.8s\n",
  "NM_H02_RANGE_NAME: range_name.............. %8.8s    ",
  "NM_H02_RANK_TOL: rank_tol................%9.2e    ",
  "NM_H02_RHS_NAME: rhs_name................ %8.8s\n",
  "NM_H02_SOLN_RES: \nFinal solution:\n",
  "NM_H02_STATE: state...................     %s\n",

/* END OF MESSAGE STRINGS */
  ""
};
#else
extern char *nag_h02mesg[];
#endif

#endif /* not NAG_H02MESG */



