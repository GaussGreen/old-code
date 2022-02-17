/*******************************************************************************
**
**              Copyright (c) 1993 PARIBAS Capital Markets Group
**
********************************************************************************
**
**      SYSTEM:         SRT     SORT  , Fixed Income 2020 Addins
**      SUB_SYSTEM:     PROB    Proba Tools
**
**      MODULE NAME:    PROBA_TOOLs
**
**      PURPOSE:        Import Header
**
**		INCLUDED BY:	proba_tools
**
**		AUTHORS:        Denis DESBIEZ
**
**      DATE:
**
**      DESCRIPTION:    Header for the library of general related option
**						functions used in the SORT program
*applications
*******************************************************************************/

/* ========================================================================== */

#ifndef PROBFNCTNS_H
#define PROBFNCTNS_H

/* ========================================================================== */

/*******************************************************************************
**                      Macros  , Typedefs and Constants
*******************************************************************************/

#define DEFAULT_PROB_TERMS 3

/******************************************************************************/

/* =============================================================================
  FUNCTION     : proba_max_fct(...)
 ============================================================================ */
double proba_max_fct(double spot, double barrier, double vol, double mat,
                     double disc);

/* =============================================================================
  FUNCTION     : proba_max_explo_fct(...)
 ============================================================================ */
double proba_max_explo_fct(double spot, double barrier, double vol, double mat,
                           double disc);

/* =============================================================================
  FUNCTION     : proba_min_fct(...)
 ============================================================================ */
double proba_min_fct(double spot, double barrier, double vol, double mat,
                     double disc);

/* =============================================================================
  FUNCTION     : proba_min_explo_fct(...)
 ============================================================================ */
double proba_min_explo_fct(double spot, double barrier, double vol, double mat,
                           double disc);

/* =============================================================================
  FUNCTION     : proba_gauss_fct(...)
 ============================================================================ */
double proba_gauss_fct(double spot, double b_up, double b_do, double vol,
                       double mat, double disc, int nb_term);

/* =============================================================================
  FUNCTION     : proba_gauss_explo_fct(...)
 ============================================================================ */
double proba_gauss_explo_fct(double spot, double b_up, double b_do, double vol,
                             double mat, double disc, int nb_term);

/* =============================================================================
  FUNCTION     : proba_yes_no_fct(...)
 ============================================================================ */
double probag_yes_no_fct(double spot, double b_yes, double b_no, double vol,
                         double mat, double disc, int nb_term);

/* =============================================================================
  FUNCTION     : proba_yes_no_explo_fct(...)
 ============================================================================ */
double probag_yes_no_explo_fct(double spot, double b_yes, double b_no,
                               double vol, double mat, double disc,
                               int nb_term);

/*	Numerical calculation of conditional correlation  ,
        i.e. correl ( X   , Y / abs ( X * Y ) > k
        where X and Y are Gaussians with correl rho */
Err proba_cond_corr(double mean1, double mean2, double std1, double std2,
                    double rho, double k, int npth, double *empave1,
                    double *empave2, double *empexpx1sqr, double *empexpx2sqr,
                    double *empexpx1x2);

#endif