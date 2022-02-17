#ifndef _SRT_H_TWOFACTREFAC_H
#define _SRT_H_TWOFACTREFAC_H

/*
        FILE NAME: SRT_F_TWOFACTREFCT.C
        PURPOSE: UTILITY FUNCTIONS FOR LGM 2D NONO TREE

        NONO TREE FUNCTIONS
        ATTENTION: THESE FUNCTIONS HAVE BEEN MODIFIED AND WILL WORK ONLY WITH
   LGM 2D TREE BACKUP OF FORMER FUNCTIONS IS AVAILABLE IN
   SRT_F_TWOFACTREFCT_BACKUP.C (NOT IN PROJECT) AS WELL  , SRT_F_TREBDT2F.C HAS
   BEEN REMOVED FROM PROJECT

        THIS FILE CONTAINS TWO PARTS:
        -	FIRST PART CONTAINS FUNCTIONS USED TO COMPUTE PV OF ASSETS AT
   SOME NODE  , GIVEN PV OF ASSETS AT ALL NODES AT NEXT TIME STEP -	SECOND
   PART CONTAINS FUNCTIONS USED TO DESIGN THE GEOMETRY OF THE TREE
*/

/*
                                                                                PART I
                                        DISCOUNTING FROM TIME STEP N+1 TO NODE N
   ,I  ,J
*/

/* Computes forwards of statevars and stores in node->forward */
/* FORWARDS ARE EXPRESSED IN THE DUAL BASIS OF THE NEXT LAYER OF STEPS */
void lgm2f_forwards(
    SrtIRMTmInf *tminf, /* Time info  , ie. sigma  , tau and rho at time step */
    SrtNonoTreeNodeInfo *node, /* Output will be returned in node->forward */
    double dt,                 /* Delta t between this step and next */
    double nxt_inverse_basis[2][2]); /* Transfer Matrix to Dual basis at next
                                        layer of steps */

/* Establishes connectivity to next layer of steps  ,
   stores results in node->son_index and node->son_statevar */
/* Indexes and Statevars are expressed in dual basis
   (statevar[MID][0] will not be X1 at mid son node  , but eigenvector1 at mid
   son node */
/* Actually  , everything is expressed in dual basis */

/* ATTENTION:	THE SON INDEXES WILL BE CAPPED AND FLOORED BY TREE GEOMETRY AT
   NEXT TIME STEP; ON THE CONTRARY  , THE SON STATEVARS WILL BE KEPT AT THEIR
   TRUE VALUES  , EVEN IF THEY ARE OUTSIDE THE GRID. THIS MEANS THAT WE
   CORRECTLY DISCRETISE THE DISTRIBUTION (PROBAS IN lgm2f_probas)  , BUT WE
   SUPPOSE A FLAT EXTRAPOLATION OF THE ASSETS (IN lgm2f_discount_to_node) (ASK
   AS IF YOU REALLY WANT TO UNDERSTAND) */
void lgm2f_connect(
    SrtDiagTwoFacTreeInfo
        *nxttrinf,            /* Info about tree structure at next time step */
    SrtNonoTreeNodeInfo *node /* Input: forwards in node->forward */
    /* Output in node->son_index and node->son_statevar */
);

/*	Compute probabilities given moments and connectivity information */
void lgm2f_probas(
    SrtNonoTreeNodeInfo
        *node); /*	Requires connectivity (->son_index and ->son_statevar)
                                                and moments (->moments)
                                                computes probas (->p) */

/*	Computes PV at some node
        using information given by moments  , connectivity and probas at this
   node */
void lgm2f_discount_to_node(SrtNonoTreeNodeInfo *node, /* Node info */
                            double ***next_assets, /* Next assets [x1][x2][i] */
                            double *cur_assets,    /* Current assets at node */
                            long nb_assets);       /* Number of assets */

/*
                                                                                PART II
                                                                CONSTRUCTION OF
   THE TREE
*/

/*	Populates trinf and maxtrinf */
/*  At each time step:
        - The dual basis is computed: eigenvectors of local cov matrix * sqrt
   (eigenvalues)
        - The grid is constructed in dual basis with spacing NUMSTDEVINSPACING
   (in srt_h_trestruct.h)
        - The grid is trimmed at CUTTREEATSTDEV (srt_h_trestruct.h) GLOBAL
   standard deviations in each direction
*/
Err lgm2f_trelim(SrtStpPtr stp, SrtDiagTwoFacTreeInfo *maxtrinf);

#endif