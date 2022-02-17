#ifndef SRT_H_ALLTREFCT_H
#define SRT_H_ALLTREFCT_H

#define TREE_DOWN 0
#define TREE_MID 1
#define TREE_UP 2

double interp_linear(double x, double x_lower, double prem_lower,
                     double x_upper, double prem_upper);

/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :find_closest
  AUTHOR          :O. Van Eyseren
  DESCRIPTION     :find guess such that val_array[index] is the closest to
target CALL            :

<%%END---------------------------------------------------------------------*/

void find_closest(double key, double *val_array, long min_index, long max_index,
                  long *guess);

/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :find_closest_strictly_below
  AUTHOR          :O. Van Eyseren
  DESCRIPTION     :find guess such that
                        val_array[guess] < key < val_array[guess+1]
  MODIFIES	  :node->mid_son_index
  CALL            :

<%%END---------------------------------------------------------------------*/
void find_closest_strictly_below(double target, double *val_array,
                                 long min_index, long max_index, long *guess);

/** assume that we are probably increasing  , so we take as initial guess
    the value in node->mid_son_index.
    want index of element of prev_x closest to (drift_sam...STATEVAR);
    values in prev_x are increasing **/
/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :srt_f_trintrexindex
  AUTHOR          :E.Auld
  DESCRIPTION     :finds the index to center on in a trinomial tree in x
  MODIFIES	  :node->mid_son_index
  CALL            :

<%%END---------------------------------------------------------------------*/
void srt_f_trintrexindex(SrtStpPtr stp,  /* current step */
                         double *prev_x, /* state vars at previous step */
                         SrtTrinTreNdInf *node,   /* current node in tree */
                         SrtGrfnParam *grfnparam, /* grfn parameter */
                         int closestmeth /* CLOSESTINX or CLOSESTINLNX */
);

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        : srt_f_trintredcntvec
  AUTHOR          : O. Van Eyseren from E.Auld code
  DESCRIPTION     : discounts assets to cur (cur is prob. weighted sum of three
                        vectors in node.)
  MODIFIES        :cur
  CALL            :

<%%END>---------------------------------------------------------------------*/
void srt_f_trintredcntvec(SrtStpPtr stp, double *cur, SrtTrinTreNdInf *node,
                          double **assets, long dim);

/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :srt_f_trintredeltax
  AUTHOR          :O. Van Eyseren (from previous source code by E. Auld)
  DESCRIPTION     :calculates the delta_x step in a trinomial tree
                   this function assumes the the corresponding underlying
                   (specified by the und_index) is of the BS type)
  MODIFIES	  :
  CALL            :
  INPUT		  : und_index corresponds to the index of the underlying
                    used (in the order of the **TmInf)

<%%END---------------------------------------------------------------------*/
void srt_f_trintredeltax(SrtStpPtr top, double *delta_x, int und_index);

/** calculate probabilities **/
/*<%%STA-----------------------------------------------------------------
  FUNCNAME        :srt_f_trintreprob
  AUTHOR          :O. Van Eyseren (from E. Auld code)
  DESCRIPTION     :calculate probabilities  , assuming that node->drift_sam and
    node->var_at_sam are correctly set  , as well as node->mid_son_index
    If the mid_son_index hits the top (bottom)  , it is decremented
    (incremented).  Probabilities are bounded to be between PROBUPB and
    PROBLUB. They are printed if PRINFO is true.
  MODIFIES	  :node->p
  CALL            :

<%%END---------------------------------------------------------------------*/
void srt_f_trintreprob(double *prev_x, SrtTrinTreNdInf *node,
                       SrtGrfnParam *grfnparam);

/*-------------------------------------------------------------------------
  FUNCNAME        :srt_f_nonotreprob
  AUTHOR          :O. Van Eyseren (from E. Auld code)
  DESCRIPTION     :calculates probabilities for a nononomial tree  ,
                                   assuming that node -> moments are correctly
set  , as well as node -> mid_son_index MODIFIES		  :node ->
p[0..2][0..2]
---------------------------------------------------------------------------*/
void srt_f_nonotreprob(SrtTwoFacTreInf *trinf, SrtNonoTreeNodeInfo *node,
                       double **next_x);

/* Generates the matrix of state variables */
void srt_f_nonotregenx(SrtStpPtr stp, double **cur_x);

/* Calculate the discount asset price of the current node.  Asset prices
   are stored in the convention asset[x][y][dim] where dim is the width of
   GRFN tableau */
void srt_f_nonotredcntvec(SrtStpPtr stp, double *cur, SrtNonoTreeNodeInfo *node,
                          double ***assets, long dim);

#endif
