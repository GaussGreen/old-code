/* =============================================================================
   FILENAME:   srt_f_product_list.h

   PURPOSE:    pre-initialized products routines
   =============================================================================
 */

#ifndef __SRT_H_PRODUCT_LIST_H__
#define __SRT_H_PRODUCT_LIST_H__

typedef struct _SrtProduct SrtProduct;

struct _SrtProduct {
  char product_name[SRTBUFSZ];
  char product_lbl[SRTBUFSZ];
  long product_ticker;

  Err (*Destroy)(SrtProduct *); /* Destructor */

  /* The following function is to be called by a pricer (GRFN) to obtain dates
  for which DFs need to be calculated and passed to the Payoff function.
  This pointer can be NULL if the product doesn't need any additional DFs to
  calculate its payoff.

  Parameter 1: "this"
  Parameter 2: fixing date
  Parameter 3: number of underlyings involved in payoff calculation is stored
  here Parameter 4: number of DFs per underlying is stored here Parameter 5: DFs
  maturities are stored here

  Parameters 4 and 5 should be allocated inside and need to be freed by caller.
  Number of underlyings returned can be 1 or 2. First one is domestic and
  the second one is foreign (if the product is quanto) */

  Err (*RequestDfDates)(SrtProduct *, long, int *, int **, long ***);

  /* The following function calculates the payoff of the product at a specified
  date.

  Parameter 1: "this"
  Parameter 2: date at which the payoff is calculated
  Parameter 3: fixing date
  Parameter 4: Requested DFs for each underlying (can be NULL if not needed)
  Parameter 5: Calculated payoff is stored here */

  Err (*Payoff)(SrtProduct *, long, long, double **, double *);

  void *spec_desc; /* product dependent description */
};

Err create_product_list(char *prod_list_name);
Err destroy_all_products();
SrtList *get_product_list();

Err add_product_to_list(SrtProduct *product);
SrtProduct *lookup_product(char *name);

#endif /* #define __SRT_H_PRODUCT_LIST_H__ */