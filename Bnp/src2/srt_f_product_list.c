/* =============================================================================
   FILENAME:   srt_f_product_list.c

   PURPOSE:    pre-initialized products routines
   =============================================================================
 */

#include "srt_h_all.h"

/* ------ The Static Pointer used to refer to the Products List ----------------
 */
static SrtList *_srt_products_list = NULL;

Err create_product_list(char *prod_list_name) {
  Err err = srt_f_lstcreate(&_srt_products_list, prod_list_name);
  if (err)
    return serror("Could not create product list %s", prod_list_name);

  return NULL;
}

Err destroy_all_products() {
  Err err = srt_f_lstfree(_srt_products_list, SRT_YES);
  if (err)
    return serror("%s in destroy_all_products", err);

  srt_free(_srt_products_list);

  return NULL;
}

SrtList *get_product_list() { return _srt_products_list; }

static char *destroy_product(void *obj) {
  SrtProduct *product = (SrtProduct *)obj;
  return product->Destroy(product);
}

Err add_product_to_list(SrtProduct *product) {
  return srt_f_lstins(get_product_list(), product->product_name, 0.0,
                      OBJ_PTR_UND, product, destroy_product,
                      &(product->product_ticker));
}

SrtProduct *lookup_product(char *name) {
  Err err = NULL;
  char copy_name[SRTBUFSZ];
  SrtListObject *obj;

  strcpy(copy_name, name);
  strupper(copy_name);
  strip_white_space(copy_name);
  rem_tick_string(copy_name, copy_name);

  err = srt_f_lstgetobj(*(get_product_list()), copy_name, 0.0, &obj);
  if (err)
    return NULL;

  return (SrtProduct *)(obj->val.pval);
}
