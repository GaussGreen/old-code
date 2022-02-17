// srt_f_hash.cxx:  Hash function

#include "srt_h_all.h"
#include "srt_h_hash.h"

Err InitHash(SHash *hash, int nfac, int maxsize) {
  hash->hash_tree = (SHashTree *)calloc(1, sizeof(SHashTree));
  if (!hash->hash_tree)
    return serror("Memory failure");

  hash->maxsize = maxsize;
  hash->nfac = nfac;
  hash->pos = 1;

  return NULL;
}

Err Hash(SHash *hash, int *imono, int *hmono) {
  SHashTree *p = hash->hash_tree;
  int i, j, k;

  for (j = hash->nfac - 1; j >= 0 && imono[j] == 0; j--)
    ;

  for (i = 0; i <= j; i++) {
    k = (imono[i] >= 0 ? imono[i] : Hash_MAXNPOW + imono[i]);
    if (!p->facpow[k]) {
      p->facpow[k] = (SHashTree *)calloc(1, sizeof(SHashTree));
      if (!p->facpow[k])
        return serror("Memory failure");
      p->facpow[k]->idx = -1;
    }
    p = p->facpow[k];
  }
  if (p->idx < 0) {
    if (hash->pos >= hash->maxsize)
      return serror("Hash maxsize exceeded");
    p->idx = hash->pos++;
  }
  *hmono = p->idx;

  return NULL;
}

void FreeHashTree(SHashTree *p) {
  int i;
  if (p)
    for (i = 0; i < Hash_MAXNPOW; i++)
      FreeHashTree(p->facpow[i]);
  free(p);
}

Err Poly_Init(SPolynomial *poly, SHash *hash, int maxmono, int nstp) {
  poly->hash = hash;
  poly->maxmono = maxmono;
  poly->nstp = nstp;
  poly->mono = imatrix(0, maxmono - 1, 0, hash->nfac - 1);
  poly->coef = dmatrix(0, maxmono - 1, 0, nstp - 1);
  poly->pool = (int *)calloc(hash->maxsize, sizeof(int));

  if (!poly->mono || !poly->coef || !poly->pool)
    return serror("Memory failure");
  memset(poly->pool, -1, hash->maxsize * sizeof(int));

  return NULL;
}

void Poly_Free(SPolynomial *poly) {
  if (poly->mono)
    free_imatrix(poly->mono, 0, poly->maxmono - 1, 0, poly->hash->nfac - 1);
  if (poly->coef)
    free_dmatrix(poly->coef, 0, poly->maxmono - 1, 0, poly->nstp - 1);
  free(poly->pool);
}

void Poly_Copy(SPolynomial *poly, SPolynomial *p) {
  poly->n_mono = p->n_mono;
  memcpy(&poly->mono[0][0], &p->mono[0][0],
         p->n_mono * p->hash->nfac * sizeof(int));
  memcpy(&poly->coef[0][0], &p->coef[0][0],
         p->n_mono * p->nstp * sizeof(double));
  memcpy(poly->pool, p->pool, p->hash->pos * sizeof(int));
}

void Poly_Clear(SPolynomial *poly) {
  memset(&poly->coef[0][0], 0, poly->n_mono * poly->nstp * sizeof(double));
  memset(poly->pool, -1, poly->hash->pos * sizeof(int));
  poly->n_mono = 0;
}

Err Poly_AddDPoly(SPolynomial *poly, SPolynomial *p, int ifac) {
  Err err;
  int j, k;
  int imono[Poly_MAXNFAC], hmono;

  for (j = 0; j < p->n_mono; j++) {
    if (p->mono[j][ifac] != 0) {
      memcpy(imono, p->mono[j], poly->hash->nfac * sizeof(int));
      imono[ifac]--;
      err = Hash(poly->hash, imono, &hmono);
      if (err)
        return err;

      if (poly->pool[hmono] < 0) {
        if (poly->n_mono >= poly->maxmono)
          return serror("Max mono exceeded");
        poly->pool[hmono] = poly->n_mono++;
        memcpy(poly->mono[poly->n_mono - 1], imono,
               poly->hash->nfac * sizeof(int));
      }
      for (k = 0; k < poly->nstp; k++)
        poly->coef[poly->pool[hmono]][k] += p->mono[j][ifac] * p->coef[j][k];
    }
  }
  return NULL;
}

Err Poly_AddConstMono(SPolynomial *poly, double c, double *pc, int *mono_idx) {
  Err err;
  int k;
  int hmono = 0;

  if (mono_idx) {
    err = Hash(poly->hash, mono_idx, &hmono);
    if (err)
      return err;
  }

  if (poly->pool[hmono] < 0) {
    if (poly->n_mono >= poly->maxmono)
      return serror("Max mono exceeded");
    poly->pool[hmono] = poly->n_mono++;
    memcpy(poly->mono[poly->n_mono - 1], mono_idx,
           poly->hash->nfac * sizeof(int));
  }
  for (k = 0; k < poly->nstp; k++)
    poly->coef[poly->pool[hmono]][k] += (pc ? c * pc[k] : c);

  return NULL;
}

Err Poly_AddConstMonoPoly(SPolynomial *poly, double c, double *pc,
                          int *mono_idx, SPolynomial *p) {
  Err err;
  int j, k;
  int imono[Poly_MAXNFAC], hmono;

  for (j = 0; j < p->n_mono; j++) {
    for (k = 0; k < poly->hash->nfac; k++)
      imono[k] = (mono_idx ? mono_idx[k] : 0) + p->mono[j][k];
    err = Hash(poly->hash, imono, &hmono);
    if (err)
      return err;

    if (poly->pool[hmono] < 0) {
      if (poly->n_mono >= poly->maxmono)
        return serror("Max mono exceeded");
      poly->pool[hmono] = poly->n_mono++;
      memcpy(poly->mono[poly->n_mono - 1], imono,
             poly->hash->nfac * sizeof(int));
    }
    for (k = 0; k < poly->nstp; k++)
      poly->coef[poly->pool[hmono]][k] += (pc ? c * pc[k] : c) * p->coef[j][k];
  }
  return NULL;
}

Err Poly_AddConstMonoPolyPoly(SPolynomial *poly, double c, double *pc,
                              int *mono_idx, SPolynomial *p1, SPolynomial *p2) {
  Err err;
  int i, j, k;
  int imono[Poly_MAXNFAC], hmono;

  for (i = 0; i < p1->n_mono; i++) {
    for (j = 0; j < p2->n_mono; j++) {
      for (k = 0; k < poly->hash->nfac; k++)
        imono[k] =
            (mono_idx ? mono_idx[k] : 0) + p1->mono[i][k] + p2->mono[j][k];
      err = Hash(poly->hash, imono, &hmono);
      if (err)
        return err;

      if (poly->pool[hmono] < 0) {
        if (poly->n_mono >= poly->maxmono)
          return serror("Max mono exceeded");
        poly->pool[hmono] = poly->n_mono++;
        memcpy(poly->mono[poly->n_mono - 1], imono,
               poly->hash->nfac * sizeof(int));
      }
      for (k = 0; k < poly->nstp; k++)
        poly->coef[poly->pool[hmono]][k] +=
            (pc ? c * pc[k] : c) * p1->coef[i][k] * p2->coef[j][k];
    }
  }
  return NULL;
}
