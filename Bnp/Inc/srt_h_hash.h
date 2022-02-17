// srt_h_hash.h:  Hash functionality

#ifndef __SRT_H_HASH__
#define __SRT_H_HASH__

#define Hash_MAXNPOW 50
#define Poly_MAXNFAC 20

typedef struct _SHashTree SHashTree;
struct _SHashTree {
	SHashTree		*facpow[Hash_MAXNPOW];
	int				idx;
};

typedef struct _SHash {
	SHashTree		*hash_tree;
	int				nfac, pos, maxsize;
} SHash;

typedef struct _SPolynomial {
	int				n_mono, maxmono;
	int				**mono;
	double			**coef;
	int				*pool;
	int				nstp;
	SHash			*hash;
} SPolynomial;

Err InitHash(SHash *hash, int nfac, int maxsize);
Err Hash(SHash *hash, int *imono, int *hmono);
void FreeHashTree(SHashTree *p);

Err Poly_Init(SPolynomial *poly, SHash *hash, int maxmono, int nstp);
void Poly_Free(SPolynomial *poly);
void Poly_Copy(SPolynomial *poly, SPolynomial *p);
void Poly_Clear(SPolynomial *poly);
Err Poly_AddDPoly(SPolynomial *poly, SPolynomial *p, int ifac);
Err Poly_AddConstMono(SPolynomial *poly, double c, double *pc, int *mono_idx);
Err Poly_AddConstMonoPoly(SPolynomial *poly, double c, double *pc, int *mono_idx, SPolynomial *p);
Err Poly_AddConstMonoPolyPoly(SPolynomial *poly, double c, double *pc, int *mono_idx,
									 SPolynomial *p1, SPolynomial *p2);


#endif  // #ifndef __SRT_H_HASH__