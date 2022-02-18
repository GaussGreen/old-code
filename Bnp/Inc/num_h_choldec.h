#if !defined(_NUMHCHOLDEC_H_)
#define _NUMHCHOLDEC_H_

void CholDec(
    void (*CholAlg)(double**, const int, double*),
    double**  corr_matr,
    double**  sqrt_matr,
    const int n);

void CholAlg(double** a, const int n, double* p);

#endif