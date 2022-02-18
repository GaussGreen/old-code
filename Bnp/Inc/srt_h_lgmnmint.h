#ifndef SRT_H_LGMNMINT_H
#define SRT_H_LGMNMINT_H

SrtErr srt_f_lgmrfltr(
    Date      start1,
    Date      enfp1,
    char*     freq,
    char*     bas,
    Date      start2,
    Date      enfp2,
    Date      fix1,
    Date      fix2,
    Date      paydt,
    double    strike,
    double    strike2,
    double    margin,
    SrtUndPtr und,
    double    rate1,
    double    rate2,
    double*   answer);

#endif
