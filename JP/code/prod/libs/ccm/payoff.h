/* $Header$ */
#include "creditPortfolio.h"

/* Spread Curve structure */
typedef struct
{
    long nbTimes;
    double *time;                   /* [nbTimes] */
    double *spread;                 /* [nbTimes] */
} SPREAD_CURVE;

/* Uniform Sim structure */
typedef struct
{
    long nbNames;
    long nbPaths;
    double *uniformDeviate;         /* [nbNames*nbPath] */
    double *weight;                 /* [nbPaths] */
} UNIFORM_SIM;

/* Scenario SIm structure */
typedef struct 
{
    long nbNames;
    long nbPaths;
    double *defaultTime;            /* [nbNames*nbPath] */
    double *weight;                 /* [nbPath] */
} SCENARIO_SIM;

double *CumulatedSpreadToFlatForward(
        double *time,               /* (I) [nbTimes] */
        double *cumSpread,          /* (I) [nbTimes] */
        long nbTimes);

UNIFORM_SIM UniformSimCreate(
        long nbPaths,
        long nbNames,
        int copulaType,
        const double *beta,         /* (I) [nbNames] */
        long freedomDegree);

void UniformSimFree(UNIFORM_SIM *us);

SPREAD_CURVE SpreadCurveCreate(
        double *time,
        double *spread,
        long nbTimes);

void SpreadCurveFree(SPREAD_CURVE *sc);

SCENARIO_SIM ScenarioSimCreate(
        UNIFORM_SIM *us,            /* (I) */
        SPREAD_CURVE *sc);          /* (I) [nbNames] */

SCENARIO_SIM ScenarioSimCreate2(
        UNIFORM_SIM *us1,           /* (I) */
        UNIFORM_SIM *us2,           /* (I) */
        SPREAD_CURVE *sc,           /* (I) [nbNames] */
        const double *alpha);       /* (I) [nbNames] */

void ScenarioSimFree(SCENARIO_SIM *ss);


double *ExpectedPayoff(
        CREDIT_PORTFOLIO *port,     /* (I) */
        SCENARIO_SIM     *s,        /* (I) */
        const double     *samplingTime,/* (I) [nbSample] */
        long              nbSample);

UNIFORM_SIM CopulaProduct(
        UNIFORM_SIM u_s1,
        UNIFORM_SIM u_s2,
        const double *alpha);       /* (I) [nbNames] */

double DefaultTimeFromUniform(
        const double *time,         /* (I) [nbTimes] */
        const double *spread,       /* (I) [nbTimes] */
        long nbTimes,
        double triggerLevel);

