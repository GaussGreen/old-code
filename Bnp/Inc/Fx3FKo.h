
#ifndef Fx3FKo_h
#define Fx3FKo_h

Err Fx3FKoDig(
char					*undName,
int						nKo,
long					*koDates,
double					*koLvl,
int						*ab,		/*	1: Above, -1: Below */
long					payDate,
long					npth,
double					*ans);

Err Fx3FKoFwd(
char					*undName,
int						nKo,
long					*koDates,
double					*koLvl,
int						*ab,		/*	1: Above, -1: Below */
long					fixDate,
long					payDate,
long					npth,
double					*ans);

Err Fx3FKoOption(
char					*undName,
int						nKo,
long					*koDates,
double					*koLvl,
int						*ab,		/*	1: Above, -1: Below */
double					strike,
int						callPut,	/*	1: Call, -1: Put */
long					fixDate,
long					payDate,
long					npth,
double					*ans);

#endif