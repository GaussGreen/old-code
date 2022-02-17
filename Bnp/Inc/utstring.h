/* ======================================================
   FILENAME:  utstring.h

   PURPOSE:   A few simple functions when dealing with
              a string
   ====================================================== */

#ifndef UTSTRING_H
#define UTSTRING_H

/* --------------------------------------------------------------------------
   VERY USEFUL:   the Tring type is a char *
   -------------------------------------------------------------------------- */

typedef char *String;

void strupper(char *s);
void strip_white_space(char *s);

char *new_string(char *s);

void add_tick_string(char *s, long ticker, char *s_tick);
void rem_tick_string(char *s_tick, char *s);

#endif
