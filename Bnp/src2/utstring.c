/* ===============================================================
   FILE_NAME:	utstring.c

   PURPOSE:     A few functions to work with a string
   =============================================================== */

#include "utallhdr.h"

void strupper(char* s)
{
    while (*s = toupper(*s))
        s++;
}

void strip_white_space(char* a)
{
    char* b = a;

    while (*a != '\0')
    {
        if (!isspace(*a))
            *b++ = *a++;
        else
            a++;
    }
    *b = '\0';
}

/**********
        allocate space for a string s,
        copy s to that space, return new space
**********/

char* new_string(char* s)
{
    char* news;
    int   len = strlen(s);
    news      = (char*)srt_calloc(len + 1, sizeof(char));
    if (!news)
        return NULL;
    strcpy(news, s);
    news[len] = '\0';
    return news;
}

/* Add a ticker to a string */
void add_tick_string(char* s, long ticker, char* s_tick)
{
    strcpy(s_tick, s);
    strcat(s_tick, ".");
    sprintf(s_tick, "%s%d\0", s_tick, ticker);
}

/* Remove a ticker into a string */
void rem_tick_string(char* s_tick, char* s)
{
    char* dot_pos;

    strcpy(s, s_tick);

    dot_pos = strrchr(s, '.');
    if (dot_pos != NULL)
        (*dot_pos) = '\0';
}
