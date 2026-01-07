/* RNA secondary structure markup in WUSS notation.
 */
#ifndef eslWUSS_INCLUDED
#define eslWUSS_INCLUDED
#include <esl_config.h>

extern int esl_wuss_IsSimple(const char *ss, int n);
extern int esl_wuss_IsFull  (const char *ss, int n);
extern int esl_wuss2ct      (const char *ss, int n, int  *ct);
extern int esl_ct2wuss      (const int  *ct, int n, char *ss);
extern int esl_ct2simplewuss(const int  *ct, int n, char *ss);

extern int esl_wuss_full    (const char *oldss, char *newss);
extern int esl_wuss_nopseudo(const char *ss1,   char *ss2);
extern int esl_wuss_reverse (const char *ss,    char *new);

extern int esl_wuss2kh(const char *ss, char *kh);
extern int esl_kh2wuss(const char *kh, char *ss);


#endif /*eslWUSS_INCLUDED*/

