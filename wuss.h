/* wuss.h
 * RNA secondary structure markup in WUSS notation.
 * 
 * SVN $Id$
 * SRE, Tue Feb 15 10:15:28 2005
 */
#ifndef eslWUSS_INCLUDED
#define eslWUSS_INCLUDED


extern int esl_wuss2ct(char *ss, int len, int *ct);
extern int esl_wuss2kh(char *ss, char *kh);
extern int esl_kh2wuss(char *kh, char *ss);
extern int esl_wuss_nopseudo(char *ss1, char *ss2);


#endif /*eslWUSS_INCLUDED*/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/
