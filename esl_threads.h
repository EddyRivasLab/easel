/* Threaded work queue.
 * 
 * SVN $Id$
 */
#ifndef ESL_THREADS_INCLUDED
#define ESL_THREADS_INCLUDED

/* ESL_THREADS */
typedef struct ESL_THREADS_ST ESL_THREADS;

typedef void *(*ESL_THREAD_FUNCTION)(void *);

ESL_THREADS *esl_threads_Create(ESL_THREAD_FUNCTION func);
void esl_threads_Destroy(ESL_THREADS *obj);

int esl_threads_AddThread(ESL_THREADS *obj, void *data);
int esl_threads_WaitForStart(ESL_THREADS *obj);
int esl_threads_WaitForFinish(ESL_THREADS *obj);

int   esl_threads_Started(ESL_THREADS *obj);
void *esl_threads_GetData(ESL_THREADS *obj);
void  esl_threads_Exit(ESL_THREADS *obj);


#endif /*ESL_THREADS_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
