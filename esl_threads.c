/* Basic threads package
 * 
 * Contents:
 *    1. Work queue routines
 *   16. Copyright and license.
 * 
 * MSF, Thu Jun 18 11:51:39 2009
 * SVN $Id$
 */
#include "esl_config.h"

#ifdef HAVE_PTHREADS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "easel.h"
#include "esl_threads.h"

/* Internal structures */

typedef struct THREAD_INFO_ST THREAD_INFO;

struct THREAD_INFO_ST {
  pthread_t        threadId;
  void            *data;
  ESL_THREADS     *obj;
  THREAD_INFO     *next;
};

struct ESL_THREADS_ST {
  int                    threadCount;

  int                    startThread;
  pthread_mutex_t        startMutex;
  pthread_cond_t         startCond;

  THREAD_INFO           *threads;

  ESL_THREAD_FUNCTION    func;
};

#define CHECK(cond,fmt)		\
  if (cond) esl_fatal("(%s:%d) - " fmt ": %s (%d)", \
		      __FUNCTION__, __LINE__, \
		      strerror (cond), cond)

/*****************************************************************
 *# 1. Thread routines
 *****************************************************************/ 

/* Function:  esl_threads_Create()
 * Synopsis:  Create a threads object.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Creates an <ESL_THREADS> object.
 *
 * Returns:   ptr to the new <ESL_THREADS> object.
 */
ESL_THREADS *esl_threads_Create(ESL_THREAD_FUNCTION fnptr)
{
  int status;

  ESL_THREADS *obj;

  ESL_ALLOC(obj, sizeof(*obj));

  obj->threadCount     = 0;
  obj->startThread     = 0;
  obj->threads         = NULL;

  obj->func            = fnptr;

  status = pthread_mutex_init (&obj->startMutex, NULL);
  CHECK (status, "Create mutex failed");

  status = pthread_cond_init (&obj->startCond, NULL);
  CHECK (status, "Create cond var failed");

  return obj;

 ERROR:
  esl_fatal("Could not allocate threads object");
  return NULL;
}

/* Function:  esl_threads_Destroy()
 * Synopsis:  Destroys an <ESL_THREADS> object.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Frees an <ESL_THREADS> object.  
 *
 *            The calling routine is responsible for freeing the
 *            memory of the thread data.
 *
 * Returns:   void
 */
void esl_threads_Destroy(ESL_THREADS *obj)
{
  int           status;
  THREAD_INFO  *thread;

  if (obj != NULL)
    {
      status = pthread_mutex_destroy (&obj->startMutex);
      CHECK (status, "Destroy mutex failed");

      status = pthread_cond_destroy (&obj->startCond);
      CHECK (status, "Destroy cond var failed");

      thread = obj->threads;
      while (thread != NULL)
	{
	  THREAD_INFO  *tmp = thread;
	  thread = thread->next;

	  free(tmp);
	}

      free(obj);
    }
}

/* Function:  esl_threads_AddThread()
 * Synopsis:  Adds a thread to the <ESL_THREADS> object.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Create a new thread for the <ESL_THREADS> object.
 *
 *            The calling routine is responsible for freeing the
 *            memory of the thread data.
 *
 * Returns:   void
 */
int esl_threads_AddThread(ESL_THREADS *obj, void *data)
{
  int           status;
  THREAD_INFO  *thread;

  if (obj == NULL) esl_fatal("Invalid thread object");

  ESL_ALLOC(thread, sizeof(*thread));

  thread->data   = data;
  thread->obj    = obj;
  thread->next   = obj->threads;

  obj->threads   = thread;
  ++obj->threadCount;

  status = pthread_create(&thread->threadId, NULL, obj->func, obj);
  CHECK(status, "Create thread failed");

  return eslOK;

 ERROR:
  esl_fatal("Could not allocate thread");
  return eslFAIL;
}

/* Function:  esl_threads_WaitForStart()
 * Synopsis:  Blocks until all threads have started.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Block until all the threads have started.  When all the threads
 *            have started and are blocking at the start mutex, release
 *            them.
 *
 * Returns:   void
 */
int esl_threads_WaitForStart(ESL_THREADS *obj)
{
  int           status;

  if (obj == NULL) esl_fatal("Invalid thread object");

  status = pthread_mutex_lock (&obj->startMutex);
  CHECK(status, "Lock mutex failed");

  /* wait for the threads to have started */
  while (obj->startThread < obj->threadCount) {
    status = pthread_cond_wait(&obj->startCond, &obj->startMutex);
    CHECK(status, "Wait cond failed");
  }

  /* release all the threads */
  obj->startThread = 0;
  status = pthread_cond_broadcast (&obj->startCond);
  CHECK(status, "Cond broadcast failed");

  status = pthread_mutex_unlock (&obj->startMutex);
  CHECK(status, "Unlock mutex failed");

  return obj->threadCount;
}

/* Function:  esl_threads_WaitForFinish()
 * Synopsis:  Blocks until all threads have complted.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Block until all the threads have completed.  Remove the
 *            thread object from the list as it completes.
 *
 * Returns:   void
 */
int esl_threads_WaitForFinish(ESL_THREADS *obj)
{
  int           status;
  THREAD_INFO   *thread;

  if (obj == NULL) esl_fatal("Invalid thread object");

  /* wait for the threads to complete */
  thread = obj->threads;
  while (thread != NULL)
    {
      THREAD_INFO  *tmp = thread;
      thread = thread->next;

      status = pthread_join(tmp->threadId, NULL);
      CHECK(status, "Join thread failed");

      status = pthread_mutex_lock (&obj->startMutex);
      CHECK(status, "Lock mutex failed");

      obj->threads = thread;
      obj->threadCount--;

      status = pthread_mutex_unlock (&obj->startMutex);
      CHECK(status, "Unlock mutex failed");

      free(tmp);
    }

  return eslOK;
}

/* Function:  esl_threads_Started()
 * Synopsis:  Blocks until all threads have started.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Block until all the threads have started.  When all the threads
 *            have started and are blocking at the start mutex, release
 *            them.
 *
 * Returns:   void
 */
int esl_threads_Started(ESL_THREADS *obj)
{
  int           status;
  pthread_t     threadId;

  THREAD_INFO  *thread;

  if (obj == NULL) esl_fatal("Invalid thread object");

  threadId = pthread_self();

  status = pthread_mutex_lock (&obj->startMutex);
  CHECK(status, "Lock mutex failed");

  /* make sure the thread has registered */
  thread = obj->threads;
  while (thread != NULL && thread->threadId != threadId)
    {
      thread = thread->next;
    }

  if (thread == NULL) esl_fatal("Thread has not registed");

  /* signal that we have started */
  ++obj->startThread;
  status = pthread_cond_broadcast (&obj->startCond);
  CHECK(status, "Cond broadcast failed");

  /* wait for the signal to start the calculations */
  while (obj->startThread) {
    status = pthread_cond_wait(&obj->startCond, &obj->startMutex);
    CHECK(status, "Cond wait failed");
  }

  status = pthread_mutex_unlock (&obj->startMutex);
  CHECK(status, "Unlock mutex failed");

  return eslOK;
}

/* Function:  esl_threads_GetData()
 * Synopsis:  Return the data associated with this thread.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Return the data associated with this thread.  The data for
 *            this thread is set in the esl_threads_AddThread() function.
 *
 * Returns:   void *
 */
void *esl_threads_GetData(ESL_THREADS *obj)
{
  int           status;
  pthread_t     threadId;

  THREAD_INFO  *thread;

  if (obj == NULL) esl_fatal("Invalid thread object");

  threadId = pthread_self();

  status = pthread_mutex_lock (&obj->startMutex);
  CHECK(status, "Lock mutex failed");

  /* make sure the thread has registered */
  thread = obj->threads;
  while (thread != NULL && thread->threadId != threadId)
    {
      thread = thread->next;
    }

  status = pthread_mutex_unlock (&obj->startMutex);
  CHECK(status, "Unlock mutex failed");

  if (thread == NULL) esl_fatal("Thread has not registed");

  return thread->data;
}

/* Function:  esl_threads_Exit()
 * Synopsis:  Terminate the thread.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Terminate the thread.
 *
 * Returns:   void
 */
void esl_threads_Exit(ESL_THREADS *obj)
{
  int           status;
  pthread_t     threadId;

  THREAD_INFO  *thread;

  if (obj == NULL) esl_fatal("Invalid thread object");

  threadId = pthread_self();

  status = pthread_mutex_lock (&obj->startMutex);
  CHECK(status, "Lock mutex failed");

  /* make sure the thread has registered */
  thread = obj->threads;
  while (thread != NULL && thread->threadId != threadId)
    {
      thread = thread->next;
    }

  status = pthread_mutex_unlock (&obj->startMutex);
  CHECK(status, "Unlock mutex failed");

  if (thread == NULL) esl_fatal("Thread has not registed");

  pthread_exit (NULL);
  return;
}

#endif /* HAVE_PTHREADS */

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
