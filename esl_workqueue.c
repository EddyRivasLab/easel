/* Threaded work queue.
 * 
 * Contents:
 *    1. Work queue routines
 *   16. Copyright and license.
 * 
 * MSF, Thu Jun 18 11:51:39 2009
 * SVN $Id: $
 */
#include "esl_config.h"

#ifdef HAVE_PTHREADS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "easel.h"
#include "esl_workqueue.h"


/* Internal structures */
struct ESL_WORK_QUEUE_ST {
  pthread_mutex_t  queueMutex;
  pthread_cond_t   readerQueueCond;
  pthread_cond_t   workerQueueCond;

  void           **readerQueue;
  int              readerQueueCnt;
  int              readerQueueHead;

  void           **workerQueue;
  int              workerQueueCnt;
  int              workerQueueHead;

  int              queueSize;
  int              pendingWorkers;
};

#define CHECK(cond,fmt)		\
  if (cond) esl_fatal("(%s:%d) - " fmt ": %s (%d)", \
		      __FUNCTION__, __LINE__, \
		      strerror (cond), cond)

/*****************************************************************
 *# 1. Work queue routines
 *****************************************************************/ 

/* Function:  esl_workqueue_Create()
 * Synopsis:  Create a work queue object.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Creates an <ESL_WORK_QUEUE> object of <size>.
 *
 * Returns:   ptr to the new <ESL_WORK_QUEUE> object.
 */
ESL_WORK_QUEUE *esl_workqueue_Create(int size)
{
  int i;
  int status;

  ESL_WORK_QUEUE *queue;

  ESL_ALLOC(queue, sizeof(*queue));

  queue->readerQueue     = NULL;
  queue->readerQueueCnt  = 0;
  queue->readerQueueHead = 0;

  queue->workerQueue     = NULL;
  queue->workerQueueCnt  = 0;
  queue->workerQueueHead = 0;

  queue->queueSize       = size;
  queue->pendingWorkers  = 0;

  status = pthread_mutex_init (&queue->queueMutex, NULL);
  CHECK (status, "Create mutex failed");

  status = pthread_cond_init (&queue->readerQueueCond, NULL);
  CHECK (status, "Create cond var failed");

  status = pthread_cond_init (&queue->workerQueueCond, NULL);
  CHECK (status, "Create cond var failed");

  ESL_ALLOC(queue->readerQueue, sizeof(void *) * size);
  ESL_ALLOC(queue->workerQueue, sizeof(void *) * size);

  for (i = 0; i < queue->queueSize; ++i)
    {
      queue->readerQueue[i] = NULL;
      queue->workerQueue[i] = NULL;
    }

  return queue;

 ERROR:
  esl_fatal("Could not allocate work queue");
  return NULL;
}

/* Function:  esl_workqueue_Destroy()
 * Synopsis:  Destroys an <ESL_WORK_QUEUE> object.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Frees an <ESL_WORK_QUEUE> object.  
 *
 *            The calling routine is responsible for freeing the
 *            memory of the actual queued objects.
 *
 * Returns:   void
 */
void esl_workqueue_Destroy(ESL_WORK_QUEUE *queue)
{
  int status;

  if (queue != NULL)
    {
      status = pthread_mutex_destroy (&queue->queueMutex);
      CHECK (status, "Destroy mutex failed");

      status = pthread_cond_destroy (&queue->readerQueueCond);
      CHECK (status, "Destroy cond var failed");

      status = pthread_cond_destroy (&queue->workerQueueCond);
      CHECK (status, "Destroy cond var failed");

      if (queue->readerQueue != NULL) free(queue->readerQueue);
      if (queue->workerQueue != NULL) free(queue->workerQueue);

      free(queue);
    }
}

/* Function:  esl_workqueue_Init()
 * Synopsis:  Adds a queued object to the producers list.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Added a queued object to the producers list checking for
 *            any errors.
 *
 * Returns:   <eslOK> on success.
 */
int esl_workqueue_Init(ESL_WORK_QUEUE *queue, void *ptr)
{
  int cnt;
  int inx;
  int status;

  int queueSize;

  if (queue == NULL) return eslFAIL;

  if (ptr == NULL) esl_fatal("Placing NULL object in reader queue");

  status = pthread_mutex_lock (&queue->queueMutex);
  CHECK(status, "Lock mutex failed");

  queueSize = queue->queueSize;

  /* check to make sure we won't overflow */
  cnt = queue->readerQueueCnt;
  if (cnt >= queueSize) esl_fatal("Reader queue overflow");

  inx = (queue->readerQueueHead + cnt) % queueSize;
  queue->readerQueue[inx] = ptr;

  ++queue->readerQueueCnt;
  if (cnt == 0)
    {
      status = pthread_cond_signal (&queue->readerQueueCond);
      CHECK (status, "Condition signal failed");
    }

  status = pthread_mutex_unlock (&queue->queueMutex);
  CHECK(status, "Unlock mutex failed");

  return eslOK;
}

/* Function:  esl_workqueue_Remove()
 * Synopsis:  Removes a queued object from the producers list.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Removes a queued object to the producers list.
 *
 * Returns:   Pointer to the object removed.  NULL if the producers
 *            list is empty.
 */
void *esl_workqueue_Remove(ESL_WORK_QUEUE *queue)
{
  int inx;
  int status;

  void *ptr = NULL;

  if (queue == NULL) return eslFAIL;

  status = pthread_mutex_lock (&queue->queueMutex);
  CHECK(status, "Lock mutex failed");

  /* check if there are any items on the readers list */
  if (queue->readerQueueCnt > 0)
    {
      inx = (queue->readerQueueHead + queue->readerQueueCnt) % queue->queueSize;
      ptr = queue->readerQueue[inx];
      queue->readerQueue[inx] = NULL;
      --queue->readerQueueCnt;
    }

  status = pthread_mutex_unlock (&queue->queueMutex);
  CHECK(status, "Unlock mutex failed");

  return ptr;
}

/* Function:  esl_workqueue_Complete()
 * Synopsis:  Signals the end of the queue.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Signal the end of the queue.  If there are any threads
 *            waiting on an object, signal them to wake up and complete
 *            their processing.
 *
 * Returns:   <eslOK> on success.
 */
int esl_workqueue_Complete(ESL_WORK_QUEUE *queue)
{
  int status;

  if (queue == NULL) return eslFAIL;

  status = pthread_mutex_lock (&queue->queueMutex);
  CHECK(status, "Lock mutex failed");

  if (queue->pendingWorkers != 0)
    {
      status = pthread_cond_broadcast (&queue->workerQueueCond);
      CHECK (status, "Condition broadcast failed");
    }

  status = pthread_mutex_unlock (&queue->queueMutex);
  CHECK(status, "Unlock mutex failed");

  return eslOK;
}

/* Function:  esl_workqueue_Reset()
 * Synopsis:  Reset the queue for another run.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Reset the queue for another run.  This is done by moving
 *            all the queued object to the reader's list (i.e. producer).
 *
 * Returns:   <eslOK> on success.
 */
int esl_workqueue_Reset(ESL_WORK_QUEUE *queue)
{
  int inx;
  int status;

  int queueSize;

  if (queue == NULL) return eslFAIL;

  status = pthread_mutex_lock (&queue->queueMutex);
  CHECK(status, "Lock mutex failed");

  queueSize = queue->queueSize;

  /* move all buffers back to the reader queue */
  while (queue->workerQueueCnt > 0) 
    {
      inx = (queue->readerQueueHead + queue->readerQueueCnt) % queueSize;
      queue->readerQueue[inx] = queue->workerQueue[queue->workerQueueHead];
      ++queue->readerQueueCnt;

      queue->workerQueue[queue->workerQueueHead] = NULL;
      queue->workerQueueHead = (queue->workerQueueHead + 1) % queueSize;
      --queue->workerQueueCnt;
    }

  queue->pendingWorkers = 0;

  status = pthread_mutex_unlock (&queue->queueMutex);
  CHECK(status, "Unlock mutex failed");

  return eslOK;
}

/* Function:  esl_workqueue_Dump()
 * Synopsis:  Print the contents of the queues.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   Print the contents of the queues and their pointers.
 *
 * Returns:   <eslOK> on success.
 */
int esl_workqueue_Dump(ESL_WORK_QUEUE *queue)
{
  int i;
  int status;

  if (queue == NULL) return eslFAIL;

  status = pthread_mutex_lock (&queue->queueMutex);
  CHECK(status, "Lock mutex failed");

  printf ("Reader head: %2d  count: %2d\n", queue->readerQueueHead, queue->readerQueueCnt);
  printf ("Worker head: %2d  count: %2d\n", queue->workerQueueHead, queue->workerQueueCnt);
  for (i = 0; i < queue->queueSize; ++i)
    {
      printf ("  %2d:  %p  %p\n", i, queue->readerQueue[i], queue->workerQueue[i]);
    }
  printf ("Pending: %2d\n\n", queue->pendingWorkers);

  status = pthread_mutex_unlock (&queue->queueMutex);
  CHECK(status, "Unlock mutex failed");

  return eslOK;
}

/* Function:  esl_workqueue_ReaderUpdate()
 * Synopsis:  Producer routine.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   The producer (i.e. Reader) places an object on the
 *            consumers (i.e. Workers) queue.
 *
 *            If the <in> object is not null, it is placed on the
 *            workers queue.  If there are any workers waiting for
 *            an object, signal them to wake up.
 *
 *            If the calling routine has supplied an <out> pointer,
 *            an object pointer is read off of the readers queue and
 *            made available to the caller.
 *
 * Returns:   <eslOK> on success.
 */
int esl_workqueue_ReaderUpdate(ESL_WORK_QUEUE *queue, void *in, void **out)
{
  int inx;
  int status;

  int queueSize;

  if (queue == NULL) return eslFAIL;

  status = pthread_mutex_lock (&queue->queueMutex);
  CHECK(status, "Lock mutex failed");

  queueSize = queue->queueSize;

  /* check if the caller is queuing up an item */
  if (in != NULL)
    {

      /* check to make sure we don't overflow */
      if (queue->workerQueueCnt >= queue->queueSize)
	{
	  esl_fatal("Work queue overflow");
	}

      inx = (queue->workerQueueHead + queue->workerQueueCnt) % queueSize;
      queue->workerQueue[inx] = in;
      ++queue->workerQueueCnt;

      if (queue->pendingWorkers != 0)
	{
	  status = pthread_cond_broadcast (&queue->workerQueueCond);
	  CHECK (status, "Condition broadcast failed");
	}
    }

  /* check if the caller is waiting for a queued item */
  if (out != NULL)
    {

      /* wait for a processed buffers to be returned */
      while (queue->readerQueueCnt == 0) 
	{
	  status = pthread_cond_wait (&queue->readerQueueCond, 
				      &queue->queueMutex);
	  CHECK (status, "Condition wait failed");
	}

      inx = queue->readerQueueHead;
      *out = queue->readerQueue[inx];
      queue->readerQueue[inx] = NULL;
      queue->readerQueueHead = (queue->readerQueueHead + 1) % queueSize;
      --queue->readerQueueCnt;
    }

  status = pthread_mutex_unlock (&queue->queueMutex);
  CHECK(status, "Unlock mutex failed");

  return eslOK;
}

/* Function:  esl_workqueue_WorkerUpdate()
 * Synopsis:  Consumer routine.
 * Incept:    MSF, Thu Jun 18 11:51:39 2009
 *
 * Purpose:   The consumer (i.e. Workers) places an object on the
 *            producers (i.e. Readers) queue.
 *
 *            If the <in> object is not null, it is placed on the
 *            readers queue.  If the reader is waiting for an object, 
 *            signal it to wake up.
 *
 *            If the calling routine has supplied an <out> pointer,
 *            an object pointer is read off of the workers queue and
 *            made available to the caller.
 *
 * Returns:   <eslOK> on success.
 */
int esl_workqueue_WorkerUpdate(ESL_WORK_QUEUE *queue, void *in, void **out)
{
  int cnt;
  int inx;
  int status;

  int queueSize;

  if (queue == NULL) return eslFAIL;

  status = pthread_mutex_lock (&queue->queueMutex);
  CHECK(status, "Lock mutex failed");

  queueSize = queue->queueSize;

  /* check if the caller is queuing up an item */
  if (in != NULL)
    {

      /* check to make sure we don't overflow */
      if (queue->readerQueueCnt >= queue->queueSize)
	{
	  esl_fatal("Reader queue overflow");
	}

      inx = (queue->readerQueueHead + queue->readerQueueCnt) % queueSize;
      queue->readerQueue[inx] = in;
      cnt = queue->readerQueueCnt++;
      if (cnt == 0)
	{
	  status = pthread_cond_signal (&queue->readerQueueCond);
	  CHECK (status, "Condition signal failed");
	}
    }

  /* check if the caller is waiting for a queued item */
  if (out != NULL)
    {

      if (queue->workerQueueCnt == 0)
	{
	  /* wait for a processed buffers to be returned */
	  ++queue->pendingWorkers;
	  while (queue->workerQueueCnt == 0)
	    {
	      status = pthread_cond_wait (&queue->workerQueueCond, 
					  &queue->queueMutex);
	      CHECK (status, "Condition wait failed");
	    }
	  --queue->pendingWorkers;
	}

      inx = queue->workerQueueHead;
      *out = queue->workerQueue[inx];
      queue->workerQueue[inx] = NULL;
      queue->workerQueueHead = (queue->workerQueueHead + 1) % queueSize;
      --queue->workerQueueCnt;
    }

  status = pthread_mutex_unlock (&queue->queueMutex);
  CHECK(status, "Unlock mutex failed");

  return eslOK;
}

#endif /* HAVE_PTHREADS */

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
