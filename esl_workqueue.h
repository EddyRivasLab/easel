/* Threaded work queue.
 * 
 * SVN $Id: $
 */
#ifndef ESL_WORKQUEUE_INCLUDED
#define ESL_WORKQUEUE_INCLUDED

/* ESL_WORKQUEUE */
typedef struct ESL_WORK_QUEUE_ST ESL_WORK_QUEUE;

ESL_WORK_QUEUE *esl_workqueue_Create(int size);
void esl_workqueue_Destroy(ESL_WORK_QUEUE *queue);

int esl_workqueue_Init(ESL_WORK_QUEUE *queue, void *ptr);
int esl_workqueue_Complete(ESL_WORK_QUEUE *queue);
int esl_workqueue_Reset(ESL_WORK_QUEUE *queue);

void * esl_workqueue_Remove(ESL_WORK_QUEUE *queue);

int esl_workqueue_ReaderUpdate(ESL_WORK_QUEUE *queue, void *in, void **out);
int esl_workqueue_WorkerUpdate(ESL_WORK_QUEUE *queue, void *in, void **out);

int esl_workqueue_Dump(ESL_WORK_QUEUE *queue);

#endif /*ESL_WORKQUEUE_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
