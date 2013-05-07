#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

/*
 *  Priority queue structure
 */
struct pqueue 
{
	int size, avail, step;
	PQDATUM *d;
};

/*
 *  pqinit: initialize the queue.
 *
 *  Parameters:
 *
 *    q           Pointer to a priority queue, or NULL if the user
 *                wishes to leave it to pqinit to allocate the queue.
 *
 *    n           Numer of queue items for which memory should be
 *                preallocated, that is, the initial size of the
 *                item array the queue uses. If you insert more than
 *                n items to the queue, another n items will
 *                be allocated automatically.
 *
 *  Return values:
 *
 *   non-NULL     Priority queue has been initialized.
 *
 *   NULL         Insufficient memory.
 */
static struct pqueue *
pqinit(struct pqueue *q, int n)
{
	struct pqueue *tmp = q;

	if (!q && !(q = malloc(sizeof(struct pqueue)))) {
		return NULL;
	}
	if (!(q->d = malloc(sizeof(PQDATUM) * n))) {
		if (!tmp) free(q);
		return NULL;
	}
	q->avail = q->step = n;
	q->size = 1;
	return q;
}

/*                  
 *  pqinsert: insert an item into the queue.
 *
 *  Parameters:
 *
 *    q           Pointer to a priority queue.
 *
 *    d           Datum to be inserted.
 *
 *  Return values:
 *
 *    1           The item has been inserted.
 *
 *    0           The item could not be appended. Either the queue i
 *                pointer provided was NULL, or the function was unable 
 *                to allocate the amount of memory needed for 
 *                the new item.
 */
static int
pqinsert(struct pqueue *q, PQDATUM d)
{
	PQDATUM *tmp;
	int i, newsize;

	if (!q) return 0;
	
	/* allocate more memory if necessary */
	if (q->size >= q->avail) {
		newsize = q->size + q->step;
		if (!(tmp = realloc(q->d, sizeof(PQDATUM) * newsize))) {
			return 0;
		};
		q->d = tmp;
		q->avail = newsize;		
	}

	/* insert item */
	i = q->size++;
	while (i > 1 && PQPRIO(q->d[i / 2]) < PQPRIO(d)) {
		q->d[i] = q->d[i / 2];
		i /= 2;
	}
	q->d[i] = d;
	return 1;	
} 

/*
 *  pqremove: remove the highest-ranking item from the queue.
 *
 *  Parameters:
 *
 *    p           Pointer to a priority queue.
 *
 *    d           Pointer to the PQDATUM variable that will hold the 
 *                datum corresponding to the queue item removed.               
 *
 *  Return values:
 *
 *    non-NULL    An item has been removed. The variable that d points
 *                to now contains the datum associated with the item
 *                in question.
 *
 *    NULL        No item could be removed. Either the queue pointer
 *                provided was NULL, or the queue was empty. The chunk
 *                of memory that d points to has not been modified.
 */
static PQDATUM *
pqremove(struct pqueue *q, PQDATUM *d)
{	
	PQDATUM tmp;
	int i = 1, j;

	if (!q || q->size == 1) return NULL;
	*d = q->d[1];
	tmp = q->d[--q->size];
	while (i <= q->size / 2) {
		j = 2 * i;
		if (j < q->size && 
			PQPRIO(q->d[j]) < PQPRIO(q->d[j + 1])) {
			j++;
		}
		if (PQPRIO(q->d[j]) <= PQPRIO(tmp)) {
			break;
		}
		q->d[i] = q->d[j];
		i = j;
	}
	q->d[i] = tmp;
	return d;	
} 

/*
 *  pqpeek: access highest-ranking item without removing it.
 *
 *  Parameters:
 *
 *    q           Pointer to a priority queue.
 *
 *    d           Pointer to the PQDATUM variable that will hold the
 *                datum corresponding to the highest-ranking item.
 *                
 *  Return values:
 *
 *    non-NULL   Success. The variable that d points to now contains
 *               the datum associated with the highest-ranking item.
 *
 *    NULL       Failure. Either the queue pointer provided was NULL,
 *               or the queue was empty. The chunk of memory that d
 *               points to has not been modified.
 */
static PQDATUM *
pqpeek(struct pqueue *q, PQDATUM *d)
{
	if (!q || q->size == 1) return NULL;
	*d = q->d[1];
	return d;
}

