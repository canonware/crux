/******************************************************************************
 *
 * <Copyright = jasone>
 * <License>
 *
 ******************************************************************************
 *
 * Version: Crux <Version = crux>
 *
 ******************************************************************************/

#include "../include/_cruxmodule.h"

static cw_bhpi_t *
bhp_p_dump(cw_bhpi_t *a_bhpi, uint32_t a_depth, cw_bhpi_t *a_last_printed)
{
    uint32_t i;

    /* Sibling. */
    if (a_bhpi->sibling != NULL)
    {
	a_last_printed = bhp_p_dump(a_bhpi->sibling, a_depth, a_last_printed);
    }
    /* Self. */
    if (a_bhpi->parent != a_last_printed)
    {
	/* Indent. */
	for (i = 0; i < (a_depth * 38); i++)
	{
	    fprintf(stderr, " ");
	}
    }
    fprintf(stderr, "[deg:%u pri:0x%8p dat:0x%8p]",
	    a_bhpi->degree, a_bhpi->priority, a_bhpi->data);
    a_last_printed = a_bhpi;

    /* Child. */
    if (a_bhpi->child != NULL)
    {
	fprintf(stderr, "-");
	a_last_printed = bhp_p_dump(a_bhpi->child, a_depth + 1, a_bhpi);
    }
    else
    {
	fprintf(stderr, "\n");
    }

    return a_last_printed;
}

/******************************************************************************
 *
 * Link two binomial heaps of the same degree (n) together into one heap
 * of degree (n + 1).  a_root points to the root of the resulting heap.
 *
 ******************************************************************************/
static void
bhp_p_bin_link(cw_bhpi_t *a_root, cw_bhpi_t *a_non_root)
{
    a_non_root->parent = a_root;
    a_non_root->sibling = a_root->child;
    a_root->child = a_non_root;
    a_root->degree++;
}

/******************************************************************************
 *
 * Merge the root lists of the two heaps specified by the arguments, in
 * monotonically increasing order.  The result is stored in a_a.
 *
 ******************************************************************************/
CW_P_INLINE void
bhp_p_merge(cw_bhp_t *a_a, cw_bhp_t *a_b)
{
    if (a_a->head == NULL)
    {
	a_a->head = a_b->head;
    }
    else if (a_b->head != NULL)
    {
	cw_bhpi_t *mark_a, *curr_a, *mark_b, *curr_b;

	/* Both heaps have contents. */

	if (a_a->head->degree > a_b->head->degree)
	{
	    cw_bhpi_t *t_bhpi;

	    /* Swap the heads to simplify the following loop.  Now we know that
	     * a_a->head is set correctly. */
	    t_bhpi = a_a->head;
	    a_a->head = a_b->head;
	    a_b->head = t_bhpi;
	}
	mark_a = NULL; /* Avoid optimization warnings about uninitialized
			* reads. */
	curr_a = a_a->head;
	curr_b = a_b->head;
	while ((curr_a->sibling != NULL) && (curr_b != NULL))
	{
	    /* Fast forward to where we need to insert from b. */
	    while ((curr_a->sibling != NULL)
		   && (curr_a->degree <= curr_b->degree))
	    {
		mark_a = curr_a;
		curr_a = curr_a->sibling;
	    }

	    /* Move forward in b. */
	    mark_b = curr_b;
	    curr_b = curr_b->sibling;

	    /* Link things together. */
	    mark_a->sibling = mark_b;
	    while ((curr_b != NULL) && (curr_b->degree <= curr_a->degree))
	    {
		mark_b = curr_b;
		curr_b = curr_b->sibling;
	    }
	    mark_a = mark_b;
	    mark_b->sibling = curr_a;
	}

	/* If there are still nodes in b, append them. */
	if (curr_b != NULL)
	{
	    /* curr_a->sibling is implicitly NULL if this is true, due to the
	     * loop exit condition above. */
	    cw_assert(curr_a->sibling == NULL);
	    cw_check_ptr(curr_b);
	    curr_a->sibling = curr_b;
	}
    }
    /* Adjust the size. */
    a_a->num_nodes += a_b->num_nodes;
}

static void
bhp_p_union(cw_bhp_t *a_a, cw_bhp_t *a_b)
{
    cw_bhpi_t *prev_node, *curr_node, *next_node;

    cw_assert(a_a->priority_compare == a_b->priority_compare);

    bhp_p_merge(a_a, a_b);

    if (a_a->head != NULL)
    {
	prev_node = NULL;
	curr_node = a_a->head;
	next_node = curr_node->sibling;
	while (next_node != NULL)
	{
	    if ((curr_node->degree != next_node->degree)
		|| ((next_node->sibling != NULL) &&
		    (next_node->sibling->degree == curr_node->degree)))
	    {
		/* Either these two roots are unequal, or we're looking at the
		 * first two of three roots of equal degree (can happen because
		 * of merge (2) plus ripple carry (1)). */
		prev_node = curr_node;
		curr_node = next_node;
	    }
	    else if (a_a->priority_compare(curr_node->priority,
					   next_node->priority) != 1)
	    {
		/* The priority of the root of curr_node is <= the priority of
		 * the root of next_node. */
		curr_node->sibling = next_node->sibling;
		bhp_p_bin_link(curr_node, next_node);
	    }
	    else
	    {
		/* The priority of the root of curr_node is > the priority of
		 * the root of next_node. */
		if (prev_node == NULL)
		{
		    a_a->head = next_node;
		}
		else
		{
		    prev_node->sibling = next_node;
		}
		bhp_p_bin_link(next_node, curr_node);
		curr_node = curr_node->parent;
	    }
	    next_node = curr_node->sibling;
	}
    }
}

cw_bhp_t *
bhp_new(cw_bhp_t *a_bhp, bhp_prio_comp_t *a_prio_comp)
{
    cw_bhp_t *retval;

    cw_check_ptr(a_prio_comp);

    if (a_bhp == NULL)
    {
	retval = (cw_bhp_t *) cw_malloc(sizeof(cw_bhp_t));
	retval->is_malloced = true;
    }
    else
    {
	retval = a_bhp;
	retval->is_malloced = false;
    }

    retval->head = NULL;
    retval->num_nodes = 0;
    retval->priority_compare = a_prio_comp;

#ifdef CW_DBG
    retval->magic = CW_BHP_MAGIC;
#endif

    return retval;
}

void
bhp_delete(cw_bhp_t *a_bhp)
{
    cw_check_ptr(a_bhp);
    cw_dassert(a_bhp->magic == CW_BHP_MAGIC);

    /* Empty the heap. */
    if (a_bhp->head != NULL)
    {
	while (a_bhp->num_nodes > 0)
	{
	    bhp_min_del(a_bhp, NULL, NULL, NULL);
	}
    }
    if (a_bhp->is_malloced)
    {
	cw_free(a_bhp);
    }
#ifdef CW_DBG
    else
    {
	memset(a_bhp, 0x5a, sizeof(cw_bhp_t));
    }
#endif
}

void
bhp_dump(cw_bhp_t *a_bhp)
{
    cw_check_ptr(a_bhp);
    cw_dassert(a_bhp->magic == CW_BHP_MAGIC);

    fprintf(stderr, "=== bhp_dump() start ==============================\n");
    
    fprintf(stderr, "num_nodes: %llu\n", a_bhp->num_nodes);
    if
    (a_bhp->head != NULL)
    {
	bhp_p_dump(a_bhp->head, 0, NULL);
    }
    fprintf(stderr, "=== bhp_dump() end ================================\n");
}

void
bhp_insert(cw_bhp_t *a_bhp, const void *a_priority, const void *a_data,
	   cw_bhpi_t *a_bhpi)
{
    cw_bhp_t temp_heap;

    cw_check_ptr(a_bhp);
    cw_dassert(a_bhp->magic == CW_BHP_MAGIC);

    /* Initialize a_bhpi. */
    if (a_bhpi != NULL)
    {
	a_bhpi->is_malloced = false;
    }
    else
    {
	a_bhpi = (cw_bhpi_t *) cw_malloc(sizeof(cw_bhpi_t));
	a_bhpi->is_malloced = true;
    }
    a_bhpi->priority = a_priority;
    a_bhpi->data = a_data;
    a_bhpi->parent = NULL;
    a_bhpi->child = NULL;
    a_bhpi->sibling = NULL;
    a_bhpi->degree = 0;

#ifdef CW_DBG
    a_bhpi->magic = CW_BHPI_MAGIC;
#endif
    
    /* Create and initialize temp_heap. */
    bhp_new(&temp_heap, a_bhp->priority_compare);
    temp_heap.head = a_bhpi;
    temp_heap.num_nodes = 1;

    /* Combine this heap and temp_heap. */
    bhp_p_union(a_bhp, &temp_heap);

    /* Destroy the old heap. */
    temp_heap.head = NULL;
    bhp_delete(&temp_heap);
}

bool
bhp_min_find(cw_bhp_t *a_bhp, void **r_priority, void **r_data)
{
    bool retval;
    cw_bhpi_t *curr_min, *curr_pos;

    cw_check_ptr(a_bhp);
    cw_dassert(a_bhp->magic == CW_BHP_MAGIC);

    if (a_bhp->head != NULL)
    {
	retval = false;

	curr_min = a_bhp->head;
	curr_pos = a_bhp->head->sibling;

	while (curr_pos != NULL)
	{
	    /* Check if curr_pos is less than curr_min For priority_compare(a,
	     * b), -1 means a < b, 0 means a == b, 1 means a > b. */
	    if (a_bhp->priority_compare(curr_pos->priority,
					curr_min->priority) == -1)
	    {
		curr_min = curr_pos;
	    }
	    curr_pos = curr_pos->sibling;
	}

	/* We've found a minimum priority item now, so point *r_priority and
	 * *r_data to it. */
	if (r_priority != NULL)
	{
	    *r_priority = (void *) curr_min->priority;
	}
	if (r_data != NULL)
	{
	    *r_data = (void *) curr_min->data;
	}
    }
    else
    {
	retval = true;
    }

    return retval;
}

bool
bhp_min_del(cw_bhp_t *a_bhp, void **r_priority, void **r_data,
	    cw_bhpi_t **r_bhpi)
{
    bool retval;
    cw_bhpi_t *prev_pos, *curr_pos, *next_pos, *before_min, *curr_min;
    cw_bhp_t temp_heap;

    cw_check_ptr(a_bhp);
    cw_dassert(a_bhp->magic == CW_BHP_MAGIC);

    if (a_bhp->num_nodes == 0)
    {
	retval = true;
    }
    else
    {
	retval = false;

	/* Find a root with minimum priority. */
	before_min = NULL;
	prev_pos = NULL;
	curr_pos = a_bhp->head;
	curr_min = curr_pos;
	while (curr_pos != NULL)
	{
	    if (a_bhp->priority_compare(curr_pos->priority,
					curr_min->priority) == -1)
	    {
		/* Found a new minimum. */
		curr_min = curr_pos;
		before_min = prev_pos;
	    }
	    prev_pos = curr_pos;
	    curr_pos = curr_pos->sibling;
	}

	/* Take the minimum root out of the list. */
	if (before_min == NULL)
	{
	    /* Minimum root is the first in the list, so move the head pointer
	     * forward. */
	    a_bhp->head = curr_min->sibling;
	}
	else
	{
	    /* Attach previous and next roots together. */
	    before_min->sibling = curr_min->sibling;
	}

	/* Reverse order of curr_min's children. */
	prev_pos = NULL;
	curr_pos = curr_min->child;
	if (curr_pos != NULL)
	{
	    next_pos = curr_pos->sibling;
	}

	while (curr_pos != NULL)
	{
	    curr_pos->parent = NULL;
	    curr_pos->sibling = prev_pos;

	    prev_pos = curr_pos;
	    curr_pos = next_pos;
	    if (next_pos != NULL)
	    {
		next_pos = next_pos->sibling;
	    }
	}

	/* Create a temporary heap and initialize it. */
	bhp_new(&temp_heap, a_bhp->priority_compare);
	temp_heap.head = prev_pos;
	bhp_p_union(a_bhp, &temp_heap);

	/* Destroy the old heap. */
	temp_heap.head = NULL;
	bhp_delete(&temp_heap);

	a_bhp->num_nodes--;

	/* Now point *r_priority and *r_data to the item and free the space
	 * taken up by the item structure. */
	if (r_priority != NULL)
	{
	    *r_priority = (void *) curr_min->priority;
	}
	if (r_data != NULL)
	{
	    *r_data = (void *) curr_min->data;
	}

	/* Clean up associated bhpi. */
	if (curr_min->is_malloced)
	{
	    cw_free(curr_min);
	}
	else
	{
#ifdef CW_DBG
	    memset(curr_min, 0x5a, sizeof(cw_bhpi_t));
#endif
	    if (r_bhpi != NULL)
	    {
		*r_bhpi = curr_min;
	    }
	}
    }

    return retval;
}

uint64_t
bhp_size_get(cw_bhp_t *a_bhp)
{
    uint64_t retval;

    cw_check_ptr(a_bhp);
    cw_dassert(a_bhp->magic == CW_BHP_MAGIC);

    retval = a_bhp->num_nodes;

    return retval;
}

void
bhp_union(cw_bhp_t *a_a, cw_bhp_t *a_b)
{
    cw_check_ptr(a_a);
    cw_dassert(a_a->magic == CW_BHP_MAGIC);
    cw_check_ptr(a_b);
    cw_dassert(a_b->magic == CW_BHP_MAGIC);

    bhp_p_union(a_a, a_b);

    /* Destroy the old heap. */
    a_b->head = NULL;
    bhp_delete(a_b);
}

int32_t
bhp_uint32_priority_compare(const void *a_a, const void *a_b)
{
    int32_t retval;
    uint32_t a = *(uint32_t *) a_a;
    uint32_t b = *(uint32_t *) a_b;

    cw_check_ptr(a_a);
    cw_check_ptr(a_b);

    if (a < b)
    {
	retval = -1;
    }
    else if (a > b)
    {
	retval = 1;
    }
    else
    {
	retval = 0;
    }

    return retval;
}

int32_t
bhp_int32_priority_compare(const void *a_a, const void *a_b)
{
    int32_t retval;
    int32_t a = *(int32_t *) a_a;
    int32_t b = *(int32_t *) a_b;

    cw_check_ptr(a_a);
    cw_check_ptr(a_b);

    if (a < b)
    {
	retval = -1;
    }
    else if (a > b)
    {
	retval = 1;
    }
    else
    {
	retval = 0;
    }

    return retval;
}

int32_t
bhp_uint64_priority_compare(const void *a_a, const void *a_b)
{
    int32_t retval;
    uint64_t a = *(uint64_t *) a_a;
    uint64_t b = *(uint64_t *) a_b;

    cw_check_ptr(a_a);
    cw_check_ptr(a_b);

    if (a < b)
    {
	retval = -1;
    }
    else if (a > b)
    {
	retval = 1;
    }
    else
    {
	retval = 0;
    }

    return retval;
}
