/* -*- mode: c ; c-file-style: "canonware-c-style" -*-
 ******************************************************************************
 *
 * <Copyright = jasone>
 * <License>
 *
 ******************************************************************************
 *
 * Version: Crux <Version = crux>
 *
 ******************************************************************************
 *
 * cw_tr_t encapsulates an (optionally) unrooted bifurcating tree, and nodes are
 * represented by cw_trn_t.
 *
 ******************************************************************************
 *
 * cw_trn_t contains a time-efficient representation of a node of an
 * (optionally) unrooted bifurcating phlyogenetic tree.
 *
 ******************************************************************************/

#include "../include/modcrux.h"

#ifdef CW_DBG
/* Validate an individual trn. */
static cw_bool_t
trn_p_validate(cw_trn_t *a_trn)
{
    cw_check_ptr(a_trn);
    cw_assert(a_trn->magic_a == CW_TRN_MAGIC_A);
    cw_assert(a_trn->magic_b == CW_TRN_MAGIC_B);

    /* A trn must either be a leaf node, an internal node, or a root node.
     *
     * This check allows an internal node to have less than 3 neighbors, since
     * multiple calls are necessary to set up all the neighbor pointers.  This
     * also allows root nodes, which have two neighbors.
     *
     * A leaf node is not required to have a neighbor. */
    if (a_trn->taxon_num != CW_TRN_TAXON_NONE)
    {
	cw_uint32_t i, j, nneighbors, nloops;

	for (i = nneighbors = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    if (a_trn->neighbors[i] != NULL)
	    {
		nneighbors++;
		for (j = nloops = 0; j < CW_TRN_MAX_NEIGHBORS; j++)
		{
		    if (a_trn->neighbors[i]->neighbors[j] == a_trn)
		    {
			nloops++;
		    }
		}
		cw_assert(nloops == 1);
	    }
	}
	cw_assert(nneighbors <= 1);
    }

    return TRUE;
}
#endif

static void
trn_p_delete_recurse(cw_trn_t *a_trn)
{
    cw_uint32_t i;
    cw_trn_t *trn;

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL)
	{
	    trn = a_trn->neighbors[i];
	    trn_detach(a_trn, trn);
	    trn_p_delete_recurse(trn);
	}
    }

    trn_delete(a_trn);
}

#ifdef CW_DBG
/* Return the number of taxa with number a_taxon_num in the subtree rooted at
 * a_trn. */
static cw_uint32_t
trn_p_validate_recurse(cw_trn_t *a_trn, cw_trn_t *a_prev,
		       cw_uint32_t a_taxon_num)
{
    cw_uint32_t retval;
    cw_uint32_t i;

    trn_p_validate(a_trn);

    if (a_trn->taxon_num != CW_TRN_TAXON_NONE)
    {
	/* Leaf node. */
	cw_assert(a_trn->neighbors[0] != NULL
		  || (a_trn->neighbors[1] == NULL
		      && a_trn->neighbors[1] == NULL));
	for (i = 1; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    cw_assert(a_trn->neighbors[i] == NULL);
	}

	if (a_trn->taxon_num == a_taxon_num)
	{
	    retval = 1;
	}
	else
	{
	    retval = 0;
	}
    }
    else
    {
	/* Internal node. */
	for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    cw_assert(a_trn->neighbors[i] != NULL);
	}

	retval = 0;
    }

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    retval += trn_p_validate_recurse(a_trn->neighbors[i],
					     a_trn, a_taxon_num);
	}
    }

    return retval;
}
#endif

/* Recursively traverse the tree and find the lowest numbered taxon. */
static cw_trn_t *
trn_p_root_get(cw_trn_t *a_trn, cw_trn_t *a_prev, cw_trn_t *a_root)
{
    cw_trn_t *retval, *root, *troot;
    cw_uint32_t i;

    if (a_trn->taxon_num != CW_TRN_TAXON_NONE
	&& (a_root == NULL || a_trn->taxon_num < a_root->taxon_num))
    {
	retval = a_trn;
	root = a_trn;
    }
    else
    {
	retval = NULL;
	root = a_root;
    }

    /* Iterate over neighbors. */
    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    troot = trn_p_root_get(a_trn->neighbors[i], a_trn, root);
	    if (troot != NULL)
	    {
		retval = troot;
		root = troot;
	    }
	}
    }

    return retval;
}

/* Recursively traverse the tree and count the number of taxa. */
static cw_uint32_t
trn_p_ntaxa_get_recurse(cw_trn_t *a_trn, cw_trn_t *a_prev)
{
    cw_uint32_t retval;
    cw_uint32_t i;

    cw_dassert(trn_p_validate(a_trn));

    if (a_trn->taxon_num != CW_TRN_TAXON_NONE)
    {
	/* Leaf node. */
	retval = 1;
    }
    else
    {
	/* Internal node. */
	retval = 0;
    }

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    retval += trn_p_ntaxa_get_recurse(a_trn->neighbors[i], a_trn);
	}
    }

    return retval;
}

#ifdef CW_DBG
/* trn. */
static cw_uint32_t
trn_p_nedges_get(cw_trn_t *a_trn)
{
    cw_uint32_t retval, ntaxa;

    cw_dassert(trn_p_validate(a_trn));

    ntaxa = trn_p_ntaxa_get_recurse(a_trn, NULL);
    if (ntaxa > 1)
    {
	retval = (ntaxa * 2) - 3;
    }
    else
    {
	retval = 0;
    }

    return retval;
}
#endif

static cw_bool_t
trn_p_edge_get_recurse(cw_trn_t *a_trn, cw_uint32_t a_edge,
		       cw_trn_t *a_prev, cw_uint32_t *r_edge_count,
		       cw_trn_t **r_trn, cw_uint32_t *r_neighbor)
{
    cw_bool_t retval;
    cw_uint32_t i;

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    /* Increment edge count before recursing.  If the edge count has
	     * reached the desired value, return this trn and neighbor index,
	     * and terminate recursion. */
	    (*r_edge_count)++;
	    if (*r_edge_count > a_edge)
	    {
		cw_assert(*r_edge_count == a_edge + 1);
		*r_trn = a_trn;
		*r_neighbor = i;

		retval = TRUE;
		goto RETURN;
	    }

	    /* Recurse into neighbor subtrees. */
	    if (trn_p_edge_get_recurse(a_trn->neighbors[i], a_edge,
				       a_trn, r_edge_count,
				       r_trn, r_neighbor))
	    {
		retval = TRUE;
		goto RETURN;
	    }
	}
    }

    retval = FALSE;
    RETURN:
    return retval;
}

CW_P_INLINE void
trn_p_edge_get(cw_trn_t *a_trn, cw_uint32_t a_edge, cw_trn_t **r_trn,
	       cw_uint32_t *r_neighbor)
{
    cw_uint32_t edge_count = 0;
#ifdef CW_DBG
    cw_bool_t found =
#endif
	trn_p_edge_get_recurse(a_trn, a_edge, NULL, &edge_count, r_trn,
			       r_neighbor);
    cw_dassert(found);
}

static void
trn_p_bisection_edge_get_recurse(cw_trn_t *a_trn, cw_trn_t *a_other,
				 cw_trn_t *a_prev, cw_uint32_t *r_edge_count,
				 cw_uint32_t *r_bisection_edge)
{
    cw_uint32_t i, prev_edge_count;

    /* Save the previous edge count, in case it ends up being the index of the
     * edge adjacent to the bisection. */
    prev_edge_count = *r_edge_count;

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL
	    && a_trn->neighbors[i] != a_prev
	    && a_trn->neighbors[i] != a_other)
	{
	    /* Increment edge count before recursing. */
	    (*r_edge_count)++;

	    /* Recurse into neighbor subtrees. */
	    trn_p_bisection_edge_get_recurse(a_trn->neighbors[i], a_other,
					     a_trn, r_edge_count,
					     r_bisection_edge);
	}
	else if (a_trn->neighbors[i] == a_other)
	{
	    /* Store the index of the edge adjacent to the bisection. */
	    if (prev_edge_count > 0)
	    {
		*r_bisection_edge = prev_edge_count - 1;
	    }
	    else
	    {
		*r_bisection_edge = CW_TRN_EDGE_NONE;
	    }
	}
    }
}

static void
trn_p_bisection_edges_get(cw_trn_t *a_trn, cw_uint32_t a_edge, cw_trt_t *r_trt)
{
    cw_uint32_t neighbor;
    cw_trn_t *trn, *adj_a, *adj_b;

    /* Count the number of edges that would be in each half of the tree, were it
     * bisected.  Also, determine which edges of the subtrees would be used to
     * reverse the bisection. */

    /* Get the trn's adjacent to the bisection edge. */
    trn_p_edge_get(a_trn, a_edge, &adj_a, &neighbor);
    adj_b = adj_a->neighbors[neighbor];

    /* Get the number of edges in the first half of the bisection, as well as
     * the index of the edge adjacent to the bisection. */
    r_trt->nedges_a = 0;
    trn_p_bisection_edge_get_recurse(a_trn, adj_b, NULL, &r_trt->nedges_a,
				     &r_trt->self_a);
    if (r_trt->nedges_a > 0)
    {
	r_trt->nedges_a--;
    }

    /* Get the lowest numbered taxon in the second half of the bisection. */
    trn = trn_p_root_get(adj_b, adj_a, NULL);

    /* Get the number of edges in the second half of the bisection, as well as
     * the index of the edge adjacent to the bisection. */
    r_trt->nedges_b = 0;
    trn_p_bisection_edge_get_recurse(trn, adj_a, NULL, &r_trt->nedges_b,
				     &r_trt->self_b);
    if (r_trt->nedges_b > 0)
    {
	r_trt->nedges_b--;
    }
}

static cw_bool_t
trn_p_edge_index_get_recurse(cw_trn_t *a_trn, cw_trn_t *a_prev,
			     cw_trn_t *a_trn_a, cw_trn_t *a_trn_b,
			     cw_uint32_t *r_edge_count)
{
    cw_bool_t retval;
    cw_uint32_t i;

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    /* If this is the desired edge, terminate recursion.  Increment edge
	     * count before recursing. */
	    if ((a_trn == a_trn_a && a_trn->neighbors[i] == a_trn_b)
		|| (a_trn == a_trn_b && a_trn->neighbors[i] == a_trn_a))
	    {
		retval = TRUE;
		goto RETURN;
	    }
	    (*r_edge_count)++;

	    /* Recurse into neighbor subtree. */
	    if (trn_p_edge_index_get_recurse(a_trn->neighbors[i], a_trn,
					     a_trn_a, a_trn_b,
					     r_edge_count))
	    {
		retval = TRUE;
		goto RETURN;
	    }
	}
    }

    retval = FALSE;
    RETURN:
    return retval;
}

static cw_uint32_t
trn_p_edge_index_get(cw_trn_t *a_trn, cw_trn_t *a_trn_a, cw_trn_t *a_trn_b)
{
    cw_uint32_t retval = 0;

    trn_p_edge_index_get_recurse(a_trn, NULL, a_trn_a, a_trn_b, &retval);

    return retval;
}

/* Convert a tree to canonical form by re-ordering the neighbors array such that
 * subtrees are in increasing order of minimum taxon number contained. */
static cw_uint32_t
trn_p_canonize(cw_trn_t *a_trn, cw_trn_t *a_prev)
{
    cw_uint32_t retval;
    cw_uint32_t i, j, t;
    cw_uint32_t subtree_mins[CW_TRN_MAX_NEIGHBORS - 1];
    cw_uint32_t subtree_inds[CW_TRN_MAX_NEIGHBORS - 1];
    cw_bool_t swapped;

    cw_assert(a_prev != NULL
	      || trn_p_root_get(a_trn, NULL, NULL) == a_trn);

    if (a_trn->taxon_num != CW_TRN_TAXON_NONE)
    {
	/* Leaf node. */
	retval = a_trn->taxon_num;
    }
    else
    {
	/* Internal node. */
	retval = CW_TRN_TAXON_NONE;
    }

    /* Iteratively canonize subtrees, keeping track of the minimum taxon
     * number seen overall, as well as for each subtree. */
    for (i = j = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    /* A neighboring subtree that hasn't been visited yet. */
	    cw_assert(j < (CW_TRN_MAX_NEIGHBORS - 1));
	    subtree_mins[j] = trn_p_canonize(a_trn->neighbors[i], a_trn);
	    if (subtree_mins[j] < retval)
	    {
		retval = subtree_mins[j];
	    }
	    subtree_inds[j] = i;
	    j++;
	}
    }

    /* Bubble sort the subtrees.  This algorithm works in the general case, and
     * in the case this code is actually designed for (bifurcating trees), it
     * only requires a couple of extra branches. */
    do
    {
	swapped = FALSE;

	for (i = 0; i + 1 < j; i++)
	{
	    if (subtree_mins[i] > subtree_mins[i + 1])
	    {
		swapped = TRUE;

		/* Swap subtrees. */
		trn_neighbors_swap(a_trn, subtree_inds[i], subtree_inds[i + 1]);

		/* Swap subtree_* arrays. */
		t = subtree_mins[i];
		subtree_mins[i] = subtree_mins[i + 1];
		subtree_mins[i + 1] = t;
	    
		t = subtree_inds[i];
		subtree_inds[i] = subtree_inds[i + 1];
		subtree_inds[i + 1] = t;
	    }
	}
    } while (swapped);

    return retval;
}

CW_P_INLINE void
trn_p_bisection_patch(cw_trn_t *a_trn, cw_trn_t **r_trn, cw_uint32_t *r_edge,
		      cw_trn_t **r_spare)
{
    /* Patch up the node adjacent to bisection.  There are two cases possible.
     * It is either a leaf node or an internal node.  For a leaf node, do
     * nothing.  For an internal node, join its neighbors together, and return
     * the node as a spare. */
    if (a_trn->taxon_num == CW_TRN_TAXON_NONE)
    {
	cw_trn_t *a, *b;
	cw_uint32_t i;

	/* Get a_trn's neighbors. */
	for (i = 0, a = b = NULL; b == NULL; i++)
	{
	    cw_assert(i < CW_TRN_MAX_NEIGHBORS);

	    if (a_trn->neighbors[i] != NULL)
	    {
		if (a == NULL)
		{
		    a = a_trn->neighbors[i];
		}
		else
		{
		    b = a_trn->neighbors[i];
		}
	    }
	}

	/* Detach. */
	trn_detach(a_trn, a);
	trn_detach(a_trn, b);

	/* Join. */
	trn_join(a, b);

	*r_trn = trn_p_root_get(a, NULL, NULL);
	trn_p_canonize(*r_trn, NULL);
	*r_edge = trn_p_edge_index_get(*r_trn, a, b);
	*r_spare = a_trn;
    }
    else
    {
	*r_trn = a_trn;
	*r_edge = CW_TRN_EDGE_NONE;
	*r_spare = NULL;
    }
}

CW_P_INLINE void
trn_p_bisect(cw_trn_t *a_trn, cw_uint32_t a_edge, cw_trn_t **r_trn_a,
	     cw_uint32_t *r_edge_a, cw_trn_t **r_trn_b,
	     cw_uint32_t *r_edge_b, cw_trn_t **r_spare_a,
	     cw_trn_t **r_spare_b)
{
    cw_trn_t *trn_a, *trn_b;
    cw_uint32_t edge;

    cw_assert(trn_p_root_get(a_trn, NULL, NULL) == a_trn);
    cw_assert(a_edge < trn_p_nedges_get(a_trn));

    /* Get the nodes to either side of the edge where the bisection will be
     * done. */
    trn_p_edge_get(a_trn, a_edge, &trn_a, &edge);
    trn_b = trn_a->neighbors[edge];

    /* Detach the two nodes. */
    trn_detach(trn_a, trn_b);

    /* Patch up the nodes adjacent to bisection. */
    trn_p_bisection_patch(trn_a, r_trn_a, r_edge_a, r_spare_a);
    trn_p_bisection_patch(trn_b, r_trn_b, r_edge_b, r_spare_b);

    /* Move *r_spare_b to *r_spare_a if *r_spare_a is NULL. */
    if (*r_spare_a == NULL)
    {
	*r_spare_a = *r_spare_b;
	*r_spare_b = NULL;
    }
}

CW_INLINE cw_trn_t *
trn_p_connection_patch(cw_trn_t *a_trn, cw_uint32_t a_edge,
		       cw_trn_t **ar_spare)
{
    cw_trn_t *retval;

    /* There are two cases possible.  a_trn is either a single leaf node, or a
     * larger subtree.  For a leaf node, do nothing.  For a larger subtree,
     * splice in the spare. */
    if (a_edge != CW_TRN_EDGE_NONE)
    {
	cw_trn_t *a, *b;
	cw_uint32_t edge;

	cw_check_ptr(*ar_spare);
	retval = *ar_spare;
	*ar_spare = NULL;

	trn_p_edge_get(a_trn, a_edge, &a, &edge);
	b = a->neighbors[edge];

	/* Detach a and b. */
	trn_detach(a, b);

	/* Join a and b to retval. */
	trn_join(retval, a);
	trn_join(retval, b);
    }
    else
    {
	retval = a_trn;
    }

    return retval;
}

CW_P_INLINE void
trn_p_connect(cw_trn_t *a_trn_a, cw_uint32_t a_edge_a,
	      cw_trn_t *a_trn_b, cw_uint32_t a_edge_b,
	      cw_trn_t **ar_spare_a, cw_trn_t **ar_spare_b,
	      cw_trn_t **r_trn, cw_uint32_t *r_edge)
{
    cw_trn_t *trn_a, *trn_b;

    cw_assert(trn_p_root_get(a_trn_a, NULL, NULL) == a_trn_a);
    cw_dassert((a_edge_a == CW_TRN_EDGE_NONE
		&& trn_p_ntaxa_get_recurse(a_trn_a, NULL) == 1)
	       || (a_edge_a < trn_p_nedges_get(a_trn_a)));
    cw_assert(trn_p_root_get(a_trn_b, NULL, NULL) == a_trn_b);
    cw_dassert((a_edge_b == CW_TRN_EDGE_NONE
		&& trn_p_ntaxa_get_recurse(a_trn_b, NULL) == 1)
	       || (a_edge_b < trn_p_nedges_get(a_trn_b)));

    /* Patch in internal nodes adjacent to where the connection will occur,
     * if necessary. */
    trn_a = trn_p_connection_patch(a_trn_a, a_edge_a, ar_spare_a);
    trn_b = trn_p_connection_patch(a_trn_b, a_edge_b,
				   ((*ar_spare_a != NULL)
				    ? ar_spare_a
				    : ar_spare_b));

    /* Join trn_a and trn_b. */
    trn_join(trn_a, trn_b);

    /* Get the tree root, and canonize the tree. */
    if (a_trn_a->taxon_num < a_trn_b->taxon_num)
    {
	*r_trn = a_trn_a;
    }
    else
    {
	*r_trn = a_trn_b;
    }
    trn_p_canonize(*r_trn, NULL);

    /* Get the edge index for the connection. */
    *r_edge = trn_p_edge_index_get(*r_trn, trn_a, trn_b);
}

void
trn_new(cw_trn_t *a_trn)
{
    cw_check_ptr(a_trn);

    memset(a_trn, 0, sizeof(cw_trn_t));
    a_trn->taxon_num = CW_TRN_TAXON_NONE;

#ifdef CW_DBG
    a_trn->magic_a = CW_TRN_MAGIC_A;
    a_trn->magic_b = CW_TRN_MAGIC_B;
#endif
}

void
trn_delete(cw_trn_t *a_trn)
{
#ifdef CW_DBG
    memset(a_trn, 0x5a, sizeof(cw_trn_t));
#endif
}

cw_uint32_t
trn_taxon_num_get(cw_trn_t *a_trn)
{
    cw_dassert(trn_p_validate(a_trn));

    return a_trn->taxon_num;
}

void
trn_taxon_num_set(cw_trn_t *a_trn, cw_uint32_t a_taxon_num)
{
    cw_dassert(trn_p_validate(a_trn));

    a_trn->taxon_num = a_taxon_num;
}

cw_trn_t *
trn_neighbor_get(cw_trn_t *a_trn, cw_uint32_t a_i)
{
    cw_dassert(trn_p_validate(a_trn));
    cw_assert(a_i < CW_TRN_MAX_NEIGHBORS);

    return a_trn->neighbors[a_i];
}

void
trn_neighbors_swap(cw_trn_t *a_trn, cw_uint32_t a_i, cw_uint32_t a_j)
{
    cw_trn_t *t_trn;

    cw_dassert(trn_p_validate(a_trn));
    cw_assert(a_i < CW_TRN_MAX_NEIGHBORS);
    cw_assert(a_j < CW_TRN_MAX_NEIGHBORS);
    cw_assert(a_i != a_j);

    t_trn = a_trn->neighbors[a_i];
    a_trn->neighbors[a_i] = a_trn->neighbors[a_j];
    a_trn->neighbors[a_j] = t_trn;
}

void
trn_join(cw_trn_t *a_a, cw_trn_t *a_b)
{
    cw_uint32_t i, j;

    cw_dassert(trn_p_validate(a_a));
    cw_dassert(trn_p_validate(a_b));
    cw_assert(a_a != a_b);
#ifdef CW_DBG
    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	cw_assert(a_a->neighbors[i] != a_b);
	cw_assert(a_b->neighbors[i] != a_a);
    }
#endif

    /* Find an empty slot in a_a. */
    for (i = 0; a_a->neighbors[i] != NULL; i++)
    {
	cw_assert(i < CW_TRN_MAX_NEIGHBORS);
    }

    /* Find an empty slot in a_b. */
    for (j = 0; a_b->neighbors[j] != NULL; j++)
    {
	cw_assert(j < CW_TRN_MAX_NEIGHBORS);
    }

    /* Join the two nodes. */
    a_a->neighbors[i] = a_b;
    a_b->neighbors[j] = a_a;

    cw_dassert(trn_p_validate(a_a));
    cw_dassert(trn_p_validate(a_b));
}

void
trn_detach(cw_trn_t *a_a, cw_trn_t *a_b)
{
    cw_uint32_t i, j;

    cw_dassert(trn_p_validate(a_a));
    cw_dassert(trn_p_validate(a_b));

    /* Find the slot in a_a that points to a_b. */
    for (i = 0; a_a->neighbors[i] != a_b; i++)
    {
	cw_assert(i < CW_TRN_MAX_NEIGHBORS);
    }

    /* Find the slot in a_b that points to a_a. */
    for (j = 0; a_b->neighbors[j] != a_a; j++)
    {
	cw_assert(j < CW_TRN_MAX_NEIGHBORS);
    }

    /* Detach the two nodes. */
    a_a->neighbors[i] = NULL;
    a_b->neighbors[j] = NULL;

    cw_dassert(trn_p_validate(a_a));
    cw_dassert(trn_p_validate(a_b));
}

void *
trn_aux_get(cw_trn_t *a_trn)
{
    cw_dassert(trn_p_validate(a_trn));

    return a_trn->aux;
}

void
trn_aux_set(cw_trn_t *a_trn, void *a_aux)
{
    cw_dassert(trn_p_validate(a_trn));

    a_trn->aux = a_aux;
}

/* tr. */
#ifdef CW_DBG
/* Validate a tree constructed with trn's. */
static cw_bool_t
tr_p_validate(cw_tr_t *a_tr, cw_bool_t a_validate_contiguous)
{
    cw_uint32_t i, j, n;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
    cw_assert(trn_p_root_get(a_tr->croot, NULL, NULL) == a_tr->croot);
    cw_assert(a_tr->ntaxa == trn_p_ntaxa_get_recurse(a_tr->croot, NULL));
    cw_assert((a_tr->rooted == FALSE
	       && a_tr->nedges == trn_p_nedges_get(a_tr->croot) - 1)
	      || a_tr->nedges == trn_p_nedges_get(a_tr->croot));

    /* Traverse the tree, and make sure that the following invariants hold:
     *
     * + Leaf nodes have a taxon number and precisely 1 neighbor.
     *
     * + Internal nodes have no taxon number and precisely 3 neighbors.
     *
     * + The root node (if there is one) has precisely 2 neighbors and no taxon
     *   number.
     *
     * + Each taxon number appears no more than once in the tree.
     *
     * These invariants allow gaps in the taxon numbering, which has the
     * potential to cause problems (not in the tr implementation but elswhere),
     * but requiring contiguous taxon numbering would make validating tree
     * bisections impossible.  Therefore, contiguous taxon numbering is only
     * optionally validated. */
    if (a_tr->rooted)
    {
	cw_trn_t *left, *right;

	/* Rooted tree. */

	/* Get left and right subtrees. */
	for (i = 0, left = right = NULL; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    if (a_tr->root->neighbors[i] != NULL)
	    {
		if (left == NULL)
		{
		    left = a_tr->root->neighbors[i];
		}
		else
		{
		    right = a_tr->root->neighbors[i];
		    break;
		}
	    }
	}

	for (i = j = 0; j < a_tr->ntaxa; i++, j += n)
	{
	    n = trn_p_validate_recurse(left, a_tr->root, i);
	    n += trn_p_validate_recurse(right, a_tr->root, i);
	    cw_assert(n <= 1);
	}
    }
    else
    {
	/* Unrooted tree. */
	for (i = j = 0; j < a_tr->ntaxa; i++, j += n)
	{
	    n = trn_p_validate_recurse(a_tr->croot, NULL, i);
	    cw_assert(n <= 1);
	}
    }

    if (a_validate_contiguous)
    {
	cw_assert(i == a_tr->ntaxa);
    }

    return TRUE;
}
#endif

static void
tr_p_tbr(cw_tr_t *a_tr, cw_uint32_t a_bisect, cw_uint32_t a_reconnect_a,
	 cw_uint32_t a_reconnect_b, cw_uint32_t *r_bisect,
	 cw_uint32_t *r_reconnect_a, cw_uint32_t *r_reconnect_b)
{
    cw_trn_t *trn, *trn_a, *trn_b, *spare_a, *spare_b;

    /* Bisect. */
    trn_p_bisect(a_tr->croot, a_bisect,
		 &trn_a, r_reconnect_a,
		 &trn_b, r_reconnect_b,
		 &spare_a, &spare_b);

    /* Reconnect. */
    trn_p_connect(trn_a, a_reconnect_a,
		  trn_b, a_reconnect_b,
		  &spare_a, &spare_b,
		  &trn, r_bisect);
    cw_assert(trn == a_tr->croot);
    cw_assert(spare_a == NULL);
    cw_assert(spare_b == NULL);

    cw_dassert(tr_p_validate(a_tr, TRUE));
}

static int
tr_p_trt_compare(const void *a_key, const void *a_val)
{
    int retval;
    const cw_trt_t *key = (const cw_trt_t *) a_key;
    const cw_trt_t *val = (const cw_trt_t *) a_val;

    if (key->offset < val->offset)
    {
	retval = -1;
    }
    else if (key->offset < (&val[1])->offset)
    {
	retval = 0;
    }
    else
    {
	retval = 1;
    }

    return retval;
}

static void
tr_p_tbr_init(cw_tr_t *a_tr, cw_bool_t a_use_ri)
{
    cw_uint32_t i, j, n, offset, a, b;

    if (a_tr->trt == NULL)
    {
	/* Build table to make neighbor iteration fast. */
	if (a_tr->nedges > 0)
	{
	    a_tr->trt = (cw_trt_t *) nxa_malloc(sizeof(cw_trt_t)
						* (a_tr->nedges + 1));
	}
    }

    for (i = j = offset = 0; i < a_tr->nedges; i++)
    {
	/* Record offset. */
	a_tr->trt[j].offset = offset;

	/* Record bisection edge. */
	a_tr->trt[j].bisect_edge = i;

	/* Record number of subtree edges. */
	a_tr->trt[j].nedges_a = 0;
	a_tr->trt[j].nedges_b = 0;

	/* Set trt[i].{nedges,self}_[ab]. */
	trn_p_bisection_edges_get(a_tr->croot, i, &a_tr->trt[j]);

	/* Update offset. */
	if (a_tr->trt[j].nedges_a != 0)
	{
	    a = a_tr->trt[j].nedges_a;
	}
	else
	{
	    a = 1;
	}

	if (a_tr->trt[j].nedges_b != 0)
	{
	    b = a_tr->trt[j].nedges_b;
	}
	else
	{
	    b = 1;
	}

	n = (a * b) - 1;
	if (n != 0)
	{
	    offset += n;
	    j++;
	}
    }
    a_tr->trt[j].offset = offset;

    /* It may be that not all bisections result in neighbors, so the table may
     * not be full, so keep track of the number of valid elements (not counting
     * the trailing one that stores the total number of TBR neighbors). */
    a_tr->trtused = j;
}

void
tr_new(cw_tr_t *a_tr, cw_trn_t *a_trn)
{
    cw_check_ptr(a_tr);
    cw_dassert(trn_p_validate(a_trn));

    memset(a_tr, 0, sizeof(cw_tr_t));

#ifdef CW_DBG
    {
	cw_uint32_t i, nneighbors;

	/* Make sure this is a rooted tree. */
	for (i = nneighbors = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    if (a_trn->neighbors[i] != NULL)
	    {
		nneighbors++;
	    }
	}
	cw_assert(nneighbors == 2);
	cw_assert(a_trn->taxon_num == CW_TRN_TAXON_NONE);
    }
#endif
    a_tr->rooted = TRUE;
    a_tr->root = a_trn;

    a_tr->croot = trn_p_root_get(a_trn, NULL, NULL);

    a_tr->ntaxa = trn_p_ntaxa_get_recurse(a_trn, NULL);

    if (a_tr->ntaxa > 1)
    {
	a_tr->nedges = (a_tr->ntaxa * 2) - 3;
    }
    else
    {
	a_tr->nedges = 0;
    }

#ifdef CW_DBG
    a_tr->magic = CW_TR_MAGIC;
#endif

    trn_p_canonize(a_tr->croot, NULL);
}

void
tr_delete(cw_tr_t *a_tr, cw_bool_t a_delete_trns)
{
    if (a_delete_trns)
    {
	cw_dassert(tr_p_validate(a_tr, TRUE));
	if (a_tr->rooted == FALSE)
	{
	    trn_delete(a_tr->root);
	}

	trn_p_delete_recurse(a_tr->croot);
    }

    if (a_tr->trt != NULL)
    {
	nxa_free(a_tr->trt, sizeof(cw_trt_t) * (a_tr->nedges + 1));
    }

#ifdef CW_DBG
    memset(a_tr, 0x5a, sizeof(cw_tr_t));
#endif
}

cw_uint32_t
tr_ntaxa_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    return a_tr->ntaxa;
}

cw_uint32_t
tr_nedges_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    return a_tr->nedges;
}

/* Do an in-order traversal of the tree and return the node that neighbors the
 * edge being sought, along with which neighbor of that node the edge is
 * between. */
void
tr_edge_get(cw_tr_t *a_tr, cw_uint32_t a_edge, cw_trn_t **r_trn,
	    cw_uint32_t *r_neighbor)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    if (a_edge != CW_TRN_EDGE_NONE)
    {
	cw_assert(a_tr->rooted == FALSE);
	cw_assert(trn_p_ntaxa_get_recurse(a_tr->croot, NULL) > 1);

	trn_p_edge_get(a_tr->croot, a_edge, r_trn, r_neighbor);
    }
}

cw_uint32_t
tr_edge_index_get(cw_tr_t *a_tr, cw_trn_t *a_trn_a, cw_trn_t *a_trn_b)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    return trn_p_edge_index_get(a_tr->croot, a_trn_a, a_trn_b);
}

cw_trn_t *
tr_croot_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    return a_tr->croot;
}

cw_trn_t *
tr_root_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    return a_tr->root;
}

cw_bool_t
tr_rooted(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    return a_tr->rooted;
}

void
tr_root(cw_tr_t *a_tr, cw_trn_t *a_trn_a, cw_trn_t *a_trn_b)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));
    cw_assert(a_tr->rooted == FALSE);

    trn_detach(a_trn_a, a_trn_b);
    trn_join(a_tr->root, a_trn_a);
    trn_join(a_tr->root, a_trn_b);

    a_tr->rooted = TRUE;
}

void
tr_unroot(cw_tr_t *a_tr)
{
    cw_uint32_t i;
    cw_trn_t *trn_a, *trn_b;

    cw_dassert(tr_p_validate(a_tr, TRUE));
    cw_assert(a_tr->rooted);

    trn_a = NULL;
    trn_b = NULL;
    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_tr->root->neighbors[i] != NULL)
	{
	    if (trn_a == NULL)
	    {
		trn_a = a_tr->root->neighbors[i];
	    }
	    else
	    {
		trn_b = a_tr->root->neighbors[i];
		break;
	    }
	}
    }
    cw_check_ptr(trn_a);
    cw_check_ptr(trn_b);

    trn_detach(a_tr->root, trn_a);
    trn_detach(a_tr->root, trn_b);
    trn_join(trn_a, trn_b);

    a_tr->rooted = FALSE;
}

void
tr_canonize(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));
    cw_assert(trn_p_root_get(a_tr->croot, NULL, NULL) == a_tr->croot);

    trn_p_canonize(a_tr->croot, NULL);
}

void
tr_tbr(cw_tr_t *a_tr, cw_uint32_t a_bisect, cw_uint32_t a_reconnect_a,
       cw_uint32_t a_reconnect_b, cw_uint32_t *r_bisect,
       cw_uint32_t *r_reconnect_a, cw_uint32_t *r_reconnect_b)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));
    cw_assert(a_tr->rooted == FALSE);

    /* Perform TBR. */
    tr_p_tbr(a_tr, a_bisect, a_reconnect_a, a_reconnect_b, r_bisect,
	     r_reconnect_a, r_reconnect_b);
}

cw_uint32_t
tr_tbr_nneighbors_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    tr_p_tbr_init(a_tr, FALSE);

    return a_tr->trt[a_tr->trtused].offset;
}

void
tr_tbr_neighbor_get(cw_tr_t *a_tr, cw_uint32_t a_neighbor,
		    cw_uint32_t *r_bisect, cw_uint32_t *r_reconnect_a,
		    cw_uint32_t *r_reconnect_b)
{
    cw_trt_t key, *trt;
    cw_uint32_t rem, nedges_a, nedges_b, a, b;

    cw_dassert(tr_p_validate(a_tr, TRUE));

    tr_p_tbr_init(a_tr, FALSE);
    cw_assert(a_neighbor < a_tr->trt[a_tr->trtused].offset);

    /* Get the bisection edge. */
    key.offset = a_neighbor;
    trt = bsearch(&key, a_tr->trt, a_tr->trtused, sizeof(cw_trt_t),
		  tr_p_trt_compare);
    cw_check_ptr(trt);
    *r_bisect = trt->bisect_edge;

    /* Get the reconnection edges. */
    rem = a_neighbor - trt->offset;

    nedges_a = trt->nedges_a;
    nedges_b = trt->nedges_b;

    /* If the reconnection edges happen to be those that would reverse the
     * bisection, instead return the last possible reconnection combination for
     * this bisection.  This results in a rather strange ordering for the
     * enumeration, but is always correct. */

    /* No edges in one or both of the subtrees must be handled specially. */
    if (nedges_a == 0)
    {
	/* {0,b}. */

	/* A 2-taxon tree has no TBR neighbors. */
	cw_assert(nedges_b != 0);

	a = CW_TRN_EDGE_NONE;
	b = rem;

	if (a == trt->self_a && b == trt->self_b)
	{
	    b = nedges_b - 1;
	}
    }
    else if (nedges_b == 0)
    {
	/* {a,0}. */
	a = rem;
	b = CW_TRN_EDGE_NONE;

	if (a == trt->self_a && b == trt->self_b)
	{
	    a = nedges_a - 1;
	}
    }
    else
    {
	/* {a,b}. */
	a = rem / nedges_b;
	b = rem % nedges_b;

	/* (a * b) must be less than (nedges_a * nedges_b), since the last
	 * combination is reserved to replace the combination that corresponds
	 * to undoing the bisection. */
	cw_assert(a * b < nedges_a * nedges_b);

	if (a == trt->self_a && b == trt->self_b)
	{
	    a = nedges_a - 1;
	    b = nedges_b - 1;
	}
    }

    *r_reconnect_a = a;
    *r_reconnect_b = b;
}

void *
tr_aux_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    return a_tr->aux;
}

void
tr_aux_set(cw_tr_t *a_tr, void *a_aux)
{
    cw_dassert(tr_p_validate(a_tr, TRUE));

    a_tr->aux = a_aux;
}
