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
 ******************************************************************************/

#include "../include/modcrux.h"

typedef struct cw_tr_ps_s cw_tr_ps_t;
typedef struct cw_trn_s cw_trn_t;
typedef struct cw_tre_s cw_tre_t;
typedef struct cw_trt_s cw_trt_t;

/* Partial parsimony score information. */
struct cw_tr_ps_s
{
    /* Parent which most recently used this node's partial score when caching
     * its results.  Both children must still point to the parent in order for
     * the cached results to be valid. */
    cw_tr_node_t parent;

    /* Sum of the subtree scores, and this node's score, given particular
     * children.  In order for this to be useful, both childrens' parent
     * pointers must still point to this node. */
    cw_uint32_t subtrees_score;
    cw_uint32_t node_score;

    /* chars points to an array of Fitch parsimony state sets.  Each element in
     * the array contains a bitmap representation of a subset of {ACGT} in the 4
     * least significant bits.  T is the least significant bit.  1 means that a
     * nucleotide is in the set. */
    cw_uint32_t *chars;
    cw_uint32_t nchars;
};

/* Tree node for an unrooted bifurcating phylogenetic tree. */
struct cw_trn_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TRN_MAGIC 0x63329478
#endif

    /* Auxiliary opaque data pointer.  This is used by the treenode wrapper code
     * for reference iteration. */
    void *aux;

    /* If CW_TR_NODE_TAXON_NONE, then the node is not a leaf node. */
    cw_uint32_t taxon_num;

    /* Pointers to neighbors.  Only the first element is used if the node is a
     * leaf node. */
    cw_tr_node_t neighbors[CW_TR_NODE_MAX_NEIGHBORS];

    /* Edge indices for the edges that correspond to neighbors. */
    cw_uint32_t edges[CW_TR_NODE_MAX_NEIGHBORS];

    /* Used for Fitch parsimony scoring. */
    cw_tr_ps_t *ps;
};

/* Tree edge information. */
struct cw_tre_s
{
    /* Nodes adjacent to this edge. */
    cw_tr_node_t node_a;
    cw_tr_node_t node_b;

    /* Used for Fitch parsimony scoring. */
    cw_tr_ps_t *ps;
};

/* TBR neighbor. */
struct cw_trt_s
{
    /* Number of neighbors that can be reached by doing TBR at edges before this
     * one.  This is also the neighbor number of the first neighbor that can be
     * reached by doing TBR on this edge. */
    cw_uint32_t offset;

    /* Bisection edge. */
    cw_uint32_t bisect_edge;

    /* Root nodes (lowest-numbered taxa) in the two subtrees.  root_a is always
     * less than root_b. */
    cw_tr_node_t root_a;
    cw_tr_node_t root_b;

    /* Nodes adjacent to the bisection edge.  These correspond to the subtrees
     * rooted at root_[ab]. */
    cw_tr_node_t adj_a;
    cw_tr_node_t adj_b;

    /* Number of edges in the two subtrees.  Note that 0 and 1 are different
     * logical cases, but the number of connections possible for those two cases
     * is the same. */
    cw_uint32_t nedges_a;
    cw_uint32_t nedges_b;

    /* Edge indices for the edges that will reverse bisection.  This is used to
     * avoid enumerating reconnections that undo the bisections. */
    cw_uint32_t self_a;
    cw_uint32_t self_b;
};

struct cw_tr_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TR_MAGIC 0x39886394
#endif

    /* Used for memory allocation. */
    cw_mema_t *mema;

    /* Auxiliary opaque data pointer.  This is used by the treenode wrapper code
     * for reference iteration. */
    void *aux;

    /* TRUE if this tree has been modified since the internal state (ntaxa,
     * nedges, trt, tre) was updated, FALSE otherwise. */
    cw_bool_t modified;

    /* Base of the tree (may or may not be set). */
    cw_tr_node_t base;

    /* Number of taxa in tree. */
    cw_uint32_t ntaxa;

    /* Number of edges in tree.  This can be derived from ntaxa, but is used
     * often enough to make storing it worthwhile. */
    cw_uint32_t nedges;

    /* Array of information about edges.  There are always nedges elements in
     * tre -- nedges and tre are always kept in sync. */
    cw_tre_t *tres;

    /* Array of triplets that store per-edge information that is used for
     * TBR-related functions.  There is one more element in trt than there are
     * edges in the tree.  This is critical to the way binary searching on the
     * array is done, and it also makes it easy to get the total number of
     * TBR neighbors this tree has (trt[nedges].offset).
     *
     * Only the first trtused elements are valid, since not all bisection edges
     * necessarily result in neighbors. */
    cw_trt_t *trt;
    cw_uint32_t trtused;

    /* Pointer to an array of trn's.  ntrns is the total number of trn's, not
     * all of which are necessarily in use.  The first element is reserved as a
     * temporary, which gets used whereever a single temporary node is briefly
     * needed. */
    cw_trn_t *trns;
    cw_uint32_t ntrns;

    /* Index of first spare trn in the spares stack.  trns[spares].neighbors[0]
     * is used for list linkage. */
    cw_tr_node_t spares;
};

/******************************************************************************/

#ifdef CW_DBG
static cw_bool_t
tr_p_validate(cw_tr_t *a_tr);

/* Validate a node. */
static cw_bool_t
tr_p_node_validate(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_trn_t *trn;
    cw_uint32_t i, j, nneighbors, nloops;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
    cw_assert(a_node < a_tr->ntrns);

    trn = &a_tr->trns[a_node];

    cw_assert(trn->magic == CW_TRN_MAGIC);

    for (i = nneighbors = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE)
	{
	    nneighbors++;

	    /* Make sure that all neighbors point back to this node. */
	    for (j = nloops = 0; j < CW_TR_NODE_MAX_NEIGHBORS; j++)
	    {
		if (a_tr->trns[trn->neighbors[i]].neighbors[j] == a_node)
		{
		    nloops++;
		}
	    }
	    cw_assert(nloops == 1);
	}
    }

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	/* Only leaf nodes can have taxon numbers.  Leaf nodes have at most
	 * 1 neighbor. */
	cw_assert(nneighbors <= 1);
    }

    return TRUE;
}
#endif

/******************************************************************************/

/* tr_ps. */

CW_P_INLINE cw_tr_ps_t *
tr_p_ps_new(cw_tr_t *a_tr)
{
    cw_tr_ps_t *retval;

    retval = (cw_tr_ps_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					    mema_arg_get(a_tr->mema),
					    sizeof(cw_tr_ps_t));

    retval->parent = CW_TR_NODE_NONE;
    retval->chars = NULL;
    retval->nchars = 0;

    return retval;
}

CW_P_INLINE void
tr_p_ps_delete(cw_tr_t *a_tr, cw_tr_ps_t *a_ps)
{
    if (a_ps->chars != NULL)
    {
	cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			  mema_arg_get(a_tr->mema),
			  a_ps->chars, sizeof(cw_uint32_t) * a_ps->nchars);
    }

    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
		      mema_arg_get(a_tr->mema),
		      a_ps, sizeof(cw_tr_ps_t));
}

CW_P_INLINE void
tr_p_ps_prepare(cw_tr_t *a_tr, cw_tr_ps_t *a_ps, cw_uint32_t a_nchars)
{
    /* Clean up old character vector if it isn't the right size for a_nchars
     * characters. */
    if (a_ps->chars != NULL && a_ps->nchars != a_nchars)
    {
	cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			  mema_arg_get(a_tr->mema),
			  a_ps->chars, sizeof(cw_uint32_t) * a_ps->nchars);
	a_ps->chars = NULL;
    }

    /* Allocate character vector if necessary. */
    if (a_ps->chars == NULL)
    {
	a_ps->chars
	    = (cw_uint32_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					      mema_arg_get(a_tr->mema),
					      sizeof(cw_uint32_t) * a_nchars);
	a_ps->nchars = a_nchars;
    }
}

/******************************************************************************/

/* tr_node. */

CW_P_INLINE void
tr_p_node_init(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_trn_t *trn;
    cw_uint32_t i;

    trn = &a_tr->trns[a_node];

    trn->aux = NULL;
    trn->taxon_num = CW_TR_NODE_TAXON_NONE;

    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	trn->neighbors[i] = CW_TR_NODE_NONE;
    }

    trn->ps = NULL;

#ifdef CW_DBG
    trn->magic = CW_TRN_MAGIC;
#endif
}

CW_P_INLINE cw_tr_node_t
tr_p_node_alloc(cw_tr_t *a_tr)
{
    cw_tr_node_t retval;

    if (a_tr->spares == CW_TR_NODE_NONE)
    {
	cw_uint32_t i, nspares;

	/* Allocate spares. */
	cw_check_ptr(a_tr->trns);
	a_tr->trns
	    = (cw_trn_t *) cw_opaque_realloc(mema_realloc_get(a_tr->mema),
					     mema_arg_get(a_tr->mema),
					     a_tr->trns,
					     sizeof(cw_trn_t)
					     * a_tr->ntrns * 2,
					     sizeof(cw_trn_t)
					     * a_tr->ntrns);
	nspares = a_tr->ntrns;
	a_tr->ntrns *= 2;

	/* Initialize last spare. */
	a_tr->spares = a_tr->ntrns - 1;
	a_tr->trns[a_tr->spares].neighbors[0] = CW_TR_NODE_NONE;

	/* Insert other spares into spares stack. */
	for (i = 1; i < nspares; i++)
	{
	    a_tr->spares--;
	    a_tr->trns[a_tr->spares].neighbors[0] = a_tr->spares + 1;
	}
    }

    /* Remove a spare from the spares stack. */
    retval = a_tr->spares;
    a_tr->spares = a_tr->trns[retval].neighbors[0];

    /* Initialize retval. */
    tr_p_node_init(a_tr, retval);

    return retval;
}

CW_P_INLINE void
tr_p_node_dealloc(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_trn_t *trn;

    trn = &a_tr->trns[a_node];
    
    if (trn->ps != NULL)
    {
	tr_p_ps_delete(a_tr, trn->ps);
    }

#ifdef CW_DBG
    memset(&a_tr->trns[a_node], 0x5a, sizeof(cw_trn_t));
#endif

    a_tr->trns[a_node].neighbors[0] = a_tr->spares;
    a_tr->spares = a_node;
}

cw_tr_node_t
tr_node_new(cw_tr_t *a_tr)
{
    cw_tr_node_t retval;
    cw_trn_t *trn;
    cw_uint32_t i;

    cw_dassert(tr_p_validate(a_tr));

    retval = tr_p_node_alloc(a_tr);
    trn = &a_tr->trns[retval];

    trn->aux = NULL;
    trn->taxon_num = CW_TR_NODE_TAXON_NONE;

    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	trn->neighbors[i] = CW_TR_NODE_NONE;
    }

    trn->ps = NULL;

#ifdef CW_DBG
    trn->magic = CW_TRN_MAGIC;
#endif

    return retval;
}

void
tr_node_delete(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    if (a_tr->trns[a_node].ps != NULL)
    {
	tr_p_ps_delete(a_tr, a_tr->trns[a_node].ps);
    }

    tr_p_node_dealloc(a_tr, a_node);
}

cw_uint32_t
tr_node_taxon_num_get(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].taxon_num;
}

void
tr_node_taxon_num_set(cw_tr_t *a_tr, cw_tr_node_t a_node,
		      cw_uint32_t a_taxon_num)
{
    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    a_tr->trns[a_node].taxon_num = a_taxon_num;

    a_tr->modified = TRUE;
}

cw_tr_node_t
tr_node_neighbor_get(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_uint32_t a_i)
{
    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].neighbors[a_i];
}

void
tr_node_neighbors_swap(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_uint32_t a_i,
		       cw_uint32_t a_j)
{
    cw_tr_node_t t_node;

    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_node));
    cw_assert(a_i < CW_TR_NODE_MAX_NEIGHBORS);
    cw_assert(a_j < CW_TR_NODE_MAX_NEIGHBORS);
    cw_assert(a_i != a_j);

    t_node = a_tr->trns[a_node].neighbors[a_i];
    a_tr->trns[a_node].neighbors[a_i] = a_tr->trns[a_node].neighbors[a_j];
    a_tr->trns[a_node].neighbors[a_j] = t_node;

    a_tr->modified = TRUE;
}

void
tr_node_join(cw_tr_t *a_tr, cw_tr_node_t a_a, cw_tr_node_t a_b)
{
    cw_trn_t *trn_a, *trn_b;
    cw_uint32_t i, j;

    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_a));
    cw_dassert(tr_p_node_validate(a_tr, a_b));
    cw_assert(a_a != a_b);

    trn_a = &a_tr->trns[a_a];
    trn_b = &a_tr->trns[a_b];

#ifdef CW_DBG
    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	cw_assert(trn_a->neighbors[i] != a_b);
	cw_assert(trn_b->neighbors[i] != a_a);
    }
#endif

    /* Find an empty slot in a_a. */
    for (i = 0; trn_a->neighbors[i] != CW_TR_NODE_NONE; i++)
    {
	cw_assert(i < CW_TR_NODE_MAX_NEIGHBORS);
    }
    
    /* Find an empty slot in a_b. */
    for (j = 0; trn_b->neighbors[j] != CW_TR_NODE_NONE; j++)
    {
	cw_assert(j < CW_TR_NODE_MAX_NEIGHBORS);
    }

    /* Join the two nodes. */
    trn_a->neighbors[i] = a_b;
    trn_b->neighbors[j] = a_a;

    a_tr->modified = TRUE;

    cw_dassert(tr_p_node_validate(a_tr, a_a));
    cw_dassert(tr_p_node_validate(a_tr, a_b));
}

void
tr_node_detach(cw_tr_t *a_tr, cw_tr_node_t a_a, cw_tr_node_t a_b)
{
    cw_trn_t *trn_a, *trn_b;
    cw_uint32_t i, j;

    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_a));
    cw_dassert(tr_p_node_validate(a_tr, a_b));

    trn_a = &a_tr->trns[a_a];
    trn_b = &a_tr->trns[a_b];

    /* Find the slot in a_a that points to a_b. */
    for (i = 0; trn_a->neighbors[i] != a_b; i++)
    {
	cw_assert(i < CW_TR_NODE_MAX_NEIGHBORS);
    }

    /* Find the slot in a_b that points to a_a. */
    for (j = 0; trn_b->neighbors[j] != a_a; j++)
    {
	cw_assert(j < CW_TR_NODE_MAX_NEIGHBORS);
    }

    /* Detach the two nodes. */
    trn_a->neighbors[i] = CW_TR_NODE_NONE;
    trn_b->neighbors[j] = CW_TR_NODE_NONE;

    a_tr->modified = TRUE;

    cw_dassert(tr_p_node_validate(a_tr, a_a));
    cw_dassert(tr_p_node_validate(a_tr, a_b));
}

void *
tr_node_aux_get(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].aux;
}

void
tr_node_aux_set(cw_tr_t *a_tr, cw_tr_node_t a_node, void *a_aux)
{
    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    a_tr->trns[a_node].aux = a_aux;
}

/******************************************************************************/

/* tr. */

/* Recursively traverse the tree and find the lowest numbered taxon. */
static cw_tr_node_t
tr_p_root_get(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
	      cw_tr_node_t a_root)
{
    cw_tr_node_t retval, root, troot;
    cw_trn_t *trn;
    cw_uint32_t i;

    cw_assert(a_node != CW_TR_NODE_NONE);

    trn = &a_tr->trns[a_node];

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE
	&& (a_root == CW_TR_NODE_NONE
	    || trn->taxon_num < a_tr->trns[a_root].taxon_num))
    {
	retval = a_node;
	root = a_node;
    }
    else
    {
	retval = CW_TR_NODE_NONE;
	root = a_root;
    }

    /* Iterate over neighbors. */
    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE && trn->neighbors[i] != a_prev)
	{
	    troot = tr_p_root_get(a_tr, trn->neighbors[i], a_node, root);
	    if (troot != CW_TR_NODE_NONE)
	    {
		retval = troot;
		root = troot;
	    }
	}
    }

    return retval;
}

/* Recursively traverse the tree, count the number of taxa, and find the lowest
 * numbered taxon. */
static cw_tr_node_t
tr_p_update_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
		    cw_uint32_t *r_ntaxa, cw_tr_node_t a_root)
{
    cw_tr_node_t retval, root, troot;
    cw_trn_t *trn;
    cw_uint32_t i;

    cw_assert(a_node != CW_TR_NODE_NONE);

    trn = &a_tr->trns[a_node];

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	/* Leaf node. */
	(*r_ntaxa)++;
    }

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE
	&& (a_root == CW_TR_NODE_NONE
	    || trn->taxon_num < a_tr->trns[a_root].taxon_num))
    {
	retval = a_node;
	root = a_node;
    }
    else
    {
	retval = CW_TR_NODE_NONE;
	root = a_root;
    }

    /* Iterate over neighbors. */
    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE && trn->neighbors[i] != a_prev)
	{
	    troot = tr_p_update_recurse(a_tr, trn->neighbors[i], a_node,
					r_ntaxa, root);
	    if (troot != CW_TR_NODE_NONE)
	    {
		retval = troot;
		root = troot;
	    }
	}
    }

    return retval;
}

#ifdef CW_DBG
/* Determine whether a_other is reachable from a_node.  This function helps to
 * make sure that no cycles are introduced. */
static cw_bool_t
tr_p_reachable(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
	       cw_tr_node_t a_other)
{
    cw_uint32_t retval, i;
    cw_trn_t *trn;

    trn = &a_tr->trns[a_node];

    if (a_node == a_other)
    {
	retval = TRUE;
	goto RETURN;
    }

    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE && trn->neighbors[i] != a_prev)
	{
	    if ((retval = tr_p_reachable(a_tr, trn->neighbors[i], a_node,
					 a_other)))
	    {
		goto RETURN;
	    }
	}
    }

    retval = FALSE;
    RETURN:
    return retval;
}

/* Return the number of taxa with number a_taxon_num in the subtree rooted at
 * a_node. */
static cw_uint32_t
tr_p_validate_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
		      cw_uint32_t a_taxon_num)
{
    cw_uint32_t retval, i;
    cw_trn_t *trn;

    tr_p_node_validate(a_tr, a_node);

    trn = &a_tr->trns[a_node];

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	cw_uint32_t nneighbors;

	/* Leaf node. */
	cw_assert(trn->neighbors[i] != CW_TR_NODE_NONE);
	for (i = nneighbors = 1; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
	{
	    if (trn->neighbors[i] != CW_TR_NODE_NONE)
	    {
		nneighbors++;
	    }
	}
	cw_assert(nneighbors == 1);

	if (trn->taxon_num == a_taxon_num)
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
	retval = 0;
    }

    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE
	    && trn->neighbors[i] != a_prev)
	{
	    retval += tr_p_validate_recurse(a_tr, trn->neighbors[i], a_node,
					    a_taxon_num);
	}
    }

    return retval;
}

/* Validate a tree. */
static cw_bool_t
tr_p_validate(cw_tr_t *a_tr)
{
    cw_uint32_t i;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);

    if (a_tr->modified == FALSE)
    {
	cw_uint32_t ntaxa;

	ntaxa = 0;
	if (a_tr->base != CW_TR_NODE_NONE)
	{
	    tr_p_update_recurse(a_tr, a_tr->base, CW_TR_NODE_NONE, &ntaxa,
				CW_TR_NODE_NONE);
	}
	cw_assert(a_tr->ntaxa == ntaxa);
    }

    /* Iterate over trns and do some basic sanity checks. */
    for (i = 0; i < a_tr->ntrns; i++)
    {
	if (a_tr->trns[i].magic == CW_TRN_MAGIC)
	{
	    tr_p_node_validate(a_tr, (cw_tr_node_t) i);
	}
	else
	{
	    /* Make sure there are no valid trn's in the free list. */
	    cw_assert(a_tr->trns[i].neighbors[0] == CW_TR_NODE_NONE
		      || a_tr->trns[a_tr->trns[i].neighbors[0]].magic
		      != CW_TRN_MAGIC);
	}
    }

    cw_assert(a_tr->spares == CW_TR_NODE_NONE
	      || a_tr->trns[a_tr->spares].magic != CW_TRN_MAGIC);

    return TRUE;
}
#endif

CW_P_INLINE void
tr_p_edge_get(cw_tr_t *a_tr, cw_uint32_t a_edge, cw_tr_node_t *r_node_a,
	      cw_tr_node_t *r_node_b)
{
    *r_node_a = a_tr->tres[a_edge].node_a;
    *r_node_b = a_tr->tres[a_edge].node_b;
}

/* Pretend that the tree is bisected a the edge between a_node and a_other.
 * Count the number of edges that are in the subtree that contains a_node.  Also
 * return the edge index for the bisection.
 *
 * This function counts both non-bisection edges attached to the node adjacent
 * to the bisection edge, so the final value of *r_edge_count must be
 * decremented, if greater than 0, in order to get the edge count, were the tree
 * actually bisected. */
static void
tr_p_bisection_edge_get_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node,
				cw_tr_node_t a_other, cw_tr_node_t a_prev,
				cw_uint32_t *r_edge_count,
				cw_uint32_t *r_bisection_edge)
{
    cw_uint32_t i, prev_edge_count;
    cw_trn_t *trn;

    cw_assert(a_node != CW_TR_NODE_NONE);

    /* Save the previous edge count, in case it ends up being the index of the
     * edge adjacent to the bisection. */
    prev_edge_count = *r_edge_count;

    trn = &a_tr->trns[a_node];

    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] == a_other)
	{
	    /* Store the index of the edge adjacent to the bisection. */
	    if (prev_edge_count > 0)
	    {
		*r_bisection_edge = prev_edge_count - 1;
	    }
	    else
	    {
		*r_bisection_edge = CW_TR_NODE_EDGE_NONE;
	    }
	}
	else if (trn->neighbors[i] != CW_TR_NODE_NONE
		 && trn->neighbors[i] != a_prev)
	{
	    /* Increment edge count before recursing. */
	    (*r_edge_count)++;

	    /* Recurse into neighbor subtrees. */
	    tr_p_bisection_edge_get_recurse(a_tr, trn->neighbors[i], a_other,
					    a_node, r_edge_count,
					    r_bisection_edge);
	}
    }
}

static void
tr_p_bisection_edges_get(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_uint32_t a_edge,
			 cw_trt_t *r_trt)
{
    cw_tr_node_t tnode;

    /* Count the number of edges that would be in each half of the tree, were it
     * bisected.  Also, determine which edges of the subtrees would be used to
     * reverse the bisection. */

    /* Get the nodes adjacent to the bisection edge. */
    tr_p_edge_get(a_tr, a_edge, &r_trt->adj_a, &r_trt->adj_b);

    /* Get the root nodes (lowest-numbered taxa) for the two subtrees. */
    r_trt->root_a = tr_p_root_get(a_tr, r_trt->adj_a, r_trt->adj_b,
				  CW_TR_NODE_NONE);
    cw_assert(r_trt->root_a != CW_TR_NODE_NONE);
    cw_assert(a_tr->trns[r_trt->root_a].taxon_num != CW_TR_NODE_TAXON_NONE);

    r_trt->root_b = tr_p_root_get(a_tr, r_trt->adj_b, r_trt->adj_a,
				  CW_TR_NODE_NONE);
    cw_assert(r_trt->root_b != CW_TR_NODE_NONE);
    cw_assert(a_tr->trns[r_trt->root_b].taxon_num != CW_TR_NODE_TAXON_NONE);

    /* Make sure that root_a is less than root_b. */
    if (a_tr->trns[r_trt->root_a].taxon_num
	> a_tr->trns[r_trt->root_b].taxon_num)
    {
	tnode = r_trt->root_a;
	r_trt->root_a = r_trt->root_b;
	r_trt->root_b = tnode;

	tnode = r_trt->adj_a;
	r_trt->adj_a = r_trt->adj_b;
	r_trt->adj_b = tnode;
    }

    /* Get the number of edges in the first half of the bisection, as well as
     * the index of the edge adjacent to the bisection. */
    r_trt->nedges_a = 0;
    tr_p_bisection_edge_get_recurse(a_tr, r_trt->root_a, r_trt->adj_b,
				    CW_TR_NODE_NONE, &r_trt->nedges_a,
				    &r_trt->self_a);
    if (r_trt->nedges_a > 0)
    {
	/* Don't count both edges adjacent to the bisection. */
	r_trt->nedges_a--;
    }

    /* Get the number of edges in the second half of the bisection, as well as
     * the index of the edge adjacent to the bisection. */
    r_trt->nedges_b = 0;
    tr_p_bisection_edge_get_recurse(a_tr, r_trt->root_b, r_trt->adj_a,
				    CW_TR_NODE_NONE, &r_trt->nedges_b,
				    &r_trt->self_b);
    if (r_trt->nedges_b > 0)
    {
	/* Don't count both edges adjacent to the bisection. */
	r_trt->nedges_b--;
    }
}

static void
tr_p_ntaxa_nedges_update(cw_tr_t *a_tr)
{
    cw_uint32_t ntaxa;

    /* Update ntaxa and nedges. */
    ntaxa = 0;
    if (a_tr->base != CW_TR_NODE_NONE)
    {
	tr_p_update_recurse(a_tr, a_tr->base, CW_TR_NODE_NONE, &ntaxa,
			    CW_TR_NODE_NONE);
    }

    a_tr->ntaxa = ntaxa;
    if (ntaxa > 1)
    {
	a_tr->nedges = (ntaxa * 2) - 3;
    }
    else
    {
	a_tr->nedges = 0;
    }
}

static void
tr_p_trt_update(cw_tr_t *a_tr, cw_uint32_t a_nedges_prev)
{
    cw_uint32_t i, j, n, offset, a, b;

    cw_assert(a_tr->modified == FALSE);

    /* Allocate/reallocate/deallocate trt. */
    if (a_tr->trt == NULL)
    {
	/* Allocate trt. */
	a_tr->trt = (cw_trt_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
						 mema_arg_get(a_tr->mema),
						 sizeof(cw_trt_t)
						 * (a_tr->nedges + 1));
    }
    else if (a_tr->nedges != a_nedges_prev)
    {
	if (a_tr->nedges > 0)
	{
	    /* Reallocate trt. */
	    a_tr->trt
		= (cw_trt_t *) cw_opaque_realloc(mema_realloc_get(a_tr->mema),
						 mema_arg_get(a_tr->mema),
						 a_tr->trt,
						 sizeof(cw_trt_t)
						 * (a_tr->nedges + 1),
						 sizeof(cw_trt_t)
						 * (a_nedges_prev + 1));
	}
	else
	{
	    /* Deallocate trt. */
	    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			      mema_arg_get(a_tr->mema),
			      a_tr->trt,
			      sizeof(cw_trt_t) * (a_nedges_prev + 1));
	    a_tr->trt = NULL;
	}
    }

    /* Iteratively fill in trt. */
    for (i = j = offset = 0; i < a_tr->nedges; i++)
    {
	/* Record offset. */
	a_tr->trt[j].offset = offset;

	/* Recorde bisection edge. */
	a_tr->trt[j].bisect_edge = i;

	/* Prepare to record number of subtree edges. */
	a_tr->trt[j].nedges_a = 0;
	a_tr->trt[j].nedges_b = 0;

	/* Set trt[i].{nedges,self}_[ab]. */
	tr_p_bisection_edges_get(a_tr, a_tr->base, i, &a_tr->trt[j]);

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
     * not be full.  Keep track of the number of valid elements (not counting
     * the trailing one that stores the total number of TBR neighbors). */
    a_tr->trtused = j;
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
tr_p_tre_update_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
			cw_uint32_t *ar_edge_count)
{
    cw_trn_t *trn, *ttrn;
    cw_uint32_t i, j;

    trn = &a_tr->trns[a_node];

    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE && trn->neighbors[i] != a_prev)
	{
	    /* Record trns edge information. */
	    trn->edges[i] = *ar_edge_count;

	    ttrn = &a_tr->trns[trn->neighbors[i]];
	    for (j = 0; j < CW_TR_NODE_MAX_NEIGHBORS; j++)
	    {
		if (ttrn->neighbors[j] == a_node)
		{
		    ttrn->edges[j] = *ar_edge_count;
		    break;
		}
	    }

	    /* Record tres edge information. */
	    a_tr->tres[*ar_edge_count].node_a = a_node;
	    a_tr->tres[*ar_edge_count].node_b = trn->neighbors[i];

	    /* Increment edge count before recursing. */
	    (*ar_edge_count)++;

	    /* Recurse into neighbor subtree. */
	    tr_p_tre_update_recurse(a_tr, trn->neighbors[i], a_node,
				    ar_edge_count);
	}
    }
}

static void
tr_p_tre_update(cw_tr_t *a_tr, cw_uint32_t a_nedges_prev)
{
    cw_uint32_t edge_count;

    cw_assert(a_tr->modified == FALSE);

    /* Make sure that the tres array is the right size. */
    if (a_tr->nedges > a_nedges_prev)
    {
	if (a_tr->tres == NULL)
	{
	    cw_assert(a_nedges_prev == 0);

	    a_tr->tres
		= (cw_tre_t *) cw_opaque_calloc(mema_calloc_get(a_tr->mema),
						mema_arg_get(a_tr->mema),
						a_tr->nedges,
						sizeof(cw_tre_t));
	}
	else
	{
	    cw_uint32_t i;

	    a_tr->tres
		= (cw_tre_t *) cw_opaque_realloc(mema_realloc_get(a_tr->mema),
						 mema_arg_get(a_tr->mema),
						 a_tr->tres,
						 sizeof(cw_tre_t)
						 * a_tr->nedges,
						 sizeof(cw_tre_t)
						 * a_nedges_prev);
	    memset(&a_tr->tres[a_nedges_prev], 0,
		   (a_tr->nedges - a_nedges_prev) * sizeof(cw_tre_t));

	    /* Initialize ps pointers for newly allocated tre's. */
	    for (i = a_nedges_prev; i < a_tr->nedges; i++)
	    {
		a_tr->tres[i].ps = NULL;
	    }
	}
    }
    else if (a_tr->nedges < a_nedges_prev)
    {
	cw_uint32_t i;

	/* Shrink the array, but first clean up ps's for the tre's at the
	 * end. */
	for (i = a_tr->nedges; i < a_nedges_prev; i++)
	{
	    if (a_tr->tres[i].ps != NULL)
	    {
		tr_p_ps_delete(a_tr, a_tr->tres[i].ps);
	    }
	}

	if (a_tr->nedges > 0)
	{
	    a_tr->tres
		= (cw_tre_t *) cw_opaque_realloc(mema_realloc_get(a_tr->mema),
						 mema_arg_get(a_tr->mema),
						 a_tr->tres,
						 sizeof(cw_tre_t)
						 * a_tr->nedges,
						 sizeof(cw_tre_t)
						 * a_nedges_prev);
	}
	else
	{
	    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			      mema_arg_get(a_tr->mema),
			      a_tr->tres,
			      sizeof(cw_tre_t) * a_nedges_prev);
	    a_tr->tres = NULL;
	}
    }

    /* Recursively traverse the tree, and initialize tres and edges in trns
     * along the way. */
    if (a_tr->nedges > 0)
    {
	edge_count = 0;
	tr_p_tre_update_recurse(a_tr, a_tr->base, CW_TR_NODE_NONE,
				&edge_count);
	cw_assert(edge_count == a_tr->nedges);
    }
}

CW_P_INLINE void
tr_p_update(cw_tr_t *a_tr)
{
    if (a_tr->modified)
    {
	cw_uint32_t nedges_prev;

	/* Store nedges before updating. */
	nedges_prev = a_tr->nedges;

	/* Update ntaxa and nedges. */
	tr_p_ntaxa_nedges_update(a_tr);

	/* Reset the modified flag. */
	a_tr->modified = FALSE;

	/* Update tr and trt. */
	tr_p_tre_update(a_tr, nedges_prev);
	tr_p_trt_update(a_tr, nedges_prev);
    }
}

// Check for malformed trees during recursion.
/* Convert a tree to canonical form by re-ordering the neighbors array such that
 * subtrees are in increasing order of minimum taxon number contained. */
static cw_uint32_t
tr_p_canonize(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev)
{
    cw_uint32_t retval;
    cw_uint32_t i, j, t;
    cw_uint32_t subtree_mins[CW_TR_NODE_MAX_NEIGHBORS - 1];
    cw_uint32_t subtree_inds[CW_TR_NODE_MAX_NEIGHBORS - 1];
    cw_bool_t swapped;
    cw_trn_t *trn;

    cw_dassert(tr_p_node_validate(a_tr, a_node));

    trn = &a_tr->trns[a_node];

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	/* Leaf node. */
	retval = trn->taxon_num;
    }
    else
    {
	/* Internal node. */
	retval = CW_TR_NODE_TAXON_NONE;
    }

    /* Iteratively canonize subtrees, keeping track of the minimum taxon number
     * seen overall, as well as for each subtree. */
    for (i = j = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE && trn->neighbors[i] != a_prev)
	{
	    /* A neighboring subtree that hasn't been visited yet. */
	    cw_assert(j < (CW_TR_NODE_MAX_NEIGHBORS - 1));
	    subtree_mins[j] = tr_p_canonize(a_tr, trn->neighbors[i], a_node);
	    if (subtree_mins[j] < retval)
	    {
		retval = subtree_mins[j];
	    }
	    subtree_inds[j] = i;
	    j++;
	}
    }

    /* Bubble sort the subtrees.  This algorithm works regardless of the value
     * of CW_TR_NODE_MAX_NEIGHBORS, but for bifurcating trees it only requires a
     * couple of extra branches. */
    do
    {
	swapped = FALSE;

	for (i = 0; i + 1 < j; i++)
	{
	    if (subtree_mins[i] > subtree_mins[i + 1])
	    {
		swapped = TRUE;

		/* Swap subtrees. */
		tr_node_neighbors_swap(a_tr, a_node, subtree_inds[i],
				       subtree_inds[i + 1]);

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

/* Extract the node adjacent to bisection and patch its neighbors together. */
CW_P_INLINE void
tr_p_bisection_patch(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_trn_t *trn;
    cw_uint32_t i;
    cw_tr_node_t a, b;

    trn = &a_tr->trns[a_node];

    cw_assert(trn->taxon_num == CW_TR_NODE_TAXON_NONE);

    /* Get trn's neighbors. */
    for (i = 0, a = b = CW_TR_NODE_NONE; b == CW_TR_NODE_NONE; i++)
    {
	cw_assert(i < CW_TR_NODE_MAX_NEIGHBORS);

	if (trn->neighbors[i] != CW_TR_NODE_NONE)
	{
	    if (a == CW_TR_NODE_NONE)
	    {
		a = trn->neighbors[i];
	    }
	    else
	    {
		b = trn->neighbors[i];
	    }
	}
    }

    /* Detach. */
    tr_node_detach(a_tr, a_node, a);
    tr_node_detach(a_tr, a_node, b);

    /* Join. */
    tr_node_join(a_tr, a, b);
}

/* Bisect a_tr at a_edge, and return the nodes adjacent to the bisection. */
CW_P_INLINE void
tr_p_bisect(cw_tr_t *a_tr, cw_uint32_t a_edge,
	    cw_tr_node_t *r_node_a, cw_tr_node_t *r_node_b)
{
    cw_assert(a_edge < a_tr->nedges);

    /* Get the nodes to either side of the edge where the bisection will be
     * done. */
    *r_node_a = a_tr->tres[a_edge].node_a;
    *r_node_b = a_tr->tres[a_edge].node_b;

#ifdef CW_DBG
    /* Assert that nodes are directly connected.  Node validation makes sure
     * that the connection is bi-directional, if it exists, so only bother
     * checking in one direction. */
    {
	cw_uint32_t i;
	cw_bool_t connected;

	for (i = 0, connected = FALSE; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
	{
	    if (a_tr->trns[*r_node_a].neighbors[i] == *r_node_b)
	    {
		connected = TRUE;
		break;
	    }
	}
	cw_assert(connected);
    }
#endif
    cw_dassert(tr_p_reachable(a_tr, *r_node_a, CW_TR_NODE_NONE, *r_node_b));

    /* Detach the two nodes. */
    tr_node_detach(a_tr, *r_node_a, *r_node_b);
}

/* Determine whether the subtree of a tree bisection is ready for
 * reconnection.  The return value has the following possible meanings:
 *
 *   0 : Not ready.
 *   1 : Ready (only one node in the subtree).
 *   2 : Ready (adjacent to a_reconnect_a).
 *   3 : Ready (adjacent to a_reconnect_b).
 */
CW_P_INLINE cw_uint32_t
tr_p_reconnect_ready(cw_tr_t *a_tr, cw_tr_node_t a_node,
		     cw_uint32_t a_reconnect_a, cw_uint32_t a_reconnect_b)
{
    cw_uint32_t retval;
    cw_uint32_t i, nneighbors, ntaxa;
    cw_trn_t *trn;

    trn = &a_tr->trns[a_node];

    /* Is there only one node in the subtree?
     * Are there only two leaf nodes in the subtree? */
    for (i = nneighbors = ntaxa = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE)
	{
	    nneighbors++;

	    if (a_tr->trns[trn->neighbors[i]].taxon_num
		!= CW_TR_NODE_TAXON_NONE)
	    {
		ntaxa++;
	    }
	}
    }
    if (nneighbors == 0 || ntaxa == 2)
    {
	retval = 1;
	goto RETURN;
    }

    /* Is a reconnection edge adjacent to the bisection? */
    if (a_reconnect_a != CW_TR_NODE_EDGE_NONE
	&& (a_tr->tres[a_reconnect_a].node_a == a_node
	    || a_tr->tres[a_reconnect_a].node_b == a_node))
    {
	retval = 2;
	goto RETURN;
    }

    if (a_reconnect_b != CW_TR_NODE_EDGE_NONE
	&& (a_tr->tres[a_reconnect_b].node_a == a_node
	    || a_tr->tres[a_reconnect_b].node_b == a_node))
    {
	retval = 3;
	goto RETURN;
    }

    retval = 0;
    RETURN:
    return retval;
}

CW_P_INLINE void
tr_p_reconnect_prepare(cw_tr_t *a_tr, cw_tr_node_t a_node,
		       cw_uint32_t a_reconnect, cw_uint32_t a_ready)
{
    if (a_ready == 0 && a_reconnect != CW_TR_NODE_EDGE_NONE)
    {
	cw_tr_node_t a, b;


	tr_p_bisection_patch(a_tr, a_node);
	tr_p_edge_get(a_tr, a_reconnect, &a, &b);

	tr_node_detach(a_tr, a, b);

	tr_node_join(a_tr, a_node, a);
	tr_node_join(a_tr, a_node, b);
    }
}

/* Starting at a_root, recursively iterate over the edges in the subtree on this
 * side of the bisection, and return the edge index of the a_edge'th edge
 * iterated over.  The two non-bisection edges of the node adjacent to the
 * bisection edge are counted as a single edge.
 */
static cw_uint32_t
tr_p_bisection_reconnect_edge_get_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node,
					  cw_tr_node_t a_other,
					  cw_tr_node_t a_prev,
					  cw_uint32_t a_edge,
					  cw_uint32_t *r_edge_index)
{
    cw_uint32_t retval, i;
    cw_trn_t *trn;
    cw_tr_node_t next;
    cw_bool_t adjacent;

    cw_assert(a_node != CW_TR_NODE_NONE);

    trn = &a_tr->trns[a_node];

    /* Find neighboring subtrees to recurse into.  If this node is attached to
     * the bisection edge, do not increment *r_edge_index. */
    for (i = 0, adjacent = FALSE; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] == a_other)
	{
	    adjacent = TRUE;
	}
	else if (trn->neighbors[i] != CW_TR_NODE_NONE
		 && trn->neighbors[i] != a_prev)
	{
	    next = trn->neighbors[i];
	}
    }

    if (adjacent)
    {
	/* Do not increment. */
	retval = tr_p_bisection_reconnect_edge_get_recurse(a_tr, next, a_other,
							   a_node, a_edge,
							   r_edge_index);
	if (retval != CW_TR_NODE_NONE)
	{
	    goto RETURN;
	}
    }
    else
    {
	for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
	{
	    if (trn->neighbors[i] != CW_TR_NODE_NONE
		&& trn->neighbors[i] != a_prev)
	    {
		cw_assert(trn->neighbors[i] != a_other);

		/* Is this the edge we're looking for? */
		if (*r_edge_index == a_edge)
		{
		    retval = trn->edges[i];
		    goto RETURN;
		}

		/* Increment edge count. */
		(*r_edge_index)++;

		/* Recurse into neighbor subtree. */
		retval
		    = tr_p_bisection_reconnect_edge_get_recurse(a_tr,
								trn->
								neighbors[i],
								a_other, a_node,
								a_edge,
								r_edge_index);
		if (retval != CW_TR_NODE_NONE)
		{
		    goto RETURN;
		}
	    }
	}
    }

    retval = CW_TR_NODE_NONE;
    RETURN:
    return retval;
}

CW_P_INLINE cw_uint32_t
tr_p_bisection_reconnect_edge_get(cw_tr_t *a_tr, cw_tr_node_t a_root,
				  cw_tr_node_t a_other, cw_uint32_t a_edge)
{
    cw_uint32_t retval;

    if (a_edge == CW_TR_NODE_EDGE_NONE)
    {
	retval = CW_TR_NODE_EDGE_NONE;
    }
    else
    {
	cw_uint32_t edge_index = 0;

	retval =  tr_p_bisection_reconnect_edge_get_recurse(a_tr, a_root,
							    a_other,
							    CW_TR_NODE_NONE,
							    a_edge,
							    &edge_index);
    }

    return retval;
}

static void
tr_p_mp_prepare_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
			cw_uint8_t *a_taxa[], cw_uint32_t a_ntaxa,
			cw_uint32_t a_nchars)
{
    cw_trn_t *trn;
    cw_tr_node_t node;
    cw_uint32_t i, taxon_num;

    trn = &a_tr->trns[a_node];

    if (trn->ps == NULL)
    {
	trn->ps = tr_p_ps_new(a_tr);
    }
    tr_p_ps_prepare(a_tr, trn->ps, a_nchars);

    /* If this is a leaf node, initialize the character state sets. */
    if ((taxon_num = tr_node_taxon_num_get(a_tr, a_node))
	!= CW_TR_NODE_TAXON_NONE)
    {
	cw_uint8_t *chars;

	cw_assert(taxon_num < a_ntaxa);

	chars = a_taxa[taxon_num];
	for (i = 0; i < a_nchars; i++)
	{
	    switch (chars[i])
	    {
		case 'N':
		case 'n':
		case 'X':
		case 'x':
		{
		    trn->ps->chars[i] = 0xf;
		    break;
		}
		case 'V':
		case 'v':
		{
		    trn->ps->chars[i] = 0xe;
		    break;
		}
		case 'H':
		case 'h':
		{
		    trn->ps->chars[i] = 0xd;
		    break;
		}
		case 'M':
		case 'm':
		{
		    trn->ps->chars[i] = 0xc;
		    break;
		}
		case 'D':
		case 'd':
		{
		    trn->ps->chars[i] = 0xb;
		    break;
		}
		case 'R':
		case 'r':
		{
		    trn->ps->chars[i] = 0xa;
		    break;
		}
		case 'W':
		case 'w':
		{
		    trn->ps->chars[i] = 0x9;
		    break;
		}
		case 'A':
		case 'a':
		{
		    trn->ps->chars[i] = 0x8;
		    break;
		}
		case 'B':
		case 'b':
		{
		    trn->ps->chars[i] = 0x7;
		    break;
		}
		case 'S':
		case 's':
		{
		    trn->ps->chars[i] = 0x6;
		    break;
		}
		case 'Y':
		case 'y':
		{
		    trn->ps->chars[i] = 0x5;
		    break;
		}
		case 'C':
		case 'c':
		{
		    trn->ps->chars[i] = 0x4;
		    break;
		}
		case 'K':
		case 'k':
		{
		    trn->ps->chars[i] = 0x3;
		    break;
		}
		case 'G':
		case 'g':
		{
		    trn->ps->chars[i] = 0x2;
		    break;
		}
		case 'T':
		case 't':
		{
		    trn->ps->chars[i] = 0x1;
		    break;
		}
		case '-':
		{
		    /* Treat gaps as uncertainty.  This isn't the only way to
		     * do things, and may need to be made configurable. */
		    trn->ps->chars[i] = 0xf;
		    break;
		}
		default:
		{
		    cw_not_reached();
		}
	    }
	}
    }

    /* Recurse into subtrees. */
    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	node = tr_node_neighbor_get(a_tr, a_node, i);
	if (node != CW_TR_NODE_NONE && node != a_prev)
	{
	    tr_p_mp_prepare_recurse(a_tr, node, a_node, a_taxa, a_ntaxa,
				    a_nchars);
	}
    }
}

CW_P_INLINE cw_uint32_t
tr_p_mp_pscore(cw_tr_t *a_tr, cw_tr_node_t a_p, cw_tr_node_t a_a,
	       cw_tr_node_t a_b)
{
    cw_tr_ps_t *ps_p, *ps_a, *ps_b;

    ps_p = a_tr->trns[a_p].ps;
    ps_a = a_tr->trns[a_a].ps;
    ps_b = a_tr->trns[a_b].ps;

    /* Calculate sum of subtree scores. */
    ps_p->subtrees_score
	= ps_a->subtrees_score + ps_a->node_score
	+ ps_b->subtrees_score + ps_b->node_score;

    /* Only calculate the parent's node score if the cached value is invalid. */
    if (ps_a->parent != a_p || ps_b->parent != a_p)
    {
	cw_uint32_t i, nchars, ns, a, b, p, c, s;

	/* (Re-)calculate node score. */
	ns = 0;

	/* Reset this node's parent pointer, to keep the parent from using an
	 * invalid cached value. */
	ps_p->parent = CW_TR_NODE_NONE;

	/* Set parent pointers, so that cached values may be used in future
	 * runs. */
	ps_a->parent = a_p;
	ps_b->parent = a_p;

	/* Calculate partial Fitch parsimony scores for each character.  The
	 * code inside the loop is written such that the compiler can optimize
	 * out all branches (at least for ia32). */
	for (i = 0, nchars = ps_p->nchars; i < nchars; i++)
	{
	    a = ps_a->chars[i];
	    b = ps_b->chars[i];

	    p = a & b;
	    s = p ? 0 : 1;
	    c = -s;
	    ps_p->chars[i] = (p | (c & (a | b)));
	    ns += s;
	}

	ps_p->node_score = ns;
    }

    return (ps_p->subtrees_score + ps_p->node_score);
}

CW_P_INLINE void
tr_p_mp_nopscore(cw_tr_t *a_tr, cw_tr_node_t a_p, cw_tr_node_t a_a,
		 cw_tr_node_t a_b)
{
    cw_tr_ps_t *ps_p, *ps_a, *ps_b;

    ps_p = a_tr->trns[a_p].ps;
    ps_a = a_tr->trns[a_a].ps;
    ps_b = a_tr->trns[a_b].ps;

    /* Clear parent pointers if necessary, in order to avoid a situation where
     * more than two nodes end up claiming the same parent, due to a combination
     * of tree transformations and short-circuited MP scores (due to max score
     * being exceeded). */
    if (ps_a->parent != a_p || ps_b->parent != a_p)
    {
	ps_p->parent = CW_TR_NODE_NONE;
	ps_a->parent = CW_TR_NODE_NONE;
	ps_b->parent = CW_TR_NODE_NONE;
    }
}

static void
tr_p_mp_score_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
		      cw_uint32_t a_maxscore, cw_bool_t *ar_maxed)
{
    cw_uint32_t i;
    cw_tr_node_t node, a, b;

    a = CW_TR_NODE_NONE;
    b = CW_TR_NODE_NONE;

    /* Recurse into subtrees. */
    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	node = tr_node_neighbor_get(a_tr, a_node, i);
	if (node != CW_TR_NODE_NONE && node != a_prev)
	{
	    if (a == CW_TR_NODE_NONE)
	    {
		a = node;
	    }
	    else
	    {
		cw_assert(b == CW_TR_NODE_NONE);
		b = node;
	    }

	    tr_p_mp_score_recurse(a_tr, node, a_node, a_maxscore, ar_maxed);
	}
    }

    /* Now calculate this node's partial score (unless this is a leaf node). */
    if (b != CW_TR_NODE_NONE)
    {
	cw_assert(a != CW_TR_NODE_NONE);
	if (*ar_maxed == FALSE)
	{
	    /* Calculate the partial score for this node. */
	    if (tr_p_mp_pscore(a_tr, a_node, a, b) >= a_maxscore)
	    {
		/* Maximum score met or exceeded; prevent further score
		 * calculations. */
		*ar_maxed = TRUE;
	    }
	}
	else
	{
	    /* Clear invalid cached parial score if necessary, but do not
	     * calculate the partial score. */
	    tr_p_mp_nopscore(a_tr, a_node, a, b);
	}
    }
}

static cw_uint32_t
tr_p_mp_score(cw_tr_t *a_tr, cw_tr_ps_t *a_ps, cw_tr_node_t a_node_a,
	      cw_tr_node_t a_node_b, cw_uint32_t a_maxscore)
{
    cw_uint32_t retval;
    cw_trn_t *trn;
    cw_bool_t maxed;

    maxed = FALSE;
    tr_p_mp_score_recurse(a_tr, a_node_a, a_node_b, a_maxscore, &maxed);
    if (maxed)
    {
	retval = CW_TR_MAXSCORE_NONE;
	goto RETURN;
    }

    tr_p_mp_score_recurse(a_tr, a_node_b, a_node_a, a_maxscore, &maxed);
    if (maxed)
    {
	retval = CW_TR_MAXSCORE_NONE;
	goto RETURN;
    }

    /* Initialize the temporary node enough so that it can be used by
     * tr_p_mp_pscore(). */
    a_tr->trns[0].ps = a_ps;

    /* Clear the parent pointers of a_node_[ab], to make sure that the score is
     * actually calculated.  This is necessary since the "root" node is always
     * the temporary node.  One artifact of this is that repeating precisely the
     * same score calculation will always result in the final score being
     * recalculated. */
    a_tr->trns[a_node_a].ps->parent = CW_TR_NODE_NONE;
    a_tr->trns[a_node_b].ps->parent = CW_TR_NODE_NONE;

    /* Calculate the final score, using the temporary node. */
    tr_p_mp_pscore(a_tr, 0, a_node_a, a_node_b);

    /* Add up the final score. */
    retval = trn->ps->subtrees_score + trn->ps->node_score;
    RETURN:
    return retval;
}

cw_tr_t *
tr_new(cw_mema_t *a_mema)
{
    cw_tr_t *retval;
    cw_opaque_alloc_t *alloc;
    void *arg;

    alloc = mema_alloc_get(a_mema);
    arg = mema_arg_get(a_mema);

    retval = (cw_tr_t *) cw_opaque_alloc(alloc, arg, sizeof(cw_tr_t));

    retval->mema = a_mema;
    retval->aux = NULL;
    retval->modified = FALSE;
    retval->base = CW_TR_NODE_NONE;
    retval->ntaxa = 0;
    retval->nedges = 0;
    retval->tres = NULL;
    retval->trt = NULL;
    retval->trtused = 0;

    /* Allocate trns with one node.  This is the spare node, which gets used in
     * tr_p_mp_score(). */
    retval->trns = (cw_trn_t *) cw_opaque_alloc(alloc, arg, sizeof(cw_trn_t));
    retval->ntrns = 1;

    /* Initialize spare node. */
    tr_p_node_init(retval, 0);

    retval->spares = CW_TR_NODE_NONE;

#ifdef CW_DBG
    retval->magic = CW_TR_MAGIC;
#endif

    return retval;
}

void
tr_delete(cw_tr_t *a_tr)
{
    cw_opaque_dealloc_t *dealloc;
    void *arg;

    cw_dassert(tr_p_validate(a_tr));

    dealloc = mema_dealloc_get(a_tr->mema);
    arg = mema_arg_get(a_tr->mema);

    /* This assumes that all nodes are deallocated before tr_delete() is
     * called. */
    cw_opaque_dealloc(dealloc, arg, a_tr->trns, sizeof(cw_trn_t) * a_tr->ntrns);

    if (a_tr->trt != NULL)
    {
	cw_opaque_dealloc(dealloc, arg, a_tr->trt,
			  sizeof(cw_trt_t) * (a_tr->nedges + 1));
    }

    if (a_tr->tres != NULL)
    {
	cw_opaque_dealloc(dealloc, arg, a_tr->tres,
			  sizeof(cw_tre_t) * a_tr->nedges);
    }

    cw_opaque_dealloc(dealloc, arg, a_tr, sizeof(cw_tr_t));
}

cw_uint32_t
tr_ntaxa_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    tr_p_update(a_tr);

    return a_tr->ntaxa;
}

cw_uint32_t
tr_nedges_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    tr_p_update(a_tr);

    return a_tr->nedges;
}

void
tr_edge_get(cw_tr_t *a_tr, cw_uint32_t a_edge, cw_tr_node_t *r_node_a,
	    cw_tr_node_t *r_node_b)
{
    cw_dassert(tr_p_validate(a_tr));
    cw_check_ptr(r_node_a);
    cw_check_ptr(r_node_b);

    tr_p_update(a_tr);

    tr_p_edge_get(a_tr, a_edge, r_node_a, r_node_b);
}

cw_uint32_t
tr_edge_index_get(cw_tr_t *a_tr, cw_tr_node_t a_node_a, cw_tr_node_t a_node_b)
{
    cw_uint32_t retval, i;
    cw_trn_t *trn;

    cw_dassert(tr_p_validate(a_tr));
    cw_dassert(tr_p_node_validate(a_tr, a_node_a));
    cw_dassert(tr_p_node_validate(a_tr, a_node_b));

    tr_p_update(a_tr);

    trn = &a_tr->trns[a_node_a];

    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] == a_node_b)
	{
	    retval = trn->edges[i];
	    goto RETURN;
	}
    }

    retval = CW_TR_NODE_EDGE_NONE;
    RETURN:
    return retval;
}

cw_tr_node_t
tr_base_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->base;
}

void
tr_base_set(cw_tr_t *a_tr, cw_tr_node_t a_base)
{
    cw_dassert(tr_p_validate(a_tr));
#ifdef CW_DBG
    if (a_base != CW_TR_NODE_NONE)
    {
	cw_dassert(tr_p_node_validate(a_tr, a_base));
    }
#endif

    a_tr->base = a_base;

    a_tr->modified = TRUE;
}

void
tr_canonize(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    /* Partially update internal state, if necessary.  Don't bother updating trt
     * or tre yet, since we will invalidate them during canonization. */
    if (a_tr->modified)
    {
	tr_p_ntaxa_nedges_update(a_tr);
	a_tr->modified = FALSE;
    }

    if (a_tr->base != CW_TR_NODE_NONE)
    {
	cw_uint32_t ntaxa;

	/* Set base to be the lowest-numbered taxon. */
	ntaxa = 0;
	a_tr->base = tr_p_update_recurse(a_tr, a_tr->base, CW_TR_NODE_NONE,
					 &ntaxa, CW_TR_NODE_NONE);

	/* Canonize the tree. */
	tr_p_canonize(a_tr, a_tr->base, CW_TR_NODE_NONE);
    }

    /* Reset the modified flag. */
    a_tr->modified = FALSE;

    /* Now update tre and trt. */
    tr_p_tre_update(a_tr, a_tr->nedges);
    tr_p_trt_update(a_tr, a_tr->nedges);
}

void
tr_tbr(cw_tr_t *a_tr, cw_uint32_t a_bisect, cw_uint32_t a_reconnect_a,
       cw_uint32_t a_reconnect_b)
{
    cw_tr_node_t node_a, node_b;
    cw_bool_t ready_a, ready_b;

    cw_dassert(tr_p_validate(a_tr));

    tr_p_update(a_tr);

    /* Bisect. */
    tr_p_bisect(a_tr, a_bisect, &node_a, &node_b);

    /* For each subtree, move the node adjacent to the bisection to the
     * reconnection edge.  However, there are three case for which no changes to
     * a subtree are necessary:
     *
     *   1) There is only one node in the subtree.
     *
     *   2) There are only two leaf nodes in the subtree.
     *
     *   3) The reconnection edge is adjacent to the bisection.
     *
     * We don't know which subtree contains which reconnection edge, and
     * figuring this out beforehand would require traversing one of the
     * subtrees.  To avoid that potentially expensive traversal, instead check
     * that none of the above 3 cases apply to a subtree, then use a subtree's
     * node that is adjacent to the bisection as a spare.  It isn't important
     * which subtree the spare comes from, as long as the node is truly a spare.
     */
    ready_a = tr_p_reconnect_ready(a_tr, node_a, a_reconnect_a, a_reconnect_b);
    if (ready_a != 0)
    {
	/* One or the other of the two subtrees must be changed; otherwise the
	 * TBR has no effect on the tree topology. */
	cw_assert(tr_p_reconnect_ready(a_tr, node_b, a_reconnect_a,
				       a_reconnect_b) == 0);

	ready_b = 0;
    }
    else
    {
	ready_b = tr_p_reconnect_ready(a_tr, node_b, a_reconnect_a,
				       a_reconnect_b);
    }

    /* If one of the reconnection edges is adjacent to the bisection, make sure
     * that node_[ab] corresponds to a_reconnect_[ab].  Correspondence doesn't
     * matter otherwise. */
    if (ready_a == 3
	|| ready_b == 2
	|| (ready_a == 1 && a_reconnect_a != CW_TR_NODE_EDGE_NONE)
	|| (ready_b == 1 && a_reconnect_b != CW_TR_NODE_EDGE_NONE)
	)
    {
	cw_uint32_t treconnect;

	treconnect = a_reconnect_a;
	a_reconnect_a = a_reconnect_b;
	a_reconnect_b = treconnect;

	/* At this point, ready_[ab] are only useful for testing zero/non-zero
	 * status. */
    }

    /* Prepare the reconnection edges. */
    tr_p_reconnect_prepare(a_tr, node_a, a_reconnect_a, ready_a);
    tr_p_reconnect_prepare(a_tr, node_b, a_reconnect_b, ready_b);

    /* Reconnect. */
    tr_node_join(a_tr, node_a, node_b);

    /* All changes since the last tr_p_update() call were related to TBR, and we
     * know that this does not impact ntaxa or nedges. */
    a_tr->modified = FALSE;

    /* Update tre and trt. */
    tr_p_tre_update(a_tr, a_tr->nedges);
    tr_p_trt_update(a_tr, a_tr->nedges);

    cw_dassert(tr_p_validate(a_tr));
}

cw_uint32_t
tr_tbr_nneighbors_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    tr_p_update(a_tr);

    return a_tr->trt[a_tr->trtused].offset;
}

void
tr_tbr_neighbor_get(cw_tr_t *a_tr, cw_uint32_t a_neighbor,
		    cw_uint32_t *r_bisect, cw_uint32_t *r_reconnect_a,
		    cw_uint32_t *r_reconnect_b)
{
    cw_trt_t key, *trt;
    cw_uint32_t rem, nedges_a, nedges_b, a, b;

    cw_dassert(tr_p_validate(a_tr));

    tr_p_update(a_tr);
    cw_assert(a_neighbor < a_tr->trt[a_tr->trtused].offset);

    /* Get the bisection edge. */
    key.offset = a_neighbor;
    trt = bsearch(&key, a_tr->trt, a_tr->trtused, sizeof(cw_trt_t),
		  tr_p_trt_compare);
    cw_check_ptr(trt);
    *r_bisect = trt->bisect_edge;

    /* Get the reconnection edges.  The indices for a and b are mapped onto the
     * edges of each subtree, in a peculiar fashion.  The actual ordering isn't
     * important, as long as it is consistent (repeatable) and correct (all TBR
     * neighbors are enumerated).  The edge index mapping can be summarized as
     * follows:
     *
     * 1) Start with a full tree.
     *
     * 2) (Pretend to) bisect the tree at the appropriate edge.
     *
     * 3) For the subtree that contains the globally lowest-numbered taxon,
     *    traverse the tree starting at the lowest-numbered taxon, counting
     *    edges along the way, until edge a is reached.  Return the index of
     *    this edge as *r_reconnect_a, where the index is the number that was
     *    set by tr_p_tre_update().
     *
     * 4) For the subtree not considered in step 3, start traversing the tree at
     *    the lowest-numbered taxon in the subtree, counting edges along the
     *    way, until edge b is reached.  Return the index of this edge as
     *    *r_reconnect_b, where the index is the number that was set for that
     *    edge by tr_p_tre_update().
     */

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

	a = CW_TR_NODE_EDGE_NONE;
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
	b = CW_TR_NODE_EDGE_NONE;

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

    /* Convert a and b to actual edge indices. */
    *r_reconnect_a = tr_p_bisection_reconnect_edge_get(a_tr, trt->root_a,
						       trt->adj_b, a);
    *r_reconnect_b = tr_p_bisection_reconnect_edge_get(a_tr, trt->root_b,
						       trt->adj_a, b);
}

void *
tr_aux_get(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->aux;
}

void
tr_aux_set(cw_tr_t *a_tr, void *a_aux)
{
    cw_dassert(tr_p_validate(a_tr));

    a_tr->aux = a_aux;
}

void
tr_mp_prepare(cw_tr_t *a_tr, cw_uint8_t *a_taxa[], cw_uint32_t a_ntaxa,
	      cw_uint32_t a_nchars)
{
    cw_uint32_t i;

    cw_dassert(tr_p_validate(a_tr));

    tr_p_update(a_tr);

    tr_p_mp_prepare_recurse(a_tr, a_tr->base, CW_TR_NODE_NONE, a_taxa, a_ntaxa,
			    a_nchars);

    /* Initialize ps's for tre's. */
    for (i = 0; i < a_tr->nedges; i++)
    {
	tr_p_ps_prepare(a_tr, a_tr->tres[i].ps, a_nchars);
    }
}

cw_uint32_t
tr_mp_score(cw_tr_t *a_tr, cw_uint32_t a_maxscore)
{
    cw_dassert(tr_p_validate(a_tr));
    cw_assert(a_tr->modified == FALSE);

    return tr_p_mp_score(a_tr, a_tr->tres[0].ps,
			 a_tr->trns[a_tr->base].neighbors[0], a_tr->base,
			 a_maxscore);
}
