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
typedef struct cw_trh_s cw_trh_t;

/* Character (in the systematics sense of the word). */
typedef char cw_trc_t;

/* Partial parsimony score information. */
struct cw_tr_ps_s
{
    /* Parent which most recently used this node's partial score when caching
     * its results.  Both children must still point to the parent in order for
     * the cached results to be valid. */
    cw_tr_ps_t *parent;

    /* If this node is adjacent to the bisection edge, then this variable
     * records the node on the other end of the bisection edge.  This is used
     * when deciding whether to push the value of parent down to the child when
     * recursively scoring. */
    cw_trn_t *other;

    /* Sum of the subtree scores, and this node's score, given particular
     * children.  In order for this to be useful, both childrens' parent
     * pointers must still point to this node. */
    cw_uint32_t subtrees_score;
    cw_uint32_t node_score;

    /* chars points to an array of Fitch parsimony state sets.  Each element in
     * the array contains a bitmap representation of a subset of {ACGT} in the 4
     * least significant bits.  T is the least significant bit.  1 means that a
     * nucleotide is in the set.
     *
     * There are nchars character state sets.
     *
     * achars is the actual allocation, which is padded in order to
     * be able to guarantee that chars is 16 byte-aligned. */
    cw_trc_t *chars;
    cw_uint32_t nchars;
    cw_trc_t *achars;
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
    cw_uint32_t pscore;
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
};


/* Held tree. */
struct cw_trh_s
{
    /* Neighbor index for the tree.  This can be passed to tr_tbr_neighbor_get()
     * to get the associated TBR parameters. */
    cw_uint32_t neighbor;

    /* Fitch parsimony score for the neighboring tree. */
    cw_uint32_t score;
};

/* Specifies different tree holding strategies. */
typedef enum
{
    TR_HOLD_BEST,
    TR_HOLD_BETTER,
    TR_HOLD_ALL
} cw_tr_hold_how_t;

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

    /* bedges is an array of edge indices that is used for enumerating the edges
     * on each side of a logical tree bisection.  The first list starts at
     * offset 0 and has nbedges_a elements.  The second list starts at offset
     * nbedges_a and has nbedges_b elements. */
    cw_uint32_t *bedges;
    cw_uint32_t nbedges_a;
    cw_uint32_t nbedges_b;

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

    /* held is an array of held neighbors.  The array is iteratively doubled as
     * necessary.  heldlen is the actual length of the array, and nheld is the
     * number of elements in use. */
    cw_trh_t *held;
    cw_uint32_t heldlen;
    cw_uint32_t nheld;
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

    retval->parent = NULL;
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
			  a_ps->achars, sizeof(cw_trc_t) * (a_ps->nchars * 8));
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
			  a_ps->achars,
			  sizeof(cw_trc_t) * (a_ps->nchars + 8));
	a_ps->chars = NULL;
    }

    /* Allocate character vector if necessary. */
    if (a_ps->chars == NULL)
    {
	a_ps->achars
	    = (cw_trc_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					   mema_arg_get(a_tr->mema),
					   sizeof(cw_trc_t) * (a_nchars + 8));

	/* Make sure that chars is 16 byte-allocated. */
	if ((((unsigned) a_ps->achars) & 0xfU) == 0)
	{
	    a_ps->chars = a_ps->achars;
	}
	else
	{
	    /* All modern systems guarantee at least 8 byte alignment, so assume
	     * that offsetting by 8 bytes is correct. */
	    cw_assert(&a_ps->achars[8]
		      == &a_ps->achars[16 - (((unsigned) a_ps->achars)
					     & 0xfU)]);
	    a_ps->chars = &a_ps->achars[8];
	}

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
    return tr_p_node_alloc(a_tr);
}

void
tr_node_delete(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    tr_p_node_dealloc(a_tr, a_node);
}

cw_uint32_t
tr_node_taxon_num_get(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].taxon_num;
}

void
tr_node_taxon_num_set(cw_tr_t *a_tr, cw_tr_node_t a_node,
		      cw_uint32_t a_taxon_num)
{
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    a_tr->trns[a_node].taxon_num = a_taxon_num;

    a_tr->modified = TRUE;
}

cw_tr_node_t
tr_node_neighbor_get(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_uint32_t a_i)
{
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].neighbors[a_i];
}

void
tr_node_neighbors_swap(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_uint32_t a_i,
		       cw_uint32_t a_j)
{
    cw_tr_node_t t_node;

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
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].aux;
}

void
tr_node_aux_set(cw_tr_t *a_tr, cw_tr_node_t a_node, void *a_aux)
{
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    a_tr->trns[a_node].aux = a_aux;
}

/******************************************************************************/

/* tr. */

/* Initialize everything except trns and spares. */
CW_P_INLINE void
tr_p_new(cw_tr_t *a_tr, cw_mema_t *a_mema)
{
    a_tr->mema = a_mema;
    a_tr->aux = NULL;
    a_tr->modified = FALSE;
    a_tr->base = CW_TR_NODE_NONE;
    a_tr->ntaxa = 0;
    a_tr->nedges = 0;
    a_tr->tres = NULL;
    a_tr->bedges = NULL;
    a_tr->nbedges_a = 0;
    a_tr->nbedges_b = 0;
    a_tr->trt = NULL;
    a_tr->trtused = 0;
    a_tr->held = NULL;
    a_tr->heldlen = 0;
    a_tr->nheld = 1;

#ifdef CW_DBG
    a_tr->magic = CW_TR_MAGIC;
#endif

}

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
    cw_uint32_t i, ntaxa;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
    cw_assert(a_tr->modified == FALSE);

    ntaxa = 0;
    if (a_tr->base != CW_TR_NODE_NONE)
    {
	tr_p_update_recurse(a_tr, a_tr->base, CW_TR_NODE_NONE, &ntaxa,
			    CW_TR_NODE_NONE);
    }
    cw_assert(a_tr->ntaxa == ntaxa);

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

/* Pretend that the tree is bisected at the edge between a_node and a_other.
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

	    /* Recurse into neighbor subtree. */
	    tr_p_bisection_edge_get_recurse(a_tr, trn->neighbors[i], a_other,
					    a_node, r_edge_count,
					    r_bisection_edge);
	}
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
tr_p_bisection_edge_list_gen_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node,
				     cw_tr_node_t a_prev, cw_uint32_t *ar_edges,
				     cw_uint32_t *ar_nedges)
{
    cw_uint32_t i;
    cw_trn_t *trn;

    trn = &a_tr->trns[a_node];

    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE
	    && trn->neighbors[i] != a_prev)
	{
	    /* Add edge to list. */
	    ar_edges[*ar_nedges] = trn->edges[i];
	    (*ar_nedges)++;

	    /* Recurse into neighbor subtree. */
	    tr_p_bisection_edge_list_gen_recurse(a_tr, trn->neighbors[i],
						 a_node, ar_edges, ar_nedges);
	}
    }
}

/* Pretend that the tree is bisected at the edge between a_node and a_other.
 * Construct a list of edges that are in the subtree that contains a_node.
 *
 * The first element in the list is always the edge that is adjacent to the
 * bisection.  This facilitates recognition of reconnections that would reverse
 * bisection. */
CW_P_INLINE void
tr_p_bisection_edge_list_gen(cw_tr_t *a_tr, cw_tr_node_t a_node,
			     cw_tr_node_t a_other, cw_uint32_t *ar_edges,
			     cw_uint32_t *ar_nedges)
{
    cw_trn_t *trn;
    cw_tr_node_t a, b;
    cw_uint32_t i, i_a;
#ifdef CW_DBG
    cw_uint32_t i_b;
#endif

    /* Initialize the length of the list before recursing. */
    *ar_nedges = 0;

    trn = &a_tr->trns[a_node];

    /* Get trn's neighbors. */
    for (i = 0, a = b = CW_TR_NODE_NONE;
	 i < CW_TR_NODE_MAX_NEIGHBORS && b == CW_TR_NODE_NONE;
	 i++)
    {
	cw_assert(i < CW_TR_NODE_MAX_NEIGHBORS);

	if (trn->neighbors[i] != CW_TR_NODE_NONE
	    && trn->neighbors[i] != a_other)
	{
	    if (a == CW_TR_NODE_NONE)
	    {
		a = trn->neighbors[i];
		i_a = i;
	    }
	    else
	    {
		b = trn->neighbors[i];
#ifdef CW_DBG
		i_b = i;
#endif
	    }
	}
    }

    /* Recurse into subtree below a, if this isn't a leaf node. */
    if (a != CW_TR_NODE_NONE)
    {
	/* Add edge to list. */
	ar_edges[*ar_nedges] = trn->edges[i_a];
	(*ar_nedges)++;

	tr_p_bisection_edge_list_gen_recurse(a_tr, a, a_node, ar_edges,
					     ar_nedges);
    }

    /* Recurse into subtree below b, if this isn't a leaf node. */
    if (b != CW_TR_NODE_NONE)
    {
	/* Do not add edge to list, since logically, this edge and the one next
	 * to a are the same for the purposes of TBR. */

	tr_p_bisection_edge_list_gen_recurse(a_tr, b, a_node, ar_edges,
					     ar_nedges);
    }

    /* If the edge list is still empty (bisection edge next to a leaf node), add
     * a single entry to the list. */
    if (*ar_nedges == 0)
    {
	ar_edges[0] = CW_TR_NODE_EDGE_NONE;
	(*ar_nedges)++;
    }
}

/* Generate lists of edges in each half of a logical bisection at edge
 * a_bisect. */
CW_P_INLINE void
tr_p_bedges_gen(cw_tr_t *a_tr, cw_uint32_t a_bisect)
{
    cw_tre_t *tre;

    cw_assert(a_bisect < a_tr->nedges);

    tre = &a_tr->tres[a_bisect];

    tr_p_bisection_edge_list_gen(a_tr, tre->node_a, tre->node_b,
				 a_tr->bedges, &a_tr->nbedges_a);
    tr_p_bisection_edge_list_gen(a_tr, tre->node_b, tre->node_a,
				 &a_tr->bedges[a_tr->nbedges_a],
				 &a_tr->nbedges_b);
}

static void
tr_p_trt_update(cw_tr_t *a_tr, cw_uint32_t a_nedges_prev)
{
    cw_uint32_t i, j, n, offset;

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

	/* Record bisection edge. */
	a_tr->trt[j].bisect_edge = i;

	/* Update offset. */
	tr_p_bedges_gen(a_tr, i);
	n = (a_tr->nbedges_a * a_tr->nbedges_b) - 1;
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

static void
tr_p_bedges_update(cw_tr_t *a_tr, cw_uint32_t a_nedges_prev)
{
    /* Allocate/reallocate/deallocate bedges.  To keep things simple, allocate
     * as big an array as there are edges, even though not quite that many are
     * ever used. */
    if (a_tr->bedges == NULL)
    {
	/* Allocate bedges. */
	a_tr->bedges
	    = (cw_uint32_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					      mema_arg_get(a_tr->mema),
					      sizeof(cw_uint32_t)
					      * a_tr->nedges);
    }
    else if (a_tr->nedges != a_nedges_prev)
    {
	if (a_tr->nedges > 0)
	{
	    /* Reallocate bedges. */
	    a_tr->bedges = (cw_uint32_t *)
		cw_opaque_realloc(mema_realloc_get(a_tr->mema),
				  mema_arg_get(a_tr->mema),
				  a_tr->bedges,
				  sizeof(cw_uint32_t) * a_tr->nedges,
				  sizeof(cw_uint32_t) * a_nedges_prev);
	}
	else
	{
	    /* Deallocate bedges. */
	    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			      mema_arg_get(a_tr->mema),
			      a_tr->bedges,
			      sizeof(cw_uint32_t) * a_nedges_prev);
	    a_tr->bedges = NULL;
	}
    }

    /* Clear nbedges_[ab]. */
    a_tr->nbedges_a = 0;
    a_tr->nbedges_b = 0;
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

	/* Update tre, bedges, and trt. */
	tr_p_tre_update(a_tr, nedges_prev);
	tr_p_bedges_update(a_tr, nedges_prev);
	tr_p_trt_update(a_tr, nedges_prev);

	/* Clear held trees. */
	a_tr->nheld = 0;
    }
}

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
#ifdef CW_DBG
    cw_uint32_t nneighbors = 0;
#endif

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
#ifdef CW_DBG
	if (trn->neighbors[i] != CW_TR_NODE_NONE)
	{
	    nneighbors++;
	}
#endif

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

    cw_dassert((trn->taxon_num != CW_TR_NODE_TAXON_NONE && nneighbors <= 1)
	       || (trn->taxon_num == CW_TR_NODE_TAXON_NONE &&
		   nneighbors == CW_TR_NODE_MAX_NEIGHBORS));

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
    cw_uint32_t i, nneighbors;
    cw_trn_t *trn;

    trn = &a_tr->trns[a_node];

    /* Is there only one node in the subtree?
     * Are there only two leaf nodes in the subtree? */
    for (i = nneighbors = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	if (trn->neighbors[i] != CW_TR_NODE_NONE)
	{
	    nneighbors++;
	}
    }
    if (nneighbors == 0)
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

    /* If this is a leaf node, initialize the character state sets and
     * scores. */
    if ((taxon_num = tr_node_taxon_num_get(a_tr, a_node))
	!= CW_TR_NODE_TAXON_NONE)
    {
	cw_uint8_t *chars;

	cw_assert(taxon_num < a_ntaxa);

	trn->ps->subtrees_score = 0;
	trn->ps->node_score = 0;

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
	node = trn->neighbors[i];
	if (node != CW_TR_NODE_NONE && node != a_prev)
	{
	    tr_p_mp_prepare_recurse(a_tr, node, a_node, a_taxa, a_ntaxa,
				    a_nchars);
	}
    }
}

static void
tr_p_mp_finish_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev)
{
    cw_trn_t *trn;
    cw_tr_node_t node;
    cw_uint32_t i;

    trn = &a_tr->trns[a_node];

    if (trn->ps != NULL)
    {
	tr_p_ps_delete(a_tr, trn->ps);
	trn->ps = NULL;
    }

    /* Recurse into subtrees. */
    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	node = trn->neighbors[i];
	if (node != CW_TR_NODE_NONE && node != a_prev)
	{
	    tr_p_mp_finish_recurse(a_tr, node, a_node);
	}
    }
}

#ifdef CW_CPU_IA32
CW_P_INLINE cw_uint32_t
tr_p_mp_ia32_pscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
		    cw_tr_ps_t *a_b)
{

    /* Only calculate the parent's node score if the cached value is invalid. */
    if (a_a->parent != a_p || a_b->parent != a_p)
    {
	cw_uint32_t i, nchars, ns;
	cw_trc_t *chars_p, *chars_a, *chars_b;

	/* Calculate sum of subtree scores. */
	a_p->subtrees_score
	    = a_a->subtrees_score + a_a->node_score
	    + a_b->subtrees_score + a_b->node_score;

	/* (Re-)calculate node score. */
	ns = 0;

	/* Reset this node's parent pointer, to keep the parent from using an
	 * invalid cached value. */
	a_p->parent = NULL;

	/* Set parent pointers, so that cached values may be used in future
	 * runs. */
	a_a->parent = a_p;
	a_b->parent = a_p;

	/* Calculate partial Fitch parsimony scores for each character. */
	chars_p = a_p->chars;
	chars_a = a_a->chars;
	chars_b = a_b->chars;

	nchars = a_p->nchars;

	/* Initialize SSE2 registers. */
	{
	    static const unsigned char ones[] =
		"\x01\x01\x01\x01\x01\x01\x01\x01"
		"\x01\x01\x01\x01\x01\x01\x01\x01";

	    asm volatile (
		/* Fill xmm7 with 16 1's. */
		"movdqu %[ones], %%xmm7;"

		/* Clear pns. */
		"pxor %%xmm5, %%xmm5;"
		:
		: [ones] "m" (*ones)
		: "%xmm5", "%xmm7"
		);
	}

	/* Use SSE2 to evaluate as many of the characters as possible.  This
	 * loop handles 16 characters per iteration. */
	for (i = 0;
	     i < (nchars ^ (nchars & 0xf));
	     i += 16)
	{
	    asm volatile (
		/* Read character data, and'ing and or'ing them together.
		 *
		 * a = *chars_a;
		 * b = *chars_b;
		 * p = a & b;
		 * d = a | b;
		 */
		"movdqa %[a], %%xmm0;"
		"movdqa %%xmm0, %%xmm1;"
		"por %[b], %%xmm0;" /* xmm0 contains d. */
		"pand %[b], %%xmm1;" /* xmm1 contains p. */

		/* Create bitmasks according to whether the character state sets
		 * are empty.
		 *
		 * c = p ? 0x00 : 0xff;
		 * e = (c & d);
		 * s = p ? 0 : 1;
		 */
		"pxor %%xmm2, %%xmm2;"
		"pcmpeqb %%xmm1, %%xmm2;" /* xmm2 contains c. */
		"pand %%xmm2, %%xmm0;" /* xmm0 contains e. */
		"pand %%xmm7, %%xmm2;" /* xmm2 contains s. */

		/* Update node score.
		 *
		 * ns += s;
		 */
		// XXX In the worst case, this only works 255 times.
		"paddusb %%xmm2, %%xmm5;"

		/* p = (p | e); */
		"por %%xmm1, %%xmm0;" /* xmm0 contains p. */

		/* Store results.
		 *
		 * *chars_p = p;
		 */
		"movdqa %%xmm0, %[p];"
		: [p] "=m" (chars_p[i])
		: [a] "m" (chars_a[i]), [b] "m" (chars_b[i])
		: "memory"
		);
	}

	/* Update ns. */
	{
	    cw_uint32_t j;
	    unsigned char pns[16];

	    asm volatile (
		"movdqu %%xmm5, %[pns];"
		: [pns] "=m" (*pns)
		:
		: "memory"
	    );

	    for (j = 0; j < 16; j++)
	    {
		ns += pns[j];
	    }
	}

	/* Evaluate the last 0-15 characters that weren't evaluated in the above
	 * loop. */
	{
	    cw_uint32_t a, b, p, c, s;

	    for (; i < nchars; i++)
	    {
		a = chars_a[i];
		b = chars_b[i];

		p = a & b;
		s = p ? 0 : 1;
		c = -s;
		chars_p[i] = (p | (c & (a | b)));
		ns += s;
	    }
	}

	a_p->node_score = ns;
    }

    return (a_p->subtrees_score + a_p->node_score);
}
#endif

static cw_uint32_t
tr_p_mp_c_pscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
		 cw_tr_ps_t *a_b)
{
//#define CW_TR_MP_PSCORE_VALIDATE
#ifdef CW_TR_MP_PSCORE_VALIDATE
    cw_bool_t cached;
#endif

    /* Only calculate the parent's node score if the cached value is invalid. */
    if (a_a->parent != a_p || a_b->parent != a_p)
#ifdef CW_TR_MP_PSCORE_VALIDATE
    {
	cached = FALSE;
    }
    else
    {
	cached = TRUE;
    }
#endif
    {
	cw_uint32_t i, nchars, ns, a, b, p, c, s;
	cw_trc_t *chars_p, *chars_a, *chars_b;

#ifdef CW_TR_MP_PSCORE_VALIDATE
	if (cached)
	{
	    if (a_p->subtrees_score
		!= (a_a->subtrees_score + a_a->node_score
		    + a_b->subtrees_score + a_b->node_score))
	    {
		fprintf(stderr,
			"%s:%d:%s(): subtrees_score %u (should be %u)\n",
			__FILE__, __LINE__, __FUNCTION__,
			a_p->subtrees_score,
			a_a->subtrees_score + a_a->node_score
			+ a_b->subtrees_score + a_b->node_score);
		abort();
	    }
	}
#endif
	/* Calculate sum of subtree scores. */
	a_p->subtrees_score
	    = a_a->subtrees_score + a_a->node_score
	    + a_b->subtrees_score + a_b->node_score;

	/* (Re-)calculate node score. */
	ns = 0;

	/* Reset this node's parent pointer, to keep the parent from using an
	 * invalid cached value. */
	a_p->parent = NULL;

	/* Set parent pointers, so that cached values may be used in future
	 * runs. */
	a_a->parent = a_p;
	a_b->parent = a_p;

	/* Calculate partial Fitch parsimony scores for each character.  The
	 * code inside the loop is written such that the compiler can optimize
	 * out all branches (at least for ia32). */
	chars_p = a_p->chars;
	chars_a = a_a->chars;
	chars_b = a_b->chars;
	for (i = 0, nchars = a_p->nchars; i < nchars; i++)
	{
	    a = chars_a[i];
	    b = chars_b[i];

	    p = a & b;
	    s = p ? 0 : 1;
	    c = -s;
	    chars_p[i] = (p | (c & (a | b)));
	    ns += s;
	}

#ifdef CW_TR_MP_PSCORE_VALIDATE
	if (cached)
	{
	    if (ns != a_p->node_score)
	    {
		fprintf(stderr, "%s:%d:%s(): node_score %u (should be %u)\n",
			__FILE__, __LINE__, __FUNCTION__,
			ns, a_p->node_score);
		abort();
	    }
	}
#endif
	a_p->node_score = ns;
    }

    return (a_p->subtrees_score + a_p->node_score);
}

CW_P_INLINE cw_uint32_t
tr_p_mp_pscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b)
{
    cw_uint32_t retval;

#ifdef CW_CPU_IA32
    if (modcrux_ia32_use_sse2)
    {
	retval = tr_p_mp_ia32_pscore(a_tr, a_p, a_a, a_b);
    }
    else
#endif
    {
	retval = tr_p_mp_c_pscore(a_tr, a_p, a_a, a_b);
    }

    return retval;
}

CW_P_INLINE void
tr_p_mp_nopscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
		 cw_tr_ps_t *a_b)
{
    /* Clear parent pointers if necessary, in order to avoid a situation where
     * more than two nodes end up claiming the same parent, due to a combination
     * of tree transformations and short-circuited MP scores (due to max score
     * being exceeded). */
    if (a_a->parent != a_p || a_b->parent != a_p)
    {
	a_p->parent = NULL;
	a_a->parent = NULL;
	a_b->parent = NULL;
    }
}

CW_P_INLINE void
tr_p_mp_passpscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a)
{
    a_p->parent = a_a->parent;
    a_p->subtrees_score = a_a->subtrees_score;
    a_p->node_score = a_a->node_score;
    memcpy(a_p->chars, a_a->chars, sizeof(cw_trc_t) * a_p->nchars);
}

static void
tr_p_mp_score_recurse(cw_tr_t *a_tr, cw_trn_t *a_node, cw_trn_t *a_prev,
		      cw_trn_t *a_other, cw_uint32_t a_maxscore,
		      cw_bool_t *ar_maxed)
{
    cw_uint32_t i;
    cw_trn_t *trn, *a, *b, *o;
    cw_tr_node_t node;

    a = NULL;
    b = NULL;
    o = NULL;

    /* Set a, b, and o, according to which neighbors are subtrees or on the
     * other side of the bisection. */
    for (i = 0; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
    {
	node = a_node->neighbors[i];
	if (node != CW_TR_NODE_NONE)
	{
	    trn = &a_tr->trns[node];
	    if (trn == a_other)
	    {
		o = trn;
	    }
	    else if (trn != a_prev)
	    {
		if (a == NULL)
		{
		    a = trn;
		}
		else
		{
		    cw_assert(b == NULL);
		    b = trn;
		}
	    }
	}
    }

    /* Recursively calculate partial scores for the subtrees, then calculate the
     * partial score for this node. */
    if (b != NULL)
    {
	/* Recurse into subtrees. */
	tr_p_mp_score_recurse(a_tr, a, a_node, a_other, a_maxscore, ar_maxed);
	tr_p_mp_score_recurse(a_tr, b, a_node, a_other, a_maxscore, ar_maxed);

	/* Clear cached other value. */
	a_node->ps->other = o;

	/* Calculate this node's partial score. */

	cw_check_ptr(a);
	if (*ar_maxed == FALSE)
	{
	    /* Calculate the partial score for this node. */
	    if (tr_p_mp_pscore(a_tr, a_node->ps, a->ps, b->ps) >= a_maxscore)
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
	    tr_p_mp_nopscore(a_tr, a_node->ps, a->ps, b->ps);
	}
    }
    else if (o != NULL)
    {
	cw_check_ptr(a);
	cw_assert(b == NULL);

	/* Recurse into subtree, and pass down cached parent value, if the last
	 * time this node was recursed through, o was the same. */
	if (a_node->ps->other == o)
	{
	    /* Pass cached parent value to the child. */
	    a->ps->parent = a_node->ps->parent;
	}
	else
	{
	    /* Clear cached other value. */
	    a_node->ps->other = o;
	}

	/* Recurse into subtree. */
	tr_p_mp_score_recurse(a_tr, a, a_node, a_other, a_maxscore, ar_maxed);

	/* Copy a's scores to this node, rather than calculating a partial
	 * score.  This node is merely a filler node, as far as scoring is
	 * concerned. */
	if (a != NULL)
	{
	    tr_p_mp_passpscore(a_tr, a_node->ps, a->ps);
	}
    }
}

static cw_uint32_t
tr_p_mp_score(cw_tr_t *a_tr, cw_tr_ps_t *a_ps, cw_tr_node_t a_node_a,
	      cw_tr_node_t a_node_b, cw_tr_node_t a_other,
	      cw_uint32_t a_maxscore)
{
    cw_uint32_t retval;
    cw_bool_t maxed;

    maxed = FALSE;
    tr_p_mp_score_recurse(a_tr, &a_tr->trns[a_node_a], &a_tr->trns[a_node_b],
			  &a_tr->trns[a_other], a_maxscore, &maxed);
    if (maxed)
    {
	retval = CW_TR_MAXSCORE_NONE;
	goto RETURN;
    }

    tr_p_mp_score_recurse(a_tr, &a_tr->trns[a_node_b], &a_tr->trns[a_node_a],
			  &a_tr->trns[a_other], a_maxscore, &maxed);
    if (maxed)
    {
	retval = CW_TR_MAXSCORE_NONE;
	goto RETURN;
    }

    /* Clear the parent pointers of a_node_[ab], to make sure that the score is
     * actually calculated.  This is necessary since the "root" node is always
     * the temporary node.  One artifact of this is that repeating precisely the
     * same score calculation will always result in the final score being
     * recalculated. */
    a_tr->trns[a_node_a].ps->parent = NULL;
    a_tr->trns[a_node_b].ps->parent = NULL;

    /* Calculate the final score, using a_ps. */
    tr_p_mp_pscore(a_tr, a_ps, a_tr->trns[a_node_a].ps,
		   a_tr->trns[a_node_b].ps);

    /* Add up the final score. */
    retval = a_ps->subtrees_score + a_ps->node_score;
    RETURN:
    return retval;
}

CW_P_INLINE void
tr_p_bisection_edge_list_mp(cw_tr_t *a_tr, cw_uint32_t *a_edges,
			    cw_uint32_t a_nedges, cw_tr_node_t a_other,
			    cw_uint32_t a_maxscore)
{
    cw_uint32_t i;
    cw_tre_t *tre;

    for (i = 0; i < a_nedges; i++)
    {
	if (a_edges[i] != CW_TR_NODE_EDGE_NONE)
	{
	    tre = &a_tr->tres[a_edges[i]];

	    tre->pscore = tr_p_mp_score(a_tr, tre->ps, tre->node_a, tre->node_b,
					a_other, a_maxscore);
	}
    }
}

/* Hold a tree.  If a_max_held is exceeded, the tree is not held.  This
 * introduces a bias in which trees are held.  There exist algorithms for making
 * this an unbiased process, but there is no need for that functionality at the
 * moment. */
CW_P_INLINE void
tr_p_hold(cw_tr_t *a_tr, cw_uint32_t a_max_hold, cw_uint32_t a_neighbor,
	  cw_uint32_t a_score)
{
    if (a_tr->nheld < a_max_hold)
    {
	cw_trh_t *trh;

	/* Make sure there is space to store another held tree. */
	if (a_tr->held == NULL)
	{
	    /* Allocate. */
	    a_tr->held
		= (cw_trh_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					       mema_arg_get(a_tr->mema),
					       sizeof(cw_trh_t));
	    a_tr->heldlen = 1;
	}
	else if (a_tr->nheld == a_tr->heldlen)
	{
	    /* Reallocate. */
	    a_tr->held
		= (cw_trh_t *) cw_opaque_realloc(mema_realloc_get(a_tr->mema),
						 mema_arg_get(a_tr->mema),
						 a_tr->held,
						 sizeof(cw_trh_t)
						 * a_tr->heldlen * 2,
						 sizeof(cw_trh_t)
						 * a_tr->heldlen);
	    a_tr->heldlen *= 2;
	}
	
	/* Hold this tree. */
	trh = &a_tr->held[a_tr->nheld];
	trh->neighbor = a_neighbor;
	trh->score = a_score;

	a_tr->nheld++;
    }
}

CW_P_INLINE void
tr_p_tbr_neighbors_mp(cw_tr_t *a_tr, cw_uint32_t a_max_hold,
		      cw_uint32_t a_maxscore, cw_tr_hold_how_t a_how)
{
    cw_uint32_t neighbor, i, j, k, edge_a, edge_b;
    cw_uint32_t score;
    cw_tre_t *tre, *tre_a, *tre_b;
    cw_tr_ps_t *ps, *ps_a, *ps_b;

    cw_dassert(tr_p_validate(a_tr));

    /* Set up tree holding data structures. */
    a_tr->nheld = 0;

    /* Iteratively (logically) bisect at each edge in the tree. */
    ps = a_tr->trns[0].ps;
    neighbor = 0;
    for (i = 0; i < a_tr->nedges; i++)
    {
	tre = &a_tr->tres[i];

	/* Determine which edges are in each subtree. */
	tr_p_bedges_gen(a_tr, i);

	cw_assert((a_tr->nbedges_a == 1 && (a_tr->nbedges_b == a_tr->nedges - 2
				     || a_tr->nbedges_b == a_tr->nedges - 4))
		  || (a_tr->nbedges_b == 1
		      && (a_tr->nbedges_a == a_tr->nedges - 2
			  || a_tr->nbedges_a == a_tr->nedges - 4))
		  || (a_tr->nbedges_a != 1 && a_tr->nbedges_b != 1
		      && a_tr->nbedges_a + a_tr->nbedges_b
		      == a_tr->nedges - 3));

	/* Calculate the partial score for each edge in the edge lists. */
	tr_p_bisection_edge_list_mp(a_tr, a_tr->bedges, a_tr->nbedges_a,
				    tre->node_b, a_maxscore);
	tr_p_bisection_edge_list_mp(a_tr, &a_tr->bedges[a_tr->nbedges_a],
				    a_tr->nbedges_b, tre->node_a, a_maxscore);

	/* Iteratively (logically) reconnect every legitimate pairing of edges
	 * between the two subtrees and calculate final parsimony scores. */
	for (j = 0; j < a_tr->nbedges_a; j++)
	{
	    edge_a = a_tr->bedges[j];
	    if (edge_a != CW_TR_NODE_EDGE_NONE)
	    {
		tre_a = &a_tr->tres[edge_a];
		ps_a = a_tr->tres[edge_a].ps;
	    }
	    else
	    {
		tre_a = NULL;
		ps_a = a_tr->trns[tre->node_a].ps;
	    }

	    /* Skip this iteration if edge_a's partial score exceeded
	     * a_maxscore. */
	    if (tre_a != NULL && tre_a->pscore == CW_TR_MAXSCORE_NONE)
	    {
		continue;
	    }

	    for (k = 0; k < a_tr->nbedges_b; k++)
	    {
		edge_b = a_tr->bedges[a_tr->nbedges_a + k];
		if (edge_b != CW_TR_NODE_EDGE_NONE)
		{
		    tre_b = &a_tr->tres[edge_b];
		    ps_b = a_tr->tres[edge_b].ps;
		}
		else
		{
		    tre_b = NULL;
		    ps_b = a_tr->trns[tre->node_b].ps;
		}

		/* Skip this iteration if the reconnection would result in
		 * reversing the bisection, or if edge_b's partial score
		 * exceeded a_maxscore. */
		if ((j == 0 && k == 0)
		    || (tre_b != NULL && tre_b->pscore == CW_TR_MAXSCORE_NONE))
		{
		    continue;
		}

		/* Calculate the final parsimony score for this reconnection.
		 * Clear the parent pointers of ps_[ab], to make sure that the
		 * score is actually calculated. */
		ps_a->parent = NULL;
		ps_b->parent = NULL;
		score = tr_p_mp_pscore(a_tr, ps, ps_a, ps_b);

		/* Hold the tree, if appropriate. */
		switch (a_how)
		{
		    case TR_HOLD_BEST:
		    {
			if (a_tr->held == NULL || score == a_tr->held[0].score)
			{
			    /* No trees held, or this tree is as good as those
			     * currently held. */
			    tr_p_hold(a_tr, a_max_hold, neighbor, score);
			}
			else if (score < a_tr->held[0].score)
			{
			    /* This tree is better than the tree(s) currently
			     * held.  Clear the held trees, then hold this
			     * one. */
			    a_tr->nheld = 0;
			    tr_p_hold(a_tr, a_max_hold, neighbor, score);
			}
			break;
		    }
		    case TR_HOLD_BETTER:
		    {
			if (score < a_maxscore)
			{
			    /* No trees held, or this (neighboring) tree is
			     * better than the tree whose neighbors are being
			     * evaluated. */
			    tr_p_hold(a_tr, a_max_hold, neighbor, score);
			}
			break;
		    }
		    case TR_HOLD_ALL:
		    {
			/* Hold all trees. */
			tr_p_hold(a_tr, a_max_hold, neighbor, score);
			break;
		    }
		    default:
		    {
			cw_not_reached();
		    }
		}

		/* Increment the neighbor index here.  Due to the possibility of
		 * loop short-circuting above, this must happen at the end of
		 * the loop body, rather than in the 'for' statement. */
		neighbor++;
	    }
	}
    }
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

    tr_p_new(retval, a_mema);

    /* Allocate trns with one node.  This is the spare node, which gets used in
     * tr_p_mp_score(). */
    retval->trns = (cw_trn_t *) cw_opaque_alloc(alloc, arg, sizeof(cw_trn_t));
    retval->ntrns = 1;

    /* Initialize spare node. */
    tr_p_node_init(retval, 0);

    retval->spares = CW_TR_NODE_NONE;

    return retval;
}

cw_tr_t *
tr_dup(cw_tr_t *a_tr)
{
    cw_tr_t *retval;
    cw_opaque_alloc_t *alloc;
    void *arg;
    cw_uint32_t i;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    alloc = mema_alloc_get(a_tr->mema);
    arg = mema_arg_get(a_tr->mema);

    retval = (cw_tr_t *) cw_opaque_alloc(alloc, arg, sizeof(cw_tr_t));

    tr_p_new(retval, a_tr->mema);

    /* Allocate trns the same size as a_tr's, then copy. */
    retval->trns = (cw_trn_t *) cw_opaque_alloc(alloc, arg,
						sizeof(cw_trn_t) * a_tr->ntrns);
    memcpy(retval->trns, a_tr->trns, sizeof(cw_trn_t) * a_tr->ntrns);
    retval->ntrns = a_tr->ntrns;

    /* Clean up the copied trn's. */
    for (i = 0; i < retval->ntrns; i++)
    {
	a_tr->trns[i].aux = NULL;
	a_tr->trns[i].ps = NULL;
    }

    /* The spares list is the same as for a_tr. */
    retval->spares = a_tr->spares;

    return retval;
}

void
tr_delete(cw_tr_t *a_tr)
{
    cw_opaque_dealloc_t *dealloc;
    void *arg;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);

    dealloc = mema_dealloc_get(a_tr->mema);
    arg = mema_arg_get(a_tr->mema);

    if (a_tr->held != NULL)
    {
	cw_opaque_dealloc(dealloc, arg, a_tr->held,
			  sizeof(cw_trh_t) * a_tr->heldlen);
    }

    /* Clean up the temporary node. */
    tr_p_node_dealloc(a_tr, 0);

    /* This assumes that all nodes are deallocated before tr_delete() is
     * called. */
    cw_opaque_dealloc(dealloc, arg, a_tr->trns, sizeof(cw_trn_t) * a_tr->ntrns);

    if (a_tr->trt != NULL)
    {
	cw_opaque_dealloc(dealloc, arg, a_tr->trt,
			  sizeof(cw_trt_t) * (a_tr->nedges + 1));
    }

    if (a_tr->bedges != NULL)
    {
	cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			  mema_arg_get(a_tr->mema),
			  a_tr->bedges, sizeof(cw_uint32_t) * a_tr->nedges);
    }

    if (a_tr->tres != NULL)
    {
	cw_uint32_t i;

	for (i = 0; i < a_tr->nedges; i++)
	{
	    tr_p_ps_delete(a_tr, a_tr->tres[i].ps);
	}

	cw_opaque_dealloc(dealloc, arg, a_tr->tres,
			  sizeof(cw_tre_t) * a_tr->nedges);
    }

    cw_opaque_dealloc(dealloc, arg, a_tr, sizeof(cw_tr_t));
}

cw_uint32_t
tr_ntaxa_get(cw_tr_t *a_tr)
{
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->ntaxa;
}

cw_uint32_t
tr_nedges_get(cw_tr_t *a_tr)
{
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->nedges;
}

void
tr_edge_get(cw_tr_t *a_tr, cw_uint32_t a_edge, cw_tr_node_t *r_node_a,
	    cw_tr_node_t *r_node_b)
{
    cw_check_ptr(r_node_a);
    cw_check_ptr(r_node_b);

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    tr_p_edge_get(a_tr, a_edge, r_node_a, r_node_b);
}

cw_uint32_t
tr_edge_index_get(cw_tr_t *a_tr, cw_tr_node_t a_node_a, cw_tr_node_t a_node_b)
{
    cw_uint32_t retval, i;
    cw_trn_t *trn;

    cw_dassert(tr_p_node_validate(a_tr, a_node_a));
    cw_dassert(tr_p_node_validate(a_tr, a_node_b));

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

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
    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);

    return a_tr->base;
}

void
tr_base_set(cw_tr_t *a_tr, cw_tr_node_t a_base)
{
    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
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
    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);

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

    cw_dassert(tr_p_validate(a_tr));
}

void
tr_tbr(cw_tr_t *a_tr, cw_uint32_t a_bisect, cw_uint32_t a_reconnect_a,
       cw_uint32_t a_reconnect_b)
{
    cw_tr_node_t node_a, node_b;
    cw_bool_t ready_a, ready_b;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

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
    if (ready_a == 3 || ready_b == 2)
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
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->trt[a_tr->trtused].offset;
}

void
tr_tbr_neighbor_get(cw_tr_t *a_tr, cw_uint32_t a_neighbor,
		    cw_uint32_t *r_bisect, cw_uint32_t *r_reconnect_a,
		    cw_uint32_t *r_reconnect_b)
{
    cw_trt_t key, *trt;
    cw_uint32_t rem;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));
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
     * 3) For each subtree, do a recursive in-order traversal and build a list
     *    of edges, not including one of the edges adjacent to the bisection (in
     *    the case of an internal node adjacent to the bisection).
     *
     * 4) Use a nested loop to iterate over neighbors, where each iteration is a
     *    combination of edges in the two subtrees.  The first combination is
     *    skipped, since it would reverse the bisection.
     */

    /* Generate the edge lists. */
    tr_p_bedges_gen(a_tr, trt->bisect_edge);

    /* Calculate the offset of the neighbor from the beginning of this edge's
     * reconnection combination enumeration. */
    rem = a_neighbor - trt->offset;

    /* Avoid the first combination, since it would reverse the bisection. */
    rem++;

    *r_reconnect_a = a_tr->bedges[rem / a_tr->nbedges_b];
    *r_reconnect_b = a_tr->bedges[a_tr->nbedges_a + (rem % a_tr->nbedges_b)];
}

void *
tr_aux_get(cw_tr_t *a_tr)
{
    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);

    return a_tr->aux;
}

void
tr_aux_set(cw_tr_t *a_tr, void *a_aux)
{
    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);

    a_tr->aux = a_aux;
}

void
tr_mp_prepare(cw_tr_t *a_tr, cw_uint8_t *a_taxa[], cw_uint32_t a_ntaxa,
	      cw_uint32_t a_nchars)
{
    cw_uint32_t i;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    /* Prepare the tree. */
    tr_p_mp_prepare_recurse(a_tr, a_tr->base, CW_TR_NODE_NONE, a_taxa, a_ntaxa,
			    a_nchars);

    /* Prepare the temporary node. */
    tr_p_mp_prepare_recurse(a_tr, 0, CW_TR_NODE_NONE, a_taxa, a_ntaxa,
			    a_nchars);

    /* Initialize ps's for tre's. */
    for (i = 0; i < a_tr->nedges; i++)
    {
	if (a_tr->tres[i].ps == NULL)
	{
	    a_tr->tres[i].ps = tr_p_ps_new(a_tr);
	}
	tr_p_ps_prepare(a_tr, a_tr->tres[i].ps, a_nchars);
    }
}

void
tr_mp_finish(cw_tr_t *a_tr)
{
    cw_uint32_t i;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    /* Clean up the tree. */
    tr_p_mp_finish_recurse(a_tr, a_tr->base, CW_TR_NODE_NONE);

    /* Clean up the temporary node. */
    tr_p_mp_finish_recurse(a_tr, 0, CW_TR_NODE_NONE);

    /* Clean up the tre's. */
    for (i = 0; i < a_tr->nedges; i++)
    {
	if (a_tr->tres[i].ps != NULL)
	{
	    tr_p_ps_delete(a_tr, a_tr->tres[i].ps);
	    a_tr->tres[i].ps = NULL;
	}
    }
}

cw_uint32_t
tr_mp_score(cw_tr_t *a_tr, cw_uint32_t a_maxscore)
{
    cw_dassert(tr_p_validate(a_tr));

    return tr_p_mp_score(a_tr, a_tr->trns[0].ps,
			 a_tr->trns[a_tr->base].neighbors[0], a_tr->base,
			 CW_TR_NODE_NONE, a_maxscore);
}

void
tr_tbr_best_neighbors_mp(cw_tr_t *a_tr, cw_uint32_t a_max_hold)
{
    cw_dassert(tr_p_validate(a_tr));

    tr_p_tbr_neighbors_mp(a_tr, a_max_hold,
			  tr_mp_score(a_tr, CW_TR_MAXSCORE_NONE),
			  TR_HOLD_BEST);
}

void
tr_tbr_better_neighbors_mp(cw_tr_t *a_tr, cw_uint32_t a_max_hold)
{
    cw_dassert(tr_p_validate(a_tr));

    tr_p_tbr_neighbors_mp(a_tr, a_max_hold,
			  tr_mp_score(a_tr, CW_TR_MAXSCORE_NONE),
			  TR_HOLD_BETTER);
}

void
tr_tbr_all_neighbors_mp(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    tr_p_tbr_neighbors_mp(a_tr, CW_TR_HOLD_ALL,
			  tr_mp_score(a_tr, CW_TR_MAXSCORE_NONE),
			  TR_HOLD_ALL);
}

void
tr_held_finish(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    if (a_tr->held != NULL)
    {
	cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			  mema_arg_get(a_tr->mema), a_tr->held,
			  sizeof(cw_trh_t) * a_tr->heldlen);
	a_tr->held = NULL;
	a_tr->heldlen = 0;
	a_tr->nheld = 0;
    }
}

cw_uint32_t
tr_nheld_get(cw_tr_t *a_tr)
{
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->nheld;
}

void
tr_held_get(cw_tr_t *a_tr, cw_uint32_t a_held, cw_uint32_t *r_neighbor,
	    cw_uint32_t *r_score)
{
    cw_trh_t *trh;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));
    cw_check_ptr(a_tr->held);
    cw_assert(a_held < a_tr->nheld);

    trh = &a_tr->held[a_held];
    *r_neighbor = trh->neighbor;
    *r_score = trh->score;
}
