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
typedef struct cw_trr_s cw_trr_t;
typedef struct cw_tre_s cw_tre_t;
typedef struct cw_trt_s cw_trt_t;
typedef struct cw_trh_s cw_trh_t;

typedef uint32_t cw_tr_ring_t;

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
    uint32_t subtrees_score;
    uint32_t node_score;

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
    uint32_t nchars;
    cw_trc_t *achars;
};

/* Tree node for an unrooted bifurcating phylogenetic tree. */
struct cw_trn_s
{
#ifdef CW_DBG
    uint32_t magic;
#define CW_TRN_MAGIC 0x63329478
#endif

    union
    {
	/* Auxiliary opaque data pointer.  This is used by the treenode wrapper
	 * code for reference iteration. */
	void *aux;

	/* Spares linkage. */
	cw_tr_node_t link;
    } u;

    /* If CW_TR_NODE_TAXON_NONE, then the node is not a leaf node. */
    uint32_t taxon_num;

    /* Ring of trr's, which are associated with tre's. */
    qli_head rings;
};

/* Tree node edge ring element. */
struct cw_trr_s
{
    /* Ring linkage. */
    qri link;

    /* Node whose ring this trr is a part of. */
    cw_tr_node_t node;

    /* Used for Fitch parsimony scoring. */
    cw_tr_ps_t *ps;
};

/* Tree edge information. */
struct cw_tre_s
{
#ifdef CW_DBG
    uint32_t magic;
#define CW_TRE_MAGIC 0xa683fa07
#endif

    union
    {
	/* Auxiliary opaque data pointer.  This is used by the treeedge wrapper
	 * code for reference iteration. */
	void *aux;

	/* Spares linkage. */
	cw_tr_edge_t link;
    } u;

    /* Index into tr->trrs of the first edge ring element associated with this
     * edge.  The second edge ring element is directly after the first in
     * tr->trrs. */
    cw_tr_ring_t rings;

    /* Edge length. */
    double length;

    /* Used for Fitch parsimony scoring. */
    cw_tr_ps_t *ps;
};

/* TBR neighbor. */
struct cw_trt_s
{
    /* Number of neighbors that can be reached by doing TBR at edges before this
     * one.  This is also the neighbor index of the first neighbor that can be
     * reached by doing TBR on this edge. */
    uint32_t offset;

    /* Bisection edge. */
    uint32_t bisect_edge;
};

/* Held neighbor tree. */
struct cw_trh_s
{
    /* Neighbor index for the tree.  This can be passed to tr_tbr_neighbor_get()
     * to get the associated TBR parameters. */
    uint32_t neighbor;

    /* Fitch parsimony score for the neighboring tree. */
    uint32_t score;
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
    uint32_t magic;
#define CW_TR_MAGIC 0x39886394
#endif

    /* Used for memory allocation. */
    cw_mema_t *mema;

    /* Auxiliary opaque data pointer.  This is used by the treenode wrapper code
     * for reference iteration. */
    void *aux;

    /* True if this tree has been modified since the internal state (ntaxa,
     * nedges, bedges, trt, trei) was updated, false otherwise. */
    bool modified;

    /* Base of the tree (may or may not be set). */
    cw_tr_node_t base;

    /* Number of taxa in tree. */
    uint32_t ntaxa;

    /* Number of edges in tree.  This can be derived from ntaxa, but is used
     * often enough to make storing it worthwhile. */
    uint32_t nedges;

    /* bedges is an array of edge indices that is used for enumerating the edges
     * on each side of a logical tree bisection (used by TBR-related functions).
     * The first list starts at offset 0 and has nbedges_a elements.  The second
     * list starts at offset nbedges_a and has nbedges_b elements. */
    uint32_t *bedges;
    uint32_t nbedges_a;
    uint32_t nbedges_b;

    /* Array of triplets that store per-edge information that is used for
     * TBR-related functions.  There is one more element in trt than there are
     * edges in the tree.  This is critical to the way binary searching on the
     * array is done, and it also makes it easy to get the total number of
     * TBR neighbors this tree has (trt[nedges].offset).
     *
     * Only the first trtused elements are valid, since not all bisection edges
     * necessarily result in neighbors.
     *
     * trei is an array of edges (indices into the tres array) that correspond
     * to the entries in trt. */
    cw_trt_t *trt;
    uint32_t trtused;
    cw_tr_edge_t *trei;

    /* Pointer to an array of trn's.  ntrns is the total number of trn's, not
     * all of which are necessarily in use.
     *
     * The first element in trns is reserved as a temporary, which is used
     * whenever a single temporary trn is briefly needed.
     *
     * sparetrns is the index of the first spare trn in the spares stack. */
    cw_trn_t *trns;
    uint32_t ntrns;
    cw_tr_node_t sparetrns;

    /* tres is a pointer to an array of tre's.  ntres is the total number of
     * tre's, not all of which are necessarily in use.
     *
     * sparetres is the index of the first spare tre in the spares stack.
     *
     * trrs is a pointer to an array of trr's.  There are always twice as many
     * trr's in trrs as there are tre's in tres, and each pair of trr's in trrs
     * is implicitly associated with the corresponding tre in tres.
     *
     *   tres[0] <==> trrs[0], trrs[1]
     *   tres[1] <==> trrs[2], trrs[3]
     *   etc.
     */
    cw_tre_t *tres;
    uint32_t ntres;
    cw_tr_edge_t sparetres;
    cw_trr_t *trrs;

    /* held is an array of held neighbors.  The array is iteratively doubled as
     * necessary.  heldlen is the actual length of the array, and nheld is the
     * number of elements in use. */
    cw_trh_t *held;
    uint32_t heldlen;
    uint32_t nheld;
};
/******************************************************************************/

/* tr_ring. */

CW_P_INLINE void
tr_p_ring_init(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    cw_trr_t *trr;

    trr = &a_tr->trrs[a_ring];

    qri_new(a_tr->trrs, a_ring, link);
    trr->node = CW_TR_NODE_NONE;
    trr->ps = NULL;
}

CW_P_INLINE cw_tr_edge_t
tr_p_ring_edge_get(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    return (a_ring >> 1);
}

CW_P_INLINE cw_tr_node_t
tr_p_ring_node_get(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    return a_tr->trrs[a_ring].node;
}

CW_P_INLINE cw_tr_ring_t
tr_p_ring_other_get(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    return (a_ring ^ 1);
}

/******************************************************************************/

/* Validation functions. */

#ifdef CW_DBG
static bool
tr_p_validate(cw_tr_t *a_tr);

/* Validate an edge. */
static bool
tr_p_edge_validate(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_error("XXX Not implemented");

    return true;
}

/* Validate a node. */
static bool
tr_p_node_validate(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_trn_t *trn;
    cw_tr_ring_t ring;
    uint32_t nneighbors;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
    cw_assert(a_node < a_tr->ntrns);

    trn = &a_tr->trns[a_node];

    cw_assert(trn->magic == CW_TRN_MAGIC);

    nneighbors = 0;
    qli_foreach(ring, &trn->rings, a_tr->trrs, link)
    {
	/* Validate edge. */
	tr_p_edge_validate(a_tr, tr_p_ring_edge_get(a_tr, ring));

	nneighbors++;
    }

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	/* Only leaf nodes can have taxon numbers.  Leaf nodes have at most
	 * 1 neighbor. */
	cw_assert(nneighbors <= 1);
    }

    return true;
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
			  a_ps->achars, sizeof(cw_trc_t) * (a_ps->nchars + 8));
    }

    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
		      mema_arg_get(a_tr->mema),
		      a_ps, sizeof(cw_tr_ps_t));
}

CW_P_INLINE void
tr_p_ps_prepare(cw_tr_t *a_tr, cw_tr_ps_t *a_ps, uint32_t a_nchars)
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

/* tr_edge. */

CW_P_INLINE void
tr_p_edge_init(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_tre_t *tre;

    tre = &a_tr->tres[a_edge];

    tre->u.aux = NULL;
    tr_p_ring_init(a_tr, (a_edge << 1));
    tr_p_ring_init(a_tr, (a_edge << 1) + 1);
    tre->length = 0.0;
    tre->ps = NULL;

#ifdef CW_DBG
    tre->magic = CW_TRE_MAGIC;
#endif
}

CW_P_INLINE cw_tr_edge_t
tr_p_edge_alloc(cw_tr_t *a_tr)
{
    cw_tr_edge_t retval;

    if (a_tr->sparetres == CW_TR_EDGE_NONE)
    {
	uint32_t i, nspares;

	if (a_tr->tres == NULL)
	{
	    a_tr->tres
		= (cw_tre_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					       mema_arg_get(a_tr->mema),
					       sizeof(cw_tre_t));
	    cw_assert(a_tr->trrs == NULL);
	    a_tr->trrs
		= (cw_trr_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					       mema_arg_get(a_tr->mema),
					       sizeof(cw_trr_t) * 2);
	    nspares = 1;
	    a_tr->ntres = 1;
	}
	else
	{
	    a_tr->tres
		= (cw_tre_t *) cw_opaque_realloc(mema_realloc_get(a_tr->mema),
						 mema_arg_get(a_tr->mema),
						 a_tr->tres,
						 sizeof(cw_tre_t)
						 * a_tr->ntres * 2,
						 sizeof(cw_tre_t)
						 * a_tr->ntres);
	    cw_check_ptr(a_tr->trrs);
	    a_tr->trrs
		= (cw_trr_t *) cw_opaque_realloc(mema_realloc_get(a_tr->mema),
						 mema_arg_get(a_tr->mema),
						 a_tr->trrs,
						 sizeof(cw_trr_t)
						 * a_tr->ntres * 4,
						 sizeof(cw_trr_t)
						 * a_tr->ntres * 2);
	    nspares = a_tr->ntres;
	    a_tr->ntres *= 2;
	}

	/* Initialize last spare. */
	a_tr->sparetres = a_tr->ntres - 1;
	a_tr->tres[a_tr->sparetres].u.link = CW_TR_EDGE_NONE;

	/* Insert other spares into spares stack. */
	for (i = 1; i < nspares; i++)
	{
	    a_tr->sparetres--;
	    a_tr->tres[a_tr->sparetres].u.link = a_tr->sparetres + 1;
	}
    }

    /* Remove a spare from the spares stack. */
    retval = a_tr->sparetres;
    a_tr->sparetres = a_tr->tres[retval].u.link;

    /* Initialize retval. */
    tr_p_edge_init(a_tr, retval);

    return retval;
}

CW_P_INLINE void
tr_p_edge_dealloc(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_tre_t *tre;

    tre = &a_tr->tres[a_edge];

    if (tre->ps != NULL)
    {
	tr_p_ps_delete(a_tr, tre->ps);
    }

#ifdef CW_DBG
    memset(&a_tr->tres[a_edge], 0x5a, sizeof(cw_tre_t));
    memset(&a_tr->trrs[a_edge << 1], 0x5a, sizeof(cw_trr_t) * 2);
#endif

    a_tr->tres[a_edge].u.link = a_tr->sparetres;
    a_tr->sparetres = a_edge;
}

cw_tr_edge_t
tr_edge_new(cw_tr_t *a_tr)
{
    return tr_p_edge_alloc(a_tr);
}

void
tr_edge_delete(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_dassert(tr_p_edge_validate(a_tr, a_edge));

    tr_p_edge_dealloc(a_tr, a_edge);
}

CW_P_INLINE cw_tr_ring_t
tr_p_edge_ring_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_end)
{
    return ((a_edge << 1) + a_end);
}

cw_tr_node_t
tr_edge_node_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_i)
{
    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_assert(a_i == 0 || a_i == 1);

    return a_tr->trrs[(a_edge << 1) + a_i].node;
}

void
tr_edge_next_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_i,
		 cw_tr_edge_t *r_next, uint32_t *r_i)
{
    cw_tr_ring_t ringind;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_assert(a_i == 0 || a_i == 1);

    ringind = qri_next(a_tr->trrs, (a_edge << 1) + a_i, link);
    *r_next = a_tr->trrs[ringind].node;
    *r_i = (ringind & 1);
}

void
tr_edge_prev_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_i,
		 cw_tr_edge_t *r_prev, uint32_t *r_i)
{
    cw_tr_ring_t ringind;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_assert(a_i == 0 || a_i == 1);

    ringind = qri_prev(a_tr->trrs, (a_edge << 1) + a_i, link);
    *r_prev = a_tr->trrs[ringind].node;
    *r_i = (ringind & 1);
}

double
tr_edge_length_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_dassert(tr_p_edge_validate(a_tr, a_edge));

    return a_tr->tres[a_edge].length;
}

void
tr_edge_length_set(cw_tr_t *a_tr, cw_tr_edge_t a_edge, double a_length)
{
    cw_dassert(tr_p_edge_validate(a_tr, a_edge));

    a_tr->tres[a_edge].length = a_length;
}

void *
tr_edge_aux_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_dassert(tr_p_edge_validate(a_tr, a_edge));

    return a_tr->tres[a_edge].u.aux;
}

void
tr_edge_aux_set(cw_tr_t *a_tr, cw_tr_edge_t a_edge, void *a_aux)
{
    cw_dassert(tr_p_edge_validate(a_tr, a_edge));

    a_tr->tres[a_edge].u.aux = a_aux;
}

#if (0) // XXX
void
tr_node_join(cw_tr_t *a_tr, cw_tr_node_t a_a, cw_tr_node_t a_b,
	     cw_tr_edge_t a_edge)
{
    cw_trn_t *trn_a, *trn_b;

    cw_dassert(tr_p_node_validate(a_tr, a_a));
    cw_dassert(tr_p_node_validate(a_tr, a_b));
    cw_assert(a_a != a_b);
    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 0) == CW_TR_NODE_NONE);
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 1) == CW_TR_NODE_NONE);

    trn_a = &a_tr->trns[a_a];
    trn_b = &a_tr->trns[a_b];

#ifdef CW_DBG
    /* Make sure that the nodes aren't already connected. */
    {
	cw_tr_ring_t ring;

	qli_foreach(ring, &trn_a->rings, a_tr->trrs, link)
	{
	    cw_assert(tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring))
		      != a_b);
	}
	qli_foreach(ring, &trn_b->rings, a_tr->trrs, link)
	{
	    cw_assert(tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring))
		      != a_a);
	}
    }
#endif

    /* Attach nodes to edge. */
    qli_tail_insert(&trn_a->rings, a_tr->trrs,
		    tr_p_edge_ring_get(a_tr, a_edge, 0), link);
    qli_tail_insert(&trn_b->rings, a_tr->trrs,
		    tr_p_edge_ring_get(a_tr, a_edge, 1), link);

    a_tr->modified = true;

    cw_dassert(tr_p_node_validate(a_tr, a_a));
    cw_dassert(tr_p_node_validate(a_tr, a_b));
}

cw_tr_edge_t
tr_node_detach(cw_tr_t *a_tr, cw_tr_node_t a_a, cw_tr_node_t a_b)
{
    cw_trn_t *trn_a, *trn_b;
    uint32_t i, j;

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

    a_tr->modified = true;

    cw_dassert(tr_p_node_validate(a_tr, a_a));
    cw_dassert(tr_p_node_validate(a_tr, a_b));
}
#endif // XXX

void
tr_edge_attach(cw_tr_t *a_tr, cw_tr_edge_t a_edge, cw_tr_node_t a_node_a,
	       cw_tr_node_t a_node_b)
{
    cw_error("XXX Not implemented");
}

void
tr_edge_detach(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_error("XXX Not implemented");
}

/******************************************************************************/

/* tr_node. */

CW_P_INLINE void
tr_p_node_init(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_trn_t *trn;

    trn = &a_tr->trns[a_node];

    trn->u.aux = NULL;
    trn->taxon_num = CW_TR_NODE_TAXON_NONE;
    qli_new(&trn->rings);

#ifdef CW_DBG
    trn->magic = CW_TRN_MAGIC;
#endif
}

CW_P_INLINE cw_tr_node_t
tr_p_node_alloc(cw_tr_t *a_tr)
{
    cw_tr_node_t retval;

    if (a_tr->sparetrns == CW_TR_NODE_NONE)
    {
	uint32_t i, nspares;

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
	a_tr->sparetrns = a_tr->ntrns - 1;
	a_tr->trns[a_tr->sparetrns].u.link = CW_TR_NODE_NONE;

	/* Insert other spares into spares stack. */
	for (i = 1; i < nspares; i++)
	{
	    a_tr->sparetrns--;
	    a_tr->trns[a_tr->sparetrns].u.link = a_tr->sparetrns + 1;
	}
    }

    /* Remove a spare from the spares stack. */
    retval = a_tr->sparetrns;
    a_tr->sparetrns = a_tr->trns[retval].u.link;

    /* Initialize retval. */
    tr_p_node_init(a_tr, retval);

    return retval;
}

CW_P_INLINE void
tr_p_node_dealloc(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_trn_t *trn;

    trn = &a_tr->trns[a_node];
    
#ifdef CW_DBG
    memset(&a_tr->trns[a_node], 0x5a, sizeof(cw_trn_t));
#endif

    a_tr->trns[a_node].u.link = a_tr->sparetrns;
    a_tr->sparetrns = a_node;
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
    // XXX Assert no edges?
    cw_assert(qli_first(&a_tr->trns[a_node].rings) == UINT_MAX);

    tr_p_node_dealloc(a_tr, a_node);
}

uint32_t
tr_node_taxon_num_get(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].taxon_num;
}

void
tr_node_taxon_num_set(cw_tr_t *a_tr, cw_tr_node_t a_node,
		      uint32_t a_taxon_num)
{
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    a_tr->trns[a_node].taxon_num = a_taxon_num;

    a_tr->modified = true;
}

void
tr_node_edge_get(cw_tr_t *a_tr, cw_tr_node_t a_node,
		 cw_tr_edge_t *r_edge, uint32_t *r_i)
{
    cw_trn_t *trn;
    cw_tr_ring_t ringind;

    cw_dassert(tr_p_node_validate(a_tr, a_node));

    trn = &a_tr->trns[a_node];

    ringind = qli_first(&trn->rings);
    if (ringind != CW_TR_EDGE_NONE)
    {
	*r_edge = (ringind >> 1);
	*r_i = (ringind & 1);
    }
    else
    {
	*r_edge = CW_TR_EDGE_NONE;
    }
}

void *
tr_node_aux_get(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].u.aux;
}

void
tr_node_aux_set(cw_tr_t *a_tr, cw_tr_node_t a_node, void *a_aux)
{
    cw_dassert(tr_p_node_validate(a_tr, a_node));

    a_tr->trns[a_node].u.aux = a_aux;
}

/******************************************************************************/

/* tr. */

/* Initialize everything except trns and sparetrns. */
CW_P_INLINE void
tr_p_new(cw_tr_t *a_tr, cw_mema_t *a_mema)
{
    a_tr->mema = a_mema;
    a_tr->aux = NULL;
    a_tr->modified = false;
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

/* Recursively traverse the tree, count the number of taxa, and find the lowest
 * numbered taxon. */
static cw_tr_node_t
tr_p_update_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
		    uint32_t *r_ntaxa, cw_tr_node_t a_root)
{
    cw_tr_node_t retval, root, troot;
    cw_tr_ring_t ring;
    cw_trn_t *trn;

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
    // XXX Get ring object that links edge to previous node.  Need to change
    // function args to make this happen.
    qri_others_foreach(ring, a_tr->trrs, XXX, link)
//    qli_foreach(ring, &trn->rings, a_tr->trrs, link)
    {
	if (tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring)) != a_prev)
	{
	    troot = tr_p_update_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
					ring, r_ntaxa, root);
//	    troot = tr_p_update_recurse(a_tr, tr_p_ring_node_get(a_tr, ring),
//					a_node, r_ntaxa, root);
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
static bool
tr_p_reachable(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
	       cw_tr_node_t a_other)
{
    uint32_t retval;
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    trn = &a_tr->trns[a_node];

    if (a_node == a_other)
    {
	retval = true;
	goto RETURN;
    }

    qli_foreach(ring, &trn->rings, a_tr->trrs, link)
    {
	if (tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring)) != a_prev)
	{
	    if ((retval = tr_p_reachable(a_tr, tr_p_ring_node_get(a_tr, ring),
					 a_node, a_other)))
	    {
		goto RETURN;
	    }
	}
    }

    retval = false;
    RETURN:
    return retval;
}

/* Return the number of taxa with number a_taxon_num in the subtree rooted at
 * a_node. */
static uint32_t
tr_p_validate_recurse(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev,
		      uint32_t a_taxon_num)
{
    uint32_t retval;
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    tr_p_node_validate(a_tr, a_node);

    trn = &a_tr->trns[a_node];

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	/* Leaf node. */
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

    qli_foreach(ring, &trn->rings, a_tr->trrs, link)
    {
	if (tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring)) != a_prev)
	{
	    retval += tr_p_validate_recurse(a_tr,
					    tr_p_ring_node_get(a_tr, ring),
					    a_node, a_taxon_num);
	}
    }

    return retval;
}

/* Validate a tree. */
static bool
tr_p_validate(cw_tr_t *a_tr)
{
    uint32_t i, ntaxa;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
    cw_assert(a_tr->modified == false);

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
	    cw_assert(a_tr->trns[i].u.link == CW_TR_NODE_NONE
		      || a_tr->trns[a_tr->trns[i].u.link].magic
		      != CW_TRN_MAGIC);
	}
    }

    cw_assert(a_tr->sparetrns == CW_TR_NODE_NONE
	      || a_tr->trns[a_tr->sparetrns].magic != CW_TRN_MAGIC);

    return true;
}
#endif

#if (0) // XXX
// XXX This function needs to be converted to use trei.
CW_P_INLINE void
tr_p_edge_get(cw_tr_t *a_tr, uint32_t a_edge, cw_tr_node_t *r_node_a,
	      cw_tr_node_t *r_node_b)
{
    *r_node_a = a_tr->tres[a_edge].node_a;
    *r_node_b = a_tr->tres[a_edge].node_b;
}
#endif // XXX

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
				uint32_t *r_edge_count,
				uint32_t *r_bisection_edge)
{
    uint32_t prev_edge_count;
    cw_trn_t *trn;
    cw_tr_ring_t ring;
    cw_tr_node_t node;

    cw_assert(a_node != CW_TR_NODE_NONE);

    /* Save the previous edge count, in case it ends up being the index of the
     * edge adjacent to the bisection. */
    prev_edge_count = *r_edge_count;

    trn = &a_tr->trns[a_node];

    qli_foreach(ring, &trn->rings, a_tr->trrs, link)
    {
	node = tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring));
	if (node == a_other)
	{
	    /* Store the index of the edge adjacent to the bisection. */
	    if (prev_edge_count > 0)
	    {
		*r_bisection_edge = prev_edge_count - 1;
	    }
	    else
	    {
		*r_bisection_edge = CW_TR_EDGE_NONE;
	    }
	}
	else if (node != a_prev)
	{
	    /* Increment edge count before recursing. */
	    (*r_edge_count)++;

	    /* Recurse into neighbor subtree. */
	    tr_p_bisection_edge_get_recurse(a_tr,
					    tr_p_ring_node_get(a_tr, ring),
					    a_other, a_node, r_edge_count,
					    r_bisection_edge);
	}
    }
}

static void
tr_p_ntaxa_nedges_update(cw_tr_t *a_tr)
{
    uint32_t ntaxa;

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
				     cw_tr_node_t a_prev, uint32_t *ar_edges,
				     uint32_t *ar_nedges)
{
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    trn = &a_tr->trns[a_node];

    qli_foreach(ring, &trn->rings, a_tr->trrs, link)
    {
	if (tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring)) != a_prev)
	{
	    /* Add edge to list. */
	    ar_edges[*ar_nedges] = tr_p_ring_edge_get(a_tr, ring);
	    (*ar_nedges)++;

	    /* Recurse into neighbor subtree. */
	    tr_p_bisection_edge_list_gen_recurse(a_tr,
						 tr_p_ring_other_get(a_tr,
								     ring),
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
			     cw_tr_node_t a_other, uint32_t *ar_edges,
			     uint32_t *ar_nedges)
{
    cw_trn_t *trn;
    cw_tr_node_t a, b;
    uint32_t i, i_a;
#ifdef CW_DBG
    uint32_t i_b;
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
tr_p_bedges_gen(cw_tr_t *a_tr, uint32_t a_bisect)
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
tr_p_trt_update(cw_tr_t *a_tr, uint32_t a_nedges_prev)
{
    uint32_t i, j, n, offset;

    cw_assert(a_tr->modified == false);

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
			uint32_t *ar_edge_count)
{
    cw_trn_t *trn, *ttrn;
    uint32_t i, j;

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
tr_p_tre_update(cw_tr_t *a_tr, uint32_t a_nedges_prev)
{
    uint32_t edge_count;

    cw_assert(a_tr->modified == false);

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
	    uint32_t i;

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
	uint32_t i;

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
tr_p_bedges_update(cw_tr_t *a_tr, uint32_t a_nedges_prev)
{
    /* Allocate/reallocate/deallocate bedges.  To keep things simple, allocate
     * as big an array as there are edges, even though not quite that many are
     * ever used. */
    if (a_tr->bedges == NULL)
    {
	/* Allocate bedges. */
	a_tr->bedges
	    = (uint32_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					      mema_arg_get(a_tr->mema),
					      sizeof(uint32_t)
					      * a_tr->nedges);
    }
    else if (a_tr->nedges != a_nedges_prev)
    {
	if (a_tr->nedges > 0)
	{
	    /* Reallocate bedges. */
	    a_tr->bedges = (uint32_t *)
		cw_opaque_realloc(mema_realloc_get(a_tr->mema),
				  mema_arg_get(a_tr->mema),
				  a_tr->bedges,
				  sizeof(uint32_t) * a_tr->nedges,
				  sizeof(uint32_t) * a_nedges_prev);
	}
	else
	{
	    /* Deallocate bedges. */
	    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			      mema_arg_get(a_tr->mema),
			      a_tr->bedges,
			      sizeof(uint32_t) * a_nedges_prev);
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
	uint32_t nedges_prev;

	/* Store nedges before updating. */
	nedges_prev = a_tr->nedges;

	/* Update ntaxa and nedges. */
	tr_p_ntaxa_nedges_update(a_tr);

	/* Reset the modified flag. */
	a_tr->modified = false;

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
static uint32_t
tr_p_canonize(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_prev)
{
    uint32_t retval;
    uint32_t i, j, t;
    uint32_t subtree_mins[CW_TR_NODE_MAX_NEIGHBORS - 1];
    uint32_t subtree_inds[CW_TR_NODE_MAX_NEIGHBORS - 1];
    bool swapped;
    cw_trn_t *trn;
#ifdef CW_DBG
    uint32_t nneighbors = 0;
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
	swapped = false;

	for (i = 0; i + 1 < j; i++)
	{
	    if (subtree_mins[i] > subtree_mins[i + 1])
	    {
		swapped = true;

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
    uint32_t i;
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
tr_p_bisect(cw_tr_t *a_tr, uint32_t a_edge,
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
	uint32_t i;
	bool connected;

	for (i = 0, connected = false; i < CW_TR_NODE_MAX_NEIGHBORS; i++)
	{
	    if (a_tr->trns[*r_node_a].neighbors[i] == *r_node_b)
	    {
		connected = true;
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
CW_P_INLINE uint32_t
tr_p_reconnect_ready(cw_tr_t *a_tr, cw_tr_node_t a_node,
		     uint32_t a_reconnect_a, uint32_t a_reconnect_b)
{
    uint32_t retval;
    uint32_t i, nneighbors;
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
		       uint32_t a_reconnect, uint32_t a_ready)
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
			char *a_taxa[], uint32_t a_ntaxa,
			uint32_t a_nchars)
{
    cw_trn_t *trn;
    cw_tr_node_t node;
    uint32_t i, taxon_num;

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
	char *chars;

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
    uint32_t i;

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
CW_P_INLINE void
tr_p_mp_ia32_pscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
		    cw_tr_ps_t *a_b)
{
    /* Only calculate the parent's node score if the cached value is invalid. */
    if (a_a->parent != a_p || a_b->parent != a_p)
    {
	uint32_t endlimit, curlimit, i, nchars, ns;
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

	/* The inner loop can be run a maximum of 255 times before the partial
	 * node score results (stored in %xmm5) are added to ns (otherwise,
	 * overflow could occur).  Therefore, the outer loop calculates the
	 * upper bound for the inner loop, thereby avoiding extra computation in
	 * the inner loop. */
	endlimit = nchars ^ (nchars & 0xf);
	curlimit = 255 * 16;
	if (curlimit > endlimit)
	{
	    curlimit = endlimit;
	}
	for (;;)
	{
	    /* Use SSE2 to evaluate as many of the characters as possible.  This
	     * loop handles 16 characters per iteration. */
	    for (i = 0; i < curlimit; i += 16)
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

		    /* Create bitmasks according to whether the character state
		     * sets are empty.
		     *
		     * c = p ? 0x00 : 0xff;
		     * e = (c & d);
		     * s = p ? 0 : 1;
		     */
		    "pxor %%xmm2, %%xmm2;"
		    "pcmpeqb %%xmm1, %%xmm2;" /* xmm2 contains c. */
		    "pand %%xmm2, %%xmm0;" /* xmm0 contains e. */
		    "pand %%xmm7, %%xmm2;" /* xmm2 contains s. */

		    /* Update node score.  Each byte in %xmm5 is only capable of
		     * holding up to 255, which is why the outer loop is
		     * necessary.
		     *
		     * ns += s;
		     */
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

	    /* Update ns and reset pns. */
	    {
		uint32_t j;
		unsigned char pns[16];

		asm volatile (
		    "movdqu %%xmm5, %[pns];"
		    "pxor %%xmm5, %%xmm5;"
		    : [pns] "=m" (*pns)
		    :
		    : "memory"
		    );

		for (j = 0; j < 16; j++)
		{
		    ns += pns[j];
		}
	    }

	    /* Break out of the loop if the bound for the inner loop was the
	     * maximum possible. */
	    if (curlimit == endlimit)
	    {
		break;
	    }
	    /* Update the bound for the inner loop, taking care not to exceed
	     * the maximum possible bound. */
	    curlimit += 255 * 16;
	    if (curlimit > endlimit)
	    {
		curlimit = endlimit;
	    }
	}

	/* Evaluate the last 0-15 characters that weren't evaluated in the above
	 * loop. */
	{
	    uint32_t a, b, p, c, s;

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
}
#endif

static void
tr_p_mp_c_pscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
		 cw_tr_ps_t *a_b)
{
//#define CW_TR_MP_PSCORE_VALIDATE
#ifdef CW_TR_MP_PSCORE_VALIDATE
    bool cached;
#endif

    /* Only calculate the parent's node score if the cached value is invalid. */
    if (a_a->parent != a_p || a_b->parent != a_p)
#ifdef CW_TR_MP_PSCORE_VALIDATE
    {
	cached = false;
    }
    else
    {
	cached = true;
    }
#endif
    {
	uint32_t i, nchars, ns, a, b, p, c, s;
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
}

CW_P_INLINE void
tr_p_mp_pscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b)
{
#ifdef CW_CPU_IA32
    if (modcrux_ia32_use_sse2)
    {
	tr_p_mp_ia32_pscore(a_tr, a_p, a_a, a_b);
    }
    else
#endif
    {
	tr_p_mp_c_pscore(a_tr, a_p, a_a, a_b);
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
		      cw_trn_t *a_other)
{
    uint32_t i;
    cw_trn_t *trn, *a, *b, *o;
    cw_tr_node_t node;

    a = NULL;
    b = NULL;
    o = NULL;

    /* Set a, b, and o, according to which neighbors are subtrees or on the
     * other side of the bisection. */
    cw_assert(CW_TR_NODE_MAX_NEIGHBORS == 3);
#define TR_P_MP_SCORE_RECURSE_NEIGHBOR_ITERATE()			\
    node = a_node->neighbors[i];					\
    if (node != CW_TR_NODE_NONE)					\
    {									\
	trn = &a_tr->trns[node];					\
	if (trn == a_other)						\
	{								\
    	o = trn;							\
	}								\
	else if (trn != a_prev)						\
	{								\
	    if (a == NULL)						\
	    {								\
		a = trn;						\
	    }								\
	    else							\
	    {								\
		cw_assert(b == NULL);					\
		b = trn;						\
	    }								\
	}								\
    }

    i = 0;
    TR_P_MP_SCORE_RECURSE_NEIGHBOR_ITERATE();
    i++;
    TR_P_MP_SCORE_RECURSE_NEIGHBOR_ITERATE();
    i++;
    TR_P_MP_SCORE_RECURSE_NEIGHBOR_ITERATE();
#undef TR_P_MP_SCORE_RECURSE_NEIGHBOR_ITERATE

    /* Recursively calculate partial scores for the subtrees, then calculate the
     * partial score for this node. */
    if (b != NULL)
    {
	/* Recurse into subtrees. */
	tr_p_mp_score_recurse(a_tr, a, a_node, a_other);
	tr_p_mp_score_recurse(a_tr, b, a_node, a_other);

	/* Clear cached other value. */
	a_node->ps->other = o;

	/* Calculate this node's partial score. */

	/* Calculate the partial score for this node. */
	cw_check_ptr(a);
	tr_p_mp_pscore(a_tr, a_node->ps, a->ps, b->ps);
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
	tr_p_mp_score_recurse(a_tr, a, a_node, a_other);

	/* Copy a's scores to this node, rather than calculating a partial
	 * score.  This node is merely a filler node, as far as scoring is
	 * concerned. */
	if (a != NULL)
	{
	    tr_p_mp_passpscore(a_tr, a_node->ps, a->ps);
	}
    }
}

static uint32_t
tr_p_mp_score(cw_tr_t *a_tr, cw_tr_ps_t *a_ps, cw_tr_node_t a_node_a,
	      cw_tr_node_t a_node_b, cw_tr_node_t a_other)
{
    uint32_t retval;
    bool maxed;

    maxed = false;
    tr_p_mp_score_recurse(a_tr, &a_tr->trns[a_node_a], &a_tr->trns[a_node_b],
			  &a_tr->trns[a_other]);
    if (maxed)
    {
	retval = CW_TR_MAXSCORE_NONE;
	goto RETURN;
    }

    tr_p_mp_score_recurse(a_tr, &a_tr->trns[a_node_b], &a_tr->trns[a_node_a],
			  &a_tr->trns[a_other]);
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
tr_p_bisection_edge_list_mp(cw_tr_t *a_tr, uint32_t *a_edges,
			    uint32_t a_nedges, cw_tr_node_t a_other)
{
    uint32_t i;
    cw_tre_t *tre;

    for (i = 0; i < a_nedges; i++)
    {
	if (a_edges[i] != CW_TR_NODE_EDGE_NONE)
	{
	    tre = &a_tr->tres[a_edges[i]];

	    tr_p_mp_score(a_tr, tre->ps, tre->node_a, tre->node_b, a_other);
	}
    }
}

/* Hold a tree.  If a_max_held is exceeded, the tree is not held.  This
 * introduces a bias in which trees are held.  There exist algorithms for making
 * this an unbiased process, but there is no need for that functionality at the
 * moment. */
CW_P_INLINE void
tr_p_hold(cw_tr_t *a_tr, uint32_t a_max_hold, uint32_t a_neighbor,
	  uint32_t a_score)
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
tr_p_tbr_neighbors_mp(cw_tr_t *a_tr, uint32_t a_max_hold,
		      uint32_t a_maxscore, cw_tr_hold_how_t a_how)
{
    uint32_t neighbor, i, j, k, edge_a, edge_b;
    uint32_t score;
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
				    tre->node_b);
	tr_p_bisection_edge_list_mp(a_tr, &a_tr->bedges[a_tr->nbedges_a],
				    a_tr->nbedges_b, tre->node_a);

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
		 * reversing the bisection. */
		if (j == 0 && k == 0)
		{
		    continue;
		}

		/* Calculate the final parsimony score for this reconnection.
		 * Clear the parent pointers of ps_[ab], to make sure that the
		 * score is actually calculated. */
		ps_a->parent = NULL;
		ps_b->parent = NULL;
		tr_p_mp_pscore(a_tr, ps, ps_a, ps_b);
		score = (ps->subtrees_score + ps->node_score);

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

    retval->sparetrns = CW_TR_NODE_NONE;

    return retval;
}

cw_tr_t *
tr_dup(cw_tr_t *a_tr)
{
    cw_tr_t *retval;
    cw_opaque_alloc_t *alloc;
    void *arg;
    uint32_t i;

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
    retval->sparetrns = a_tr->sparetrns;

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
			  a_tr->bedges, sizeof(uint32_t) * a_tr->nedges);
    }

    if (a_tr->tres != NULL)
    {
	uint32_t i;

	for (i = 0; i < a_tr->nedges; i++)
	{
	    if (a_tr->tres[i].ps != NULL)
	    {
		tr_p_ps_delete(a_tr, a_tr->tres[i].ps);
	    }
	}

	cw_opaque_dealloc(dealloc, arg, a_tr->tres,
			  sizeof(cw_tre_t) * a_tr->nedges);
    }

    cw_opaque_dealloc(dealloc, arg, a_tr, sizeof(cw_tr_t));
}

uint32_t
tr_ntaxa_get(cw_tr_t *a_tr)
{
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->ntaxa;
}

uint32_t
tr_nedges_get(cw_tr_t *a_tr)
{
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->nedges;
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

    a_tr->modified = true;
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
	a_tr->modified = false;
    }

    if (a_tr->base != CW_TR_NODE_NONE)
    {
	uint32_t ntaxa;

	/* Set base to be the lowest-numbered taxon. */
	ntaxa = 0;
	a_tr->base = tr_p_update_recurse(a_tr, a_tr->base, CW_TR_NODE_NONE,
					 &ntaxa, CW_TR_NODE_NONE);

	/* Canonize the tree. */
	tr_p_canonize(a_tr, a_tr->base, CW_TR_NODE_NONE);
    }

    /* Reset the modified flag. */
    a_tr->modified = false;

    /* Now update tre and trt. */
    tr_p_tre_update(a_tr, a_tr->nedges);
    tr_p_trt_update(a_tr, a_tr->nedges);

    cw_dassert(tr_p_validate(a_tr));
}

void
tr_tbr(cw_tr_t *a_tr, uint32_t a_bisect, uint32_t a_reconnect_a,
       uint32_t a_reconnect_b)
{
    cw_tr_node_t node_a, node_b;
    uint32_t ready_a, ready_b;

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
	uint32_t treconnect;

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
    a_tr->modified = false;

    /* Update tre and trt. */
    tr_p_tre_update(a_tr, a_tr->nedges);
    tr_p_trt_update(a_tr, a_tr->nedges);

    cw_dassert(tr_p_validate(a_tr));
}

uint32_t
tr_tbr_nneighbors_get(cw_tr_t *a_tr)
{
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->trt[a_tr->trtused].offset;
}

void
tr_tbr_neighbor_get(cw_tr_t *a_tr, uint32_t a_neighbor,
		    uint32_t *r_bisect, uint32_t *r_reconnect_a,
		    uint32_t *r_reconnect_b)
{
    cw_trt_t key, *trt;
    uint32_t rem;

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
tr_mp_prepare(cw_tr_t *a_tr, char *a_taxa[], uint32_t a_ntaxa,
	      uint32_t a_nchars)
{
    uint32_t i;

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
    uint32_t i;

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

uint32_t
tr_mp_score(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    return tr_p_mp_score(a_tr, a_tr->trns[0].ps,
			 a_tr->trns[a_tr->base].neighbors[0], a_tr->base,
			 CW_TR_NODE_NONE);
}

void
tr_tbr_best_neighbors_mp(cw_tr_t *a_tr, uint32_t a_max_hold)
{
    cw_dassert(tr_p_validate(a_tr));

    tr_p_tbr_neighbors_mp(a_tr, a_max_hold,
			  CW_TR_MAXSCORE_NONE,
			  TR_HOLD_BEST);
}

void
tr_tbr_better_neighbors_mp(cw_tr_t *a_tr, uint32_t a_max_hold)
{
    cw_dassert(tr_p_validate(a_tr));

    tr_p_tbr_neighbors_mp(a_tr, a_max_hold,
			  tr_mp_score(a_tr),
			  TR_HOLD_BETTER);
}

void
tr_tbr_all_neighbors_mp(cw_tr_t *a_tr)
{
    cw_dassert(tr_p_validate(a_tr));

    tr_p_tbr_neighbors_mp(a_tr, CW_TR_HOLD_ALL,
			  CW_TR_MAXSCORE_NONE,
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

uint32_t
tr_nheld_get(cw_tr_t *a_tr)
{
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    return a_tr->nheld;
}

void
tr_held_get(cw_tr_t *a_tr, uint32_t a_held, uint32_t *r_neighbor,
	    uint32_t *r_score)
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
