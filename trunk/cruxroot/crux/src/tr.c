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

// XXX Remove spare node (trns[0]); it isn't needed at all.
#include "../include/modcrux.h"

typedef struct cw_tr_ps_s cw_tr_ps_t;
typedef struct cw_trn_s cw_trn_t;
typedef struct cw_trr_s cw_trr_t;
typedef struct cw_tre_s cw_tre_t;
typedef struct cw_trt_s cw_trt_t;
typedef struct cw_trh_s cw_trh_t;

typedef uint32_t cw_tr_ring_t;
#define CW_TR_RING_NONE UINT_MAX

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
     * records the bisection edge.  This is used when deciding whether to push
     * the value of parent down to the child when recursively scoring. */
    cw_tr_edge_t bisect;

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

    // XXX Remove once MP scoring uses multiple views.
    cw_tr_ps_t *ps;
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
     * nedges, bedges, trt, trti) was updated, false otherwise. */
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
     * trti is an array of edges (indices into the tres array) that correspond
     * to the entries in trt. */
    cw_trt_t *trt;
    uint32_t trtused;
    cw_tr_edge_t *trti;

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

/* Validate a ring. */
static bool
tr_p_ring_validate(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    cw_trr_t *trr;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
    cw_assert(a_ring < a_tr->ntres);

    trr = &a_tr->trrs[a_ring];

    cw_assert(trr->node < a_tr->ntrns);

    return true;
}

/* Validate an edge. */
static bool
tr_p_edge_validate(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_tre_t *tre;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
    cw_assert(a_edge < a_tr->ntres);

    tre = &a_tr->tres[a_edge];

    cw_assert(tre->magic == CW_TRE_MAGIC);

    tr_p_ring_validate(a_tr, (a_edge << 1));
    tr_p_ring_validate(a_tr, (a_edge << 1) + 1);

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

void
tr_edge_attach(cw_tr_t *a_tr, cw_tr_edge_t a_edge, cw_tr_node_t a_node_a,
	       cw_tr_node_t a_node_b)
{
    cw_trn_t *trn_a, *trn_b;
    cw_tr_ring_t ring;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 0) == CW_TR_NODE_NONE);
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 1) == CW_TR_NODE_NONE);
    cw_dassert(tr_p_node_validate(a_tr, a_node_a));
    cw_dassert(tr_p_node_validate(a_tr, a_node_b));
    cw_assert(a_node_a != a_node_b);

    trn_a = &a_tr->trns[a_node_a];
    trn_b = &a_tr->trns[a_node_b];

#ifdef CW_DBG
    /* Make sure that the nodes aren't already connected. */
    qli_foreach(ring, &trn_a->rings, a_tr->trrs, link)
    {
	cw_assert(tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring))
		  != a_node_b);
    }
    qli_foreach(ring, &trn_b->rings, a_tr->trrs, link)
    {
	cw_assert(tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring))
		  != a_node_a);
    }
#endif

    /* First end. */
    ring = tr_p_edge_ring_get(a_tr, a_edge, 0);
    qli_tail_insert(&trn_a->rings, a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = a_node_a;

    /* Second end. */
    ring = tr_p_edge_ring_get(a_tr, a_edge, 1);
    qli_tail_insert(&trn_b->rings, a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = a_node_b;

    /* Mark tree as modified. */
    a_tr->modified = true;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_dassert(tr_p_node_validate(a_tr, a_node_a));
    cw_dassert(tr_p_node_validate(a_tr, a_node_b));
}

void
tr_edge_detach(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_tr_ring_t ring;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 0) != CW_TR_NODE_NONE);
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 1) != CW_TR_NODE_NONE);

    ring = tr_p_edge_ring_get(a_tr, a_edge, 0);
    qri_remove(a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = CW_TR_NODE_NONE;

    ring = tr_p_edge_ring_get(a_tr, a_edge, 1);
    qri_remove(a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = CW_TR_NODE_NONE;

    /* Mark tree as modified. */
    a_tr->modified = true;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
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
    trn->ps = NULL;

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
    cw_assert(qli_first(&a_tr->trns[a_node].rings) == CW_TR_RING_NONE);

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
    a_tr->trti = NULL;
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
tr_p_lowest_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring,
		    uint32_t *r_ntaxa, cw_tr_node_t a_root)
{
    cw_tr_node_t retval, node, root, troot;
    cw_tr_ring_t ring;
    cw_trn_t *trn;

    node = tr_p_ring_node_get(a_tr, a_ring);
    trn = &a_tr->trns[node];

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	/* Leaf node. */
	(*r_ntaxa)++;
    }

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE
	&& (a_root == CW_TR_NODE_NONE
	    || trn->taxon_num < a_tr->trns[a_root].taxon_num))
    {
	retval = node;
	root = node;
    }
    else
    {
	retval = CW_TR_NODE_NONE;
	root = a_root;
    }

    /* Iterate over neighbors. */
    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	troot = tr_p_lowest_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				    r_ntaxa, root);
	if (troot != CW_TR_NODE_NONE)
	{
	    retval = troot;
	    root = troot;
	}
    }

    return retval;
}

/* Recursively traverse the tree, count the number of taxa, and find the lowest
 * numbered taxon. */
static cw_tr_node_t
tr_p_lowest(cw_tr_t *a_tr, cw_tr_node_t a_node, uint32_t *r_ntaxa,
	    cw_tr_node_t a_root)
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
    qli_foreach(ring, &trn->rings, a_tr->trrs, link)
    {
	troot = tr_p_lowest_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				    r_ntaxa, root);
	if (troot != CW_TR_NODE_NONE)
	{
	    retval = troot;
	    root = troot;
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
	tr_p_lowest(a_tr, a_tr->base, &ntaxa, CW_TR_NODE_NONE);
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

CW_P_INLINE void
tr_p_trti_nodes_get(cw_tr_t *a_tr, uint32_t a_trti, cw_tr_node_t *r_node_a,
		    cw_tr_node_t *r_node_b)
{
    *r_node_a = tr_edge_node_get(a_tr, a_tr->trti[a_trti], 0);
    *r_node_b = tr_edge_node_get(a_tr, a_tr->trti[a_trti], 1);
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
	tr_p_lowest(a_tr, a_tr->base, &ntaxa, CW_TR_NODE_NONE);
    }

    // XXX This doesn't work for multifurcating trees.  Actually count the
    // edges in tr_p_lowest() recursion.
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
tr_p_bisection_edge_list_gen_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring,
				     uint32_t *ar_edges, uint32_t *ar_nedges)
{
    cw_tr_ring_t ring;

    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	/* Add edge to list. */
	ar_edges[*ar_nedges] = tr_p_ring_edge_get(a_tr, ring);
	(*ar_nedges)++;

	/* Recurse into neighbor subtree. */
	tr_p_bisection_edge_list_gen_recurse(a_tr,
					     tr_p_ring_other_get(a_tr, ring),
					     ar_edges, ar_nedges);
    }
}

/* Pretend that the tree is bisected at the edge that contains a_ring.
 * Construct a list of edges that are in the subtree that contains a_ring.
 *
 * The first element in the list is always the edge that is adjacent to the
 * bisection.  This facilitates recognition of reconnections that would reverse
 * bisection. */
CW_P_INLINE void
tr_p_bisection_edge_list_gen(cw_tr_t *a_tr, cw_tr_ring_t a_ring,
			     cw_tr_edge_t *ar_edges, uint32_t *ar_nedges)
{
    uint32_t degree;
    cw_tr_ring_t ring;

    /* Initialize the length of the list before recursing. */
    *ar_nedges = 0;

    /* Get the degree of the node adjacent to the bisection edge. */
    degree = 1;
    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	degree++;
    }

    switch (degree)
    {
	case 1:
	{
	    /* A subtree that is composed of a single node has no edges.  Add a
	     * single entry to the list. */
	    ar_edges[0] = CW_TR_EDGE_NONE;
	    (*ar_nedges)++;
	    break;
	}
	case 2:
	{
	    /* A tree should never have nodes of degree 2. */
	    cw_not_reached();
	}
	case 3:
	{
	    /* Take care to add only one of the edges that is connected to the
	     * node, since from the perspective of TBR, the node does not exist.
	     * (A node of degree 2 is a superfluous node.) */

	    /* First edge. */
	    ring = qri_next(a_tr->trrs, a_ring, link);
	    ar_edges[0] = tr_p_ring_edge_get(a_tr, ring);
	    (*ar_nedges)++;

	    /* First subtree. */
	    tr_p_bisection_edge_list_gen_recurse(a_tr,
						 tr_p_ring_other_get(a_tr,
								     ring),
						 ar_edges, ar_nedges);

	    /* Second subtree. */
	    ring = qri_next(a_tr->trrs, ring, link);
	    tr_p_bisection_edge_list_gen_recurse(a_tr,
						 tr_p_ring_other_get(a_tr,
								     ring),
						 ar_edges, ar_nedges);
	    break;
	}
	default:
	{
	    /* Add all edges in the subtree.  Removing the bisection edge still
	     * leaves enough edges attached to the node for the node to have
	     * relevance. */
	    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
	    {
		tr_p_bisection_edge_list_gen_recurse(a_tr,
						     tr_p_ring_other_get(a_tr,
									 ring),
						     ar_edges, ar_nedges);
	    }
	    break;
	}
    }
}

/* Generate lists of edges in each half of a logical bisection at edge
 * a_bisect. */
CW_P_INLINE void
tr_p_bedges_gen(cw_tr_t *a_tr, cw_tr_edge_t a_bisect)
{
    cw_dassert(tr_p_edge_validate(a_tr, a_bisect));

    tr_p_bisection_edge_list_gen(a_tr, tr_p_edge_ring_get(a_tr, a_bisect, 0),
				 a_tr->bedges, &a_tr->nbedges_a);
    tr_p_bisection_edge_list_gen(a_tr, tr_p_edge_ring_get(a_tr, a_bisect, 1),
				 &a_tr->bedges[a_tr->nbedges_a],
				 &a_tr->nbedges_b);
}

static void
tr_p_trti_update_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring,
			 uint32_t *ar_edge_count)
{
    cw_tr_ring_t ring;

    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	/* Record edge. */
	a_tr->trti[*ar_edge_count] = tr_p_ring_edge_get(a_tr, a_ring);

	/* Recurse into neighbor subtree. */
	tr_p_trti_update_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				 ar_edge_count);
    }
}

CW_P_INLINE void
tr_p_trti_update(cw_tr_t *a_tr, uint32_t a_nedges_prev)
{
    uint32_t edge_count;

    cw_assert(a_tr->modified == false);

    /* Make sure that the trti array is the right size. */
    if (a_tr->nedges > a_nedges_prev)
    {
	// XXX Why calloc/memset?
	if (a_tr->trti == NULL)
	{
	    cw_assert(a_nedges_prev == 0);

	    a_tr->trti = (cw_tr_edge_t *)
		cw_opaque_calloc(mema_calloc_get(a_tr->mema),
				 mema_arg_get(a_tr->mema),
				 a_tr->nedges,
				 sizeof(cw_tr_edge_t));
	}
	else
	{
	    a_tr->trti = (cw_tr_edge_t *)
		cw_opaque_realloc(mema_realloc_get(a_tr->mema),
				  mema_arg_get(a_tr->mema),
				  a_tr->trti,
				  sizeof(cw_tr_edge_t) * a_tr->nedges,
				  sizeof(cw_tr_edge_t) * a_nedges_prev);
	    memset(&a_tr->trti[a_nedges_prev], 0,
		   (a_tr->nedges - a_nedges_prev) * sizeof(cw_tr_edge_t));
	}
    }
    else if (a_tr->nedges < a_nedges_prev)
    {
	/* Shrink the array. */
	if (a_tr->nedges > 0)
	{
	    a_tr->trti = (cw_tr_edge_t *)
		cw_opaque_realloc(mema_realloc_get(a_tr->mema),
				  mema_arg_get(a_tr->mema),
				  a_tr->trti,
				  sizeof(cw_tr_edge_t) * a_tr->nedges,
				  sizeof(cw_tr_edge_t) * a_nedges_prev);
	}
	else
	{
	    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			      mema_arg_get(a_tr->mema),
			      a_tr->trti,
			      sizeof(cw_tr_edge_t) * a_nedges_prev);
	    a_tr->trti = NULL;
	}
    }

    /* Recursively traverse the tree, and initialize trti along the way. */
    if (a_tr->nedges > 0)
    {
	edge_count = 0;
	tr_p_trti_update_recurse(a_tr,
				 tr_p_ring_other_get(a_tr,
						     qli_first(&a_tr->trns
							       [a_tr->base]
							       .rings)),
				 &edge_count);
	cw_assert(edge_count == a_tr->nedges);
    }
}

static void
tr_p_trt_update(cw_tr_t *a_tr, uint32_t a_nedges_prev)
{
    uint32_t i, j, n, offset;

    cw_assert(a_tr->modified == false);

    /* Update trti. */
    tr_p_trti_update(a_tr, a_nedges_prev);

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

	/* Update bedges and trt. */
	tr_p_bedges_update(a_tr, nedges_prev);
	tr_p_trt_update(a_tr, nedges_prev);

	/* Clear held trees. */
	a_tr->nheld = 0;
    }
}

struct cw_tr_canonize_s
{
    cw_tr_ring_t ring;
    uint32_t min_taxon;
};

static int
tr_p_canonize_compar(const void *a_a, const void *a_b)
{
    const struct cw_tr_canonize_s *a = (const struct cw_tr_canonize_s *) a_a;
    const struct cw_tr_canonize_s *b = (const struct cw_tr_canonize_s *) a_b;

    if (a->min_taxon < b->min_taxon)
    {
	return -1;
    }
    else
    {
	cw_assert(a->min_taxon > b->min_taxon);
	return 1;
    }
}

/* Convert a tree to canonical form by re-ordering the ring such that subtrees
 * are in increasing order of minimum taxon number contained. */
static uint32_t
tr_p_canonize(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    uint32_t retval, degree, i, min_taxon;
    cw_tr_ring_t ring;
    struct cw_tr_canonize_s *canonize;

    /* Get taxon number (an internal node has CW_TR_NODE_TAXON_NONE). */
    cw_dassert(tr_p_node_validate(a_tr, tr_p_ring_node_get(a_tr, a_ring)));
    retval = tr_node_taxon_num_get(a_tr, tr_p_ring_node_get(a_tr, a_ring));

    /* Get the degree of the node that this ring is a part of. */
    degree = 1;
    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	degree++;
    }

    /* Allocate space for a temporary array that can be used to sort the
     * ring. */
    canonize = (struct cw_tr_canonize_s *)
	cw_opaque_alloc(mema_alloc_get(a_tr->mema),
			mema_arg_get(a_tr->mema),
			sizeof(struct cw_tr_canonize_s) * (degree - 1));

    /* Iteratively canonize subtrees, keeping track of the minimum taxon number
     * seen overall, as well as for each subtree. */
    for (i = 0, ring = a_ring; i < (degree - 1); i++)
    {
	ring = qri_next(a_tr->trrs, ring, link);

	min_taxon = tr_p_canonize(a_tr, tr_p_ring_other_get(a_tr, ring));
	if (min_taxon < retval)
	{
	    retval = min_taxon;
	}

	canonize[i].ring = ring;
	canonize[i].min_taxon = min_taxon;
    }

    /* Sort the subtrees. */
    qsort(canonize, degree - 1, sizeof(struct cw_tr_canonize_s),
	  tr_p_canonize_compar);

    /* Re-arrange the ring.  The first element can be skipped, since the
     * removal/re-insertion of all other elements eventually leaves the first
     * element in the proper location. */
    for (i = 1; i < (degree - 1); i++)
    {
	qri_remove(a_tr->trrs, canonize[i].ring, link);
	qri_before_insert(a_tr->trrs, a_ring, canonize[i].ring, link);
    }

    /* Clean up. */
    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
		      mema_arg_get(a_tr->mema),
		      canonize, sizeof(struct cw_tr_canonize_s) * (degree - 1));

    return retval;
}

CW_P_INLINE cw_tr_node_t
tr_p_tbr_node_extract(cw_tr_t *a_tr, cw_tr_node_t a_node,
		      cw_tr_edge_t a_reconnect_a, cw_tr_edge_t a_reconnect_b,
		      cw_tr_edge_t *ar_tedges, uint32_t *ar_ntedges,
		      cw_tr_node_t *ar_tnodes, uint32_t *ar_ntnodes)
{
    cw_tr_node_t retval;
    cw_tr_ring_t ring;
    cw_trn_t *trn;
    uint32_t degree;

    trn = &a_tr->trns[a_node];
    degree = 0;
    qli_foreach(ring, &trn->rings, a_tr->trrs, link)
    {
	degree++;
    }
    switch (degree)
    {
	case 0:
	{
	    /* This node is the only node remaining in the subtree.  It must be
	     * directly reconnected to, so return it. */
	    retval = a_node;
	    break;
	}
	case 1:
	{
	    cw_not_reached();
	}
	case 2:
	{
	    cw_tr_ring_t ring;
	    cw_tr_edge_t edge_lose, edge_keep;
	    cw_tr_node_t tnode_a, tnode_b;

	    ring = qli_first(&trn->rings);
	    edge_lose = tr_p_ring_edge_get(a_tr, ring);
	    tnode_a = tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring));

	    ring = qri_next(a_tr->trrs, ring, link);
	    edge_keep = tr_p_ring_edge_get(a_tr, ring);
	    tnode_b = tr_p_ring_node_get(a_tr, tr_p_ring_other_get(a_tr, ring));

	    if (tr_p_ring_edge_get(a_tr, ring) == a_reconnect_a
		|| tr_p_ring_edge_get(a_tr, ring) == a_reconnect_b)
	    {
		cw_tr_edge_t tedge;

		/* This edge is a reconnection edge; lose the other edge. */
		tedge = edge_keep;
		edge_keep = edge_lose;
		edge_lose = tedge;
	    }

	    /* Detach. */
	    tr_edge_detach(a_tr, edge_keep);
	    tr_edge_detach(a_tr, edge_lose);

	    /* Reattach. */
	    tr_edge_attach(a_tr, edge_keep, tnode_a, tnode_b);

	    /* Store spares. */
	    ar_tedges[*ar_ntedges] = edge_lose;
	    (*ar_ntedges)++;
	    ar_tnodes[*ar_ntnodes] = a_node;
	    (*ar_ntnodes)++;

	    retval = CW_TR_NODE_NONE;
	    break;
	}
	default:
	{
	    /* Do nothing. */
	    retval = CW_TR_NODE_NONE;
	}
    }

    return retval;
}

CW_P_INLINE cw_tr_node_t
tr_p_tbr_node_splice(cw_tr_t *a_tr, cw_tr_edge_t a_edge, 
		     cw_tr_edge_t *ar_tedges, uint32_t *ar_ntedges,
		     cw_tr_node_t *ar_tnodes, uint32_t *ar_ntnodes,
		     cw_tr_edge_t (*a_edge_alloc_callback)(cw_tr_t *, void *),
		     cw_tr_node_t (*a_node_alloc_callback)(cw_tr_t *, void *),
		     void *a_arg)
{
    cw_tr_node_t retval, node_a, node_b;
    cw_tr_edge_t edge;

    node_a = tr_edge_node_get(a_tr, a_edge, 0);
    cw_assert(node_a != CW_TR_NODE_NONE);
    node_b = tr_edge_node_get(a_tr, a_edge, 1);
    cw_assert(node_b != CW_TR_NODE_NONE);

    /* Get an edge. */
    if (*ar_ntedges > 0)
    {
	edge = ar_tedges[*ar_ntedges];
	(*ar_ntedges)--;
    }
    else
    {
	edge = a_edge_alloc_callback(a_tr, a_arg);
    }

    /* Get a node. */
    if (*ar_ntnodes > 0)
    {
	retval = ar_tnodes[*ar_ntnodes];
	(*ar_ntnodes)--;
    }
    else
    {
	retval = a_node_alloc_callback(a_tr, a_arg);
    }

    /* Detach. */
    tr_edge_detach(a_tr, a_edge);

    /* Reattach. */
    tr_edge_attach(a_tr, a_edge, retval, node_a);
    tr_edge_attach(a_tr, edge, retval, node_b);
    
    return retval;
}

static void
tr_p_mp_trn_prepare(cw_tr_t *a_tr, cw_trn_t *a_trn, char *a_taxa[],
		    uint32_t a_ntaxa, uint32_t a_nchars)
{
    uint32_t i, taxon_num;

    if (a_trn->ps == NULL)
    {
	a_trn->ps = tr_p_ps_new(a_tr);
    }
    tr_p_ps_prepare(a_tr, a_trn->ps, a_nchars);

    /* If this is a leaf node, initialize the character state sets and
     * scores. */
    if (a_trn->taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	char *chars;

	cw_assert(taxon_num < a_ntaxa);

	a_trn->ps->subtrees_score = 0;
	a_trn->ps->node_score = 0;

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
		    a_trn->ps->chars[i] = 0xf;
		    break;
		}
		case 'V':
		case 'v':
		{
		    a_trn->ps->chars[i] = 0xe;
		    break;
		}
		case 'H':
		case 'h':
		{
		    a_trn->ps->chars[i] = 0xd;
		    break;
		}
		case 'M':
		case 'm':
		{
		    a_trn->ps->chars[i] = 0xc;
		    break;
		}
		case 'D':
		case 'd':
		{
		    a_trn->ps->chars[i] = 0xb;
		    break;
		}
		case 'R':
		case 'r':
		{
		    a_trn->ps->chars[i] = 0xa;
		    break;
		}
		case 'W':
		case 'w':
		{
		    a_trn->ps->chars[i] = 0x9;
		    break;
		}
		case 'A':
		case 'a':
		{
		    a_trn->ps->chars[i] = 0x8;
		    break;
		}
		case 'B':
		case 'b':
		{
		    a_trn->ps->chars[i] = 0x7;
		    break;
		}
		case 'S':
		case 's':
		{
		    a_trn->ps->chars[i] = 0x6;
		    break;
		}
		case 'Y':
		case 'y':
		{
		    a_trn->ps->chars[i] = 0x5;
		    break;
		}
		case 'C':
		case 'c':
		{
		    a_trn->ps->chars[i] = 0x4;
		    break;
		}
		case 'K':
		case 'k':
		{
		    a_trn->ps->chars[i] = 0x3;
		    break;
		}
		case 'G':
		case 'g':
		{
		    a_trn->ps->chars[i] = 0x2;
		    break;
		}
		case 'T':
		case 't':
		{
		    a_trn->ps->chars[i] = 0x1;
		    break;
		}
		case '-':
		{
		    /* Treat gaps as uncertainty.  This isn't the only way to
		     * do things, and may need to be made configurable. */
		    a_trn->ps->chars[i] = 0xf;
		    break;
		}
		default:
		{
		    cw_not_reached();
		}
	    }
	}
    }
}

static void
tr_p_mp_prepare_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring,
			char *a_taxa[], uint32_t a_ntaxa, uint32_t a_nchars)
{
    cw_trn_t *trn;
    cw_tr_ring_t ring;
    cw_tre_t *tre;

    trn = &a_tr->trns[tr_p_ring_node_get(a_tr, a_ring)];

    tr_p_mp_trn_prepare(a_tr, trn, a_taxa, a_ntaxa, a_nchars);

    /* Recurse into subtrees. */
    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	/* Prepare edge before recursing. */
	tre = &a_tr->tres[tr_p_ring_edge_get(a_tr, ring)];
	if (tre->ps == NULL)
	{
	    tre->ps = tr_p_ps_new(a_tr);
	}
	tr_p_ps_prepare(a_tr, tre->ps, a_nchars);

	/* Recurse. */
	tr_p_mp_prepare_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				a_taxa, a_ntaxa, a_nchars);
    }
}

static void
tr_p_mp_trn_finish(cw_tr_t *a_tr, cw_trn_t *a_trn)
{
    if (a_trn->ps != NULL)
    {
	tr_p_ps_delete(a_tr, a_trn->ps);
	a_trn->ps = NULL;
    }
}

static void
tr_p_mp_finish_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    cw_trn_t *trn;
    cw_tr_ring_t ring;
    cw_tre_t *tre;

    trn = &a_tr->trns[tr_p_ring_node_get(a_tr, a_ring)];

    tr_p_mp_trn_finish(a_tr, trn);

    /* Recurse into subtrees. */
    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	/* Clean up edge before recursing. */
	tre = &a_tr->tres[tr_p_ring_edge_get(a_tr, ring)];
	if (tre->ps != NULL)
	{
	    tr_p_ps_delete(a_tr, tre->ps);
	    tre->ps = NULL;
	}

	/* Recurse. */
	tr_p_mp_finish_recurse(a_tr, tr_p_ring_other_get(a_tr, ring));
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
tr_p_mp_score_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring, cw_tr_edge_t a_bisect)
{
    uint32_t degree;
    bool adjacent;
    cw_tr_ring_t ring;

    /* Get the degree of the node.  Don't count the bisection edge (only an
     * issue if this node is adjacent to the bisection). */
    degree = 1;
    adjacent = false;
    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	if (tr_p_ring_edge_get(a_tr, ring) != a_bisect)
	{
	    degree++;
	}
	else
	{
	    adjacent = true;
	}
    }

    switch (degree)
    {
	case 1:
	{
	    /* Leaf node.  Do nothing. */
	    break;
	}
	case 2:
	{
	    cw_tr_ps_t *ps, *ps_c;
	    cw_tr_ring_t cring;

	    /* This is a trifurcating node that is adjacent to the bisection. */

	    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
	    {
		if (tr_p_ring_edge_get(a_tr, ring) != a_bisect)
		{
		    cring = ring;
		    break;
		}
	    }

	    ps_c = a_tr->trns[tr_p_ring_node_get(a_tr, cring)].ps;

	    /* Pass down the cached parent value, if the last time this node was
	     * recursed through, the bisection edge was the same. */
	    // XXX Does TBR break this?
	    ps = a_tr->trns[tr_p_ring_node_get(a_tr, a_ring)].ps;
	    if (ps->bisect == a_bisect)
	    {
		/* Pass cached parent value to the child. */
		ps_c->parent = ps->parent;
	    }
	    else
	    {
		/* Clear cached bisection edge. */
		ps->bisect = CW_TR_EDGE_NONE;
	    }

	    /* Recurse. */
	    tr_p_mp_score_recurse(a_tr, tr_p_ring_other_get(a_tr, cring),
				  a_bisect);

	    /* Copy child's scores to this node, rather than calculating a
	     * partial score.  This node is merely a filler node, as far as
	     * scoring is concerned. */
	    tr_p_mp_passpscore(a_tr, ps, ps_c);
	    break;
	}
	case 3:
	{
	    if (adjacent == false)
	    {
		cw_tr_ps_t *ps, *ps_a, *ps_b;

		/* This is a normal trifurcating node.  This is the common case,
		 * and is handled separately from the code below for performance
		 * reasons. */

		/* Recursively calculate partial scores for the subtrees. */
		ring = qri_next(a_tr->trrs, a_ring, link);
		ps_a = a_tr->trns[tr_p_ring_node_get(a_tr, ring)].ps;
		tr_p_mp_score_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				      a_bisect);

		ring = qri_next(a_tr->trrs, ring, link);
		ps_b = a_tr->trns[tr_p_ring_node_get(a_tr, ring)].ps;
		tr_p_mp_score_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				      a_bisect);

		/* Clear cached bisection edge. */
		ps = a_tr->trns[tr_p_ring_node_get(a_tr, a_ring)].ps;
		ps->bisect = CW_TR_EDGE_NONE;

		/* Calculate the partial score for this node. */
		tr_p_mp_pscore(a_tr, ps, ps_a, ps_b);

		break;
	    }
	    /* Fall through if this node is adjacent to the bisection. */
	}
	default:
	{
	    /* This is a multifurcating node. */
	    cw_error("XXX Not implemented");
	}
    }
}

static void
tr_p_mp_score(cw_tr_t *a_tr, cw_tr_edge_t a_root, cw_tr_edge_t a_bisect)
{
    cw_tr_edge_t edge_a, edge_b;

    edge_a = tr_p_edge_ring_get(a_tr, a_root, 0);
    tr_p_mp_score_recurse(a_tr, edge_a, a_bisect);

    edge_b = tr_p_edge_ring_get(a_tr, a_root, 1);
    tr_p_mp_score_recurse(a_tr, edge_b, a_bisect);

    /* Calculate the final score. */
    tr_p_mp_pscore(a_tr, a_tr->tres[a_root].ps,
		   a_tr->tres[edge_a].ps, a_tr->tres[edge_b].ps);
}

CW_P_INLINE void
tr_p_bisection_edge_list_mp(cw_tr_t *a_tr, uint32_t *a_edges,
			    uint32_t a_nedges, cw_tr_node_t a_other)
{
    uint32_t i;
    cw_tre_t *tre;

    for (i = 0; i < a_nedges; i++)
    {
	if (a_edges[i] != CW_TR_EDGE_NONE)
	{
	    tre = &a_tr->tresXXX[a_edges[i]];

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
	tre = &a_tr->tresXXX[i];

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
	    if (edge_a != CW_TR_EDGE_NONE)
	    {
		tre_a = &a_tr->tresXXX[edge_a];
		ps_a = a_tr->tresXXX[edge_a].ps;
	    }
	    else
	    {
		tre_a = NULL;
		ps_a = a_tr->trns[tre->node_a].ps;
	    }

	    for (k = 0; k < a_tr->nbedges_b; k++)
	    {
		edge_b = a_tr->bedges[a_tr->nbedges_a + k];
		if (edge_b != CW_TR_EDGE_NONE)
		{
		    tre_b = &a_tr->tresXXX[edge_b];
		    ps_b = a_tr->tresXXX[edge_b].ps;
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
	cw_opaque_dealloc(dealloc, arg, a_tr->trti,
			  sizeof(uint32_t) * (a_tr->nedges + 1));
    }

    if (a_tr->bedges != NULL)
    {
	cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			  mema_arg_get(a_tr->mema),
			  a_tr->bedges, sizeof(uint32_t) * a_tr->nedges);
    }

    if (a_tr->tresXXX != NULL)
    {
	uint32_t i;

	for (i = 0; i < a_tr->nedges; i++)
	{
	    if (a_tr->tresXXX[i].ps != NULL)
	    {
		tr_p_ps_delete(a_tr, a_tr->tresXXX[i].ps);
	    }
	}

	cw_opaque_dealloc(dealloc, arg, a_tr->tresXXX,
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
     * yet, since we will invalidate it during canonization. */
    if (a_tr->modified)
    {
	tr_p_ntaxa_nedges_update(a_tr);
	a_tr->modified = false;
    }

    if (a_tr->base != CW_TR_NODE_NONE)
    {
	uint32_t ntaxa;
	cw_tr_ring_t ring;

	/* Set base to be the lowest-numbered taxon. */
	ntaxa = 0;
	a_tr->base = tr_p_lowest(a_tr, a_tr->base, &ntaxa, CW_TR_NODE_NONE);

	/* Get base's ring. */
	ring = qli_first(&a_tr->tns[a_tr->base].rings);
	if (ring != CW_TR_RING_NONE)
	{
	    /* Canonize the tree. */
	    tr_p_canonize(a_tr, tr_p_ring_other_get(a_tr, ring));
	}
    }

    /* Now update trt. */
    tr_p_trt_update(a_tr, a_tr->nedges);

    cw_dassert(tr_p_validate(a_tr));
}

// XXX TBR can change the number of nodes and edges in a multifurcating tree.
// This throws a monkey wrench in things, since nodes and edges need to be
// allocated in Onyx code.
//
// One partial solution would be to pass spares in.  However, what do we do with
// unused or discarded nodes/edges?  We can't deallocate them (that's the
// GC's job), but leaving them laying around makes this code unsuitable for use
// outside the context of GC.  (Well, there is never a case where discarding
// nodes/edges is necessary.  Additionally, it's possible to return whether
// spares were actually used.)
//
// Another solution is to provide a C-level TBR API that returns pointers to any
// newly allocated nodes/edges.  These can then be wrapped by the glue code
// after the fact.  The main disadvantage I see for this is that there are GC
// races; a GC must not occur while there are nodes/edges without corresponding
// Onyx objects, since reference iteration becomes *much* harder to get right.
//
// Is there a reasonable way to provide an allocation callback function?  Yes,
// this can be done.  This is the cleanest solution, from an Onyx-level
// perspective.  It also avoids extra work in the common case (trifurcating
// trees).  It's rather gross (and slow) though.  A C function will have to call
// (assuming tree is on ostack):
//
//   cw_onyx_code(a_thread, "dup treenode:new");
//
void
tr_tbr(cw_tr_t *a_tr, cw_tr_edge_t a_bisect, cw_tr_edge_t a_reconnect_a,
       cw_tr_edge_t a_reconnect_b,
       cw_tr_edge_t (*a_edge_alloc_callback)(cw_tr_t *, void *),
       cw_tr_node_t (*a_node_alloc_callback)(cw_tr_t *, void *),
       void *a_arg)
{
    cw_tr_node_t node_a, node_b;
    cw_tr_edge_t tedges[3];
    uint32_t ntedges = 0
    cw_tr_node_t tnodes[2];
    uint32_t ntnodes = 0;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    /* Get the nodes to either side of the edge where the bisection will be
     * done. */
    node_a = tr_edge_node_get(a_tr, a_edge, 0);
    node_b = tr_edge_node_get(a_tr, a_edge, 1);

    /* Bisect and save edge as a spare. */
    tr_edge_detach(a_tr, a_edge);
    tedges[ntedges] = a_edge;
    ntedges++;

    /* For node_[ab], extract the node if it has only two neighbors.  If one of
     * the adjacent edges happens to be a reconnection edge, preserve it and
     * discard the other edge, so that reconnection to the edge works.
     *
     * The return value of these calls is CW_TR_NODE_NONE, unless there is only
     * one node in the subtree, in which case that node is returned so that it
     * can be used directly during reconnection. */
    node_a = tr_p_tbr_node_extract(a_tr, node_a, a_reconnect_a, a_reconnect_b,
				   tedges, &ntedges, tnodes &ntnodes);
    node_b = tr_p_tbr_node_extract(a_tr, node_b, a_reconnect_a, a_reconnect_b,
				   tedges, &ntedges, tnodes &ntnodes);

    /* For each reconnection edge, splice a node into the edge (if the subtree
     * has more than one node). */
    if (node_a == CW_TR_NODE_NONE)
    {
	node_a = tr_p_tbr_node_splice(a_tr, a_reconnect_a,
				      tedges, &ntedges, tnodes &ntnodes,
				      a_edge_alloc_callback,
				      a_node_alloc_callback,
				      a_arg);
    }
    if (node_b == CW_TR_NODE_NONE)
    {
	node_b = tr_p_tbr_node_splice(a_tr, a_reconnect_b,
				      tedges, &ntedges, tnodes &ntnodes,
				      a_edge_alloc_callback,
				      a_node_alloc_callback,
				      a_arg);
    }

    /* Attach the two spliced-in nodes. */
    if (ntedges > 0)
    {
	edge = tedges[ntedges - 1];
	ntedges--;
    }
    else
    {
	edge = a_edge_alloc_callback(a_tr, a_arg);
    }
    tr_edge_attach(a_tr, edge, node_a, node_b);

    /* Update. */
    tr_p_update(a_tr);
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
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    trn = a_tr->trns[a_tr->base];

    /* Prepare the tree. */
    tr_p_mp_trn_prepare(a_tr, trn, a_taxa, a_ntaxa, a_nchars);
    if (qli_first(&trn->rings) != CW_TR_RING_NONE)
    {
	qli_foreach(ring, &trn->rings, a_tr->trrs, link)
	{
	    tr_p_mp_prepare_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				    a_taxa, a_ntaxa, a_nchars);
	}
    }

    /* Prepare the temporary node. */
    tr_p_mp_trn_prepare(a_tr, &a_tr->trns[0], a_taxa, a_ntaxa, a_nchars);
}

void
tr_mp_finish(cw_tr_t *a_tr)
{
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    trn = a_tr->trns[a_tr->base];

    /* Clean up the tree. */
    tr_p_mp_trn_finish(a_tr, trn);
    if (qli_first(&trn->rings) != CW_TR_RING_NONE)
    {
	qli_foreach(ring, &trn->rings, a_tr->trrs, link)
	{
	    tr_p_mp_finish_recurse(a_tr, tr_p_ring_other_get(a_tr, ring))
	}
    }

    /* Clean up the temporary node. */
    tr_p_mp_trn_finish(a_tr, trn);
}

uint32_t
tr_mp_score(cw_tr_t *a_tr)
{
    uint32_t retval;
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    cw_dassert(tr_p_validate(a_tr));

    trn = &a_tr->trns[a_tr->base];

    if ((ring = qli_first(&trn->rings)) != CW_TR_RING_NONE)
    {
	cw_tr_edge_t edge;
	cw_tr_ps_t *ps;

	edge = tr_p_ring_edge_get(a_tr, ring)
	tr_p_mp_score(a_tr, edge, CW_TR_RING_NONE);

	ps = &a_tr->tres[edge].ps;
	retval = ps->subtrees_score + ps->node_score;
    }
    else
    {
	retval = 0;
    }

    return retval;
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
