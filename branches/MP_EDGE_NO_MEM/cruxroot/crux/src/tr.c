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
 * tr implements multifurcating phylogenetic trees and various operations on
 * them.  Nodes can be manipulated via the tr_node_*() APIs, and edges can be
 * manipulated via the tr_edge_*() APIs.
 *
 * Nodes and edges are allocated from per-tr internal arrays.  This allows array
 * indices to be used for links between nodes and edges.  The primary benefit,
 * however, is that tr_dup() does not need to traverse the tree.
 *
 * Each edge has two ring objects associated with it, one for each end of the
 * edge.  These ring objects are also allocated from an internal array that has
 * a direct correspondence to the edge array.  This provides for constant-time
 * conversion between edges and rings.
 *
 * Each node optionally references a ring, which contains all the edges that
 * connect the node to other nodes.  Each ring object in the ring refers to the
 * node, which allows constant-time conversion from rings to nodes.
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
#define CW_TR_RING_NONE UINT_MAX

/* Character (in the systematics sense of the word). */
typedef char cw_trc_t;

/* Partial parsimony score information. */
struct cw_tr_ps_s
{
    /* Sum of the subtree scores and this node's score, given particular
     * children. */
    uint32_t score;

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
	/* Auxiliary opaque data pointer. */
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
	/* Auxiliary opaque data pointer. */
	void *aux;

	/* Spares linkage. */
	cw_tr_edge_t link;
    } u;

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
    cw_tr_edge_t bisect_edge;
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

    /* Auxiliary opaque data pointer. */
    void *aux;

    /* True if this tree has been modified since the internal state (ntaxa,
     * nedges, bedges, trt, trti) was updated, false otherwise. */
    bool modified;

    /* Base of the tree (may or may not be set). */
    cw_tr_node_t base;

    /* Number of taxa in tree. */
    uint32_t ntaxa;

    /* Number of edges in tree.  This has to be calculated separately from
     * ntaxa, since there is no simple formula for nedges in multifurcating
     * trees (unlike for strictly trifurcating trees). */
    uint32_t nedges;

    /* bedges is an array of edges that is used for enumerating the edges on
     * each side of a logical tree bisection (used by TBR/MP-related functions).
     * The first list starts at offset 0 and has nbedges_a elements.  The second
     * list starts at offset nbedges_a and has nbedges_b elements. */
    cw_tr_edge_t *bedges;
    uint32_t nbedges_a;
    uint32_t nbedges_b;

    /* trt is an array of elements that store per-edge information that is used
     * for TBR-related functions.  There is one more element in trt than there
     * are edges in the tree.  This is critical to the way binary searching on
     * the array is done, and it also makes it easy to get the total number of
     * TBR neighbors this tree has (trt[nedges].offset).
     *
     * Only the first trtused elements are valid, since not all bisection edges
     * necessarily result in neighbors.
     *
     * trti is an array of edges that correspond to the entries in trt. */
    cw_trt_t *trt;
    uint32_t trtused;
    cw_tr_edge_t *trti;

    /* Pointer to an array of trn's.  ntrns is the total number of trn's, not
     * all of which are necessarily in use.
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

    /* held is an array of held TBR neighbors.  The array is iteratively doubled
     * as necessary.  heldlen is the actual length of the array, and nheld is
     * the number of elements in use. */
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

    cw_assert((a_ring >> 1) < a_tr->ntres);

    trr = &a_tr->trrs[a_ring];

    qri_new(a_tr->trrs, a_ring, link);
    trr->node = CW_TR_NODE_NONE;
    trr->ps = NULL;
}

CW_P_INLINE cw_tr_edge_t
tr_p_ring_edge_get(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    cw_assert((a_ring >> 1) < a_tr->ntres);

    return (a_ring >> 1);
}

CW_P_INLINE cw_tr_node_t
tr_p_ring_node_get(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    cw_assert((a_ring >> 1) < a_tr->ntres);

    return a_tr->trrs[a_ring].node;
}

CW_P_INLINE cw_tr_ring_t
tr_p_ring_other_get(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    cw_assert((a_ring >> 1) < a_tr->ntres);

    return (a_ring ^ 1);
}

/******************************************************************************/

/* Validation functions. */

#ifdef CW_DBG
/* Validate a ring. */
static bool
tr_p_ring_validate(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    cw_trr_t *trr;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
    cw_assert((a_ring >> 1) < a_tr->ntres);

    trr = &a_tr->trrs[a_ring];

    cw_assert(trr->node < a_tr->ntrns || trr->node == CW_TR_NODE_NONE);

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

    retval->chars = NULL;

    return retval;
}

CW_P_INLINE void
tr_p_ps_delete(cw_tr_t *a_tr, cw_tr_ps_t *a_ps)
{
    if (a_ps->chars != NULL)
    {
	cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			  mema_arg_get(a_tr->mema),
			  a_ps->achars,
			  sizeof(cw_trc_t) * ((a_ps->nchars >> 1)) + 8);
    }

    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
		      mema_arg_get(a_tr->mema),
		      a_ps, sizeof(cw_tr_ps_t));
}

#if (0) /* Unused (so far). */
CW_P_INLINE cw_trc_t
tr_p_ps_char_get(cw_tr_t *a_tr, cw_tr_ps_t *a_ps, uint32_t a_offset)
{
    cw_trc_t retval;

    cw_check_ptr(a_ps->chars);
    cw_assert(a_offset < a_ps->nchars);

    retval = a_ps->chars[a_offset >> 1];
    retval >>= ((a_offset & 1) * 4);

    return retval;
}
#endif

CW_P_INLINE void
tr_p_ps_char_set(cw_tr_t *a_tr, cw_tr_ps_t *a_ps, cw_trc_t a_char,
		 uint32_t a_offset)
{
    cw_check_ptr(a_ps->chars);
    cw_assert((a_char & 0xfU) == a_char);
    cw_assert(a_offset < a_ps->nchars);

    if ((a_offset & 1) == 0)
    {
	a_ps->chars[a_offset >> 1]
	    = (a_char << 4)
	    | (a_ps->chars[a_offset >> 1] & 0xfU);
    }
    else
    {
	a_ps->chars[a_offset >> 1]
	    = (a_ps->chars[a_offset >> 1] & 0xf0U)
	    | a_char;
    }
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
			  sizeof(cw_trc_t) * ((a_ps->nchars >> 1) + 8));
	a_ps->chars = NULL;
    }

    /* Allocate character vector if necessary. */
    if (a_ps->chars == NULL)
    {
	uint32_t npad, i;

	/* Calculate the number of pad characters to append, such that the total
	 * number of characters is a multiple of 16. */
	npad = (32 - (a_nchars & 0x1f)) & 0xfU;
	cw_assert(((a_nchars + npad) & 0x1f) == 0);

	/* Tack on 8 bytes; all modern systems provide at least 8 byte
	 * alignment. */
	a_ps->achars
	    = (cw_trc_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					   mema_arg_get(a_tr->mema),
					   sizeof(cw_trc_t)
					   * (((a_nchars + npad) >> 1)) + 8);

	/* Make sure that chars is 16 byte-aligned.  Assume that chars is at
	 * least 8 byte-aligned. */
	a_ps->chars = &a_ps->achars[((unsigned) a_ps->achars) & 0xfU];

	/* Set all pad characters to {ACGT}.  This allows the pad characters to
	 * be calculated along with the actual characters, without affecting the
	 * score. */
	for (i = a_nchars; i < a_nchars + npad; i++)
	{
	    tr_p_ps_char_set(a_tr, a_ps, 0xfU, i);
	}

	a_ps->nchars = a_nchars + npad;
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
#ifdef CW_DBG
    cw_tre_t *tre;
#endif
    cw_trr_t *trr;

#ifdef CW_DBG
    tre = &a_tr->tres[a_edge];
    memset(tre, 0x5a, sizeof(cw_tre_t));
#endif

    trr = &a_tr->trrs[a_edge << 1];
    if (trr->ps != NULL)
    {
	tr_p_ps_delete(a_tr, trr->ps);
    }
#ifdef CW_DBG
    memset(trr, 0x5a, sizeof(cw_trr_t));
#endif

    trr = &a_tr->trrs[(a_edge << 1) + 1];
    if (trr->ps != NULL)
    {
	tr_p_ps_delete(a_tr, trr->ps);
    }
#ifdef CW_DBG
    memset(trr, 0x5a, sizeof(cw_trr_t));
#endif

    a_tr->tres[a_edge].u.link = a_tr->sparetres;
    a_tr->sparetres = a_edge;
}

CW_P_INLINE cw_tr_ring_t
tr_p_edge_ring_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_end)
{
    return ((a_edge << 1) + a_end);
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
    cw_assert(tr_p_ring_node_get(a_tr, tr_p_edge_ring_get(a_tr, a_edge, 0))
	      == CW_TR_NODE_NONE);
    cw_assert(tr_p_ring_node_get(a_tr, tr_p_edge_ring_get(a_tr, a_edge, 1))
	      == CW_TR_NODE_NONE);

    tr_p_edge_dealloc(a_tr, a_edge);
}

cw_tr_node_t
tr_edge_node_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_end)
{
    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_assert(a_end == 0 || a_end == 1);

    return tr_p_ring_node_get(a_tr, tr_p_edge_ring_get(a_tr, a_edge, a_end));
}

void
tr_edge_next_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_end,
		 cw_tr_edge_t *r_next, uint32_t *r_end)
{
    cw_tr_ring_t ringind;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_assert(a_end == 0 || a_end == 1);

    ringind = qri_next(a_tr->trrs, tr_p_edge_ring_get(a_tr, a_edge, a_end),
		       link);
    *r_next = tr_p_ring_edge_get(a_tr, ringind);
    *r_end = (ringind & 1);
}

void
tr_edge_prev_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_end,
		 cw_tr_edge_t *r_prev, uint32_t *r_end)
{
    cw_tr_ring_t ringind;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_assert(a_end == 0 || a_end == 1);

    ringind = qri_prev(a_tr->trrs, tr_p_edge_ring_get(a_tr, a_edge, a_end),
		       link);
    *r_prev = tr_p_ring_edge_get(a_tr, ringind);
    *r_end = (ringind & 1);
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
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 0) == CW_TR_NODE_NONE);
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 1) == CW_TR_NODE_NONE);
    cw_dassert(tr_p_node_validate(a_tr, a_node_a));
    cw_dassert(tr_p_node_validate(a_tr, a_node_b));
    cw_assert(a_node_a != a_node_b);
    cw_assert(tr_node_distance(a_tr, a_node_a, a_node_b) == 0);

    /* First end. */
    ring = tr_p_edge_ring_get(a_tr, a_edge, 0);
    trn = &a_tr->trns[a_node_a];
    qli_tail_insert(&trn->rings, a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = a_node_a;

    /* Second end. */
    ring = tr_p_edge_ring_get(a_tr, a_edge, 1);
    trn = &a_tr->trns[a_node_b];
    qli_tail_insert(&trn->rings, a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = a_node_b;

    /* Mark tree as modified. */
    a_tr->modified = true;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_dassert(tr_p_node_validate(a_tr, a_node_a));
    cw_dassert(tr_p_node_validate(a_tr, a_node_b));
    cw_assert(tr_node_distance(a_tr, a_node_a, a_node_b) == 1);
}

void
tr_edge_detach(cw_tr_t *a_tr, cw_tr_edge_t a_edge)
{
    cw_tr_ring_t ring;
    cw_trn_t *trn;

    cw_dassert(tr_p_edge_validate(a_tr, a_edge));
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 0) != CW_TR_NODE_NONE);
    cw_dassert(tr_edge_node_get(a_tr, a_edge, 1) != CW_TR_NODE_NONE);
    cw_assert(tr_node_distance(a_tr, tr_edge_node_get(a_tr, a_edge, 0),
			       tr_edge_node_get(a_tr, a_edge, 1)) == 1);

    /* Detach from neighboring nodes.  Use qli_remove() to make sure that the
     * nodes still point to their rings. */
    ring = tr_p_edge_ring_get(a_tr, a_edge, 0);
    trn = &a_tr->trns[tr_p_ring_node_get(a_tr, ring)];
    qli_remove(&trn->rings, a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = CW_TR_NODE_NONE;

    ring = tr_p_edge_ring_get(a_tr, a_edge, 1);
    trn = &a_tr->trns[tr_p_ring_node_get(a_tr, ring)];
    qli_remove(&trn->rings, a_tr->trrs, ring, link);
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
	if (a_tr->trns == NULL)
	{
	    a_tr->trns
		= (cw_trn_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					       mema_arg_get(a_tr->mema),
					       sizeof(cw_trn_t));
	    nspares = 1;
	    a_tr->ntrns = 1;
	}
	else
	{
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
	}

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
    memset(trn, 0x5a, sizeof(cw_trn_t));
#endif

    a_tr->trns[a_node].u.link = a_tr->sparetrns;
    a_tr->sparetrns = a_node;
}

/* Calculate the number of edges connected to the node that a_ring is connected
 * to. */
CW_P_INLINE uint32_t
tr_p_node_degree(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    uint32_t retval;
    cw_tr_ring_t ring;

    retval = 1;
    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	retval++;
    }

    return retval;
}

/* Calculate the number of edges between two nodes.  A distance of 0 means that
 * there is no path between the two nodes. */
static uint32_t
tr_p_node_distance(cw_tr_t *a_tr, cw_tr_ring_t a_ring, cw_tr_node_t a_other,
		   uint32_t a_distance)
{
    uint32_t retval;
    cw_tr_ring_t ring;

    if (tr_p_ring_node_get(a_tr, a_ring) == a_other)
    {
	retval = a_distance;
	goto RETURN;
    }

    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	if ((retval = tr_p_node_distance(a_tr, tr_p_ring_other_get(a_tr, ring),
					 a_other, a_distance + 1)) != 0)
	{
	    goto RETURN;
	}
    }

    retval = 0;
    RETURN:
    return retval;
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
		 cw_tr_edge_t *r_edge, uint32_t *r_end)
{
    cw_trn_t *trn;
    cw_tr_ring_t ringind;

    cw_dassert(tr_p_node_validate(a_tr, a_node));

    trn = &a_tr->trns[a_node];

    ringind = qli_first(&trn->rings);
    if (ringind != CW_TR_EDGE_NONE)
    {
	*r_edge = tr_p_ring_edge_get(a_tr, ringind);
	if (r_end != NULL)
	{
	    *r_end = (ringind & 1);
	}
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

uint32_t
tr_node_degree(cw_tr_t *a_tr, cw_tr_node_t a_node)
{
    uint32_t retval;
    cw_tr_ring_t ring;

    cw_dassert(tr_p_node_validate(a_tr, a_node));

    ring = qli_first(&a_tr->trns[a_node].rings);
    if (ring != CW_TR_RING_NONE)
    {
	retval = tr_p_node_degree(a_tr, ring);
    }
    else
    {
	retval = 0;
    }

    return retval;
}

uint32_t
tr_node_distance(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_other)
{
    uint32_t retval;
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    cw_dassert(tr_p_node_validate(a_tr, a_node));
    cw_dassert(tr_p_node_validate(a_tr, a_other));
    cw_assert(a_node != a_other);

    trn = &a_tr->trns[a_node];

    ring = qli_first(&trn->rings);
    if (ring != CW_TR_RING_NONE)
    {
	qli_foreach(ring, &trn->rings, a_tr->trrs, link)
	{
	    if ((retval = tr_p_node_distance(a_tr,
					     tr_p_ring_other_get(a_tr, ring),
					     a_other, 1)) != 0)
	    {
		break;
	    }
	}
    }
    else
    {
	retval = 0;
    }

    return retval;
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
    a_tr->bedges = NULL;
    a_tr->nbedges_a = 0;
    a_tr->nbedges_b = 0;
    a_tr->trt = NULL;
    a_tr->trtused = 0;
    a_tr->trti = NULL;
    a_tr->trns = NULL;
    a_tr->ntrns = 0;
    a_tr->sparetrns = CW_TR_NODE_NONE;
    a_tr->tres = NULL;
    a_tr->ntres = 0;
    a_tr->sparetres = CW_TR_EDGE_NONE;
    a_tr->trrs = NULL;
    a_tr->held = NULL;
    a_tr->heldlen = 0;
    a_tr->nheld = 0;

#ifdef CW_DBG
    a_tr->magic = CW_TR_MAGIC;
#endif
}

/* Recursively traverse the tree, count the number of taxa, and find the lowest
 * numbered taxon. */
static cw_tr_node_t
tr_p_lowest_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring, uint32_t *r_ntaxa,
		    uint32_t *r_nedges, cw_tr_node_t a_root)
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
	/* Count edge. */
	(*r_nedges)++;

	troot = tr_p_lowest_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				    r_ntaxa, r_nedges, root);
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
	    uint32_t *r_nedges)
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

    if (trn->taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	retval = a_node;
	root = a_node;
    }
    else
    {
	retval = CW_TR_NODE_NONE;
	root = CW_TR_NODE_NONE;
    }

    /* Iterate over neighbors. */
    qli_foreach(ring, &trn->rings, a_tr->trrs, link)
    {
	/* Count edge. */
	(*r_nedges)++;
	
	troot = tr_p_lowest_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				    r_ntaxa, r_nedges, root);
	if (troot != CW_TR_NODE_NONE)
	{
	    retval = troot;
	    root = troot;
	}
    }

    return retval;
}

#ifdef CW_DBG
/* Validate a tree. */
static bool
tr_p_validate(cw_tr_t *a_tr)
{
    uint32_t i, ntaxa, nedges;

    cw_check_ptr(a_tr);
    cw_assert(a_tr->magic == CW_TR_MAGIC);
    cw_assert(a_tr->modified == false);

    ntaxa = 0;
    nedges = 0;
    if (a_tr->base != CW_TR_NODE_NONE)
    {
	tr_p_lowest(a_tr, a_tr->base, &ntaxa, &nedges);
    }
    cw_assert(a_tr->ntaxa == ntaxa);
    cw_assert(a_tr->nedges == nedges);

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

static void
tr_p_ntaxa_nedges_update(cw_tr_t *a_tr)
{
    uint32_t ntaxa, nedges;

    /* Update ntaxa and nedges. */
    ntaxa = 0;
    nedges = 0;
    if (a_tr->base != CW_TR_NODE_NONE)
    {
	tr_p_lowest(a_tr, a_tr->base, &ntaxa, &nedges);
    }

    a_tr->ntaxa = ntaxa;
    a_tr->nedges = nedges;
}

static void
tr_p_bisection_edge_list_gen_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring,
				     cw_tr_edge_t *ar_edges,
				     uint32_t *ar_nedges)
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
 * bisection.
 *
 * If the list is empty (bisection adjacent to a leaf node), return the single
 * node, so that it can be accessed directly (there's no edge logically attached
 * to it). */
CW_P_INLINE cw_tr_node_t
tr_p_bisection_edge_list_gen(cw_tr_t *a_tr, cw_tr_ring_t a_ring,
			     cw_tr_edge_t *ar_edges, uint32_t *ar_nedges)
{
    cw_tr_node_t retval;
    cw_tr_ring_t ring;

    /* Initialize the length of the list before recursing. */
    *ar_nedges = 0;

    switch (tr_p_node_degree(a_tr, a_ring))
    {
	case 1:
	{
	    /* A subtree that is composed of a single node has no edges.  Add a
	     * single entry to the list, and return the node. */
	    ar_edges[0] = CW_TR_EDGE_NONE;
	    (*ar_nedges)++;
	    retval = tr_p_ring_node_get(a_tr, a_ring);
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

	    retval = CW_TR_NODE_NONE;
	    break;
	}
	default:
	{
	    /* Add all edges in the subtree.  Removing the bisection edge still
	     * leaves enough edges attached to the node for the node to have
	     * relevance. */
	    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
	    {
		/* Add edge to list. */
		ar_edges[*ar_nedges] = tr_p_ring_edge_get(a_tr, ring);
		(*ar_nedges)++;

		tr_p_bisection_edge_list_gen_recurse(a_tr,
						     tr_p_ring_other_get(a_tr,
									 ring),
						     ar_edges, ar_nedges);
	    }

	    retval = CW_TR_NODE_NONE;
	    break;
	}
    }

    return retval;
}

/* Generate lists of edges in each half of a logical bisection at edge
 * a_bisect. */
CW_P_INLINE void
tr_p_bedges_gen(cw_tr_t *a_tr, cw_tr_edge_t a_bisect, cw_tr_node_t *r_node_a,
		cw_tr_node_t *r_node_b)
{
    cw_tr_node_t node_a, node_b;

    cw_dassert(tr_p_edge_validate(a_tr, a_bisect));

    node_a = tr_p_bisection_edge_list_gen(a_tr,
					  tr_p_edge_ring_get(a_tr, a_bisect, 0),
					  a_tr->bedges, &a_tr->nbedges_a);
    node_b = tr_p_bisection_edge_list_gen(a_tr,
					  tr_p_edge_ring_get(a_tr, a_bisect, 1),
					  &a_tr->bedges[a_tr->nbedges_a],
					  &a_tr->nbedges_b);

    if (r_node_a != NULL)
    {
	*r_node_a = node_a;
    }
    if (r_node_b != NULL)
    {
	*r_node_b = node_b;
    }
}

static void
tr_p_trti_update_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring,
			 uint32_t *ar_edge_count)
{
    cw_tr_ring_t ring;

    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	/* Record edge. */
	a_tr->trti[*ar_edge_count] = tr_p_ring_edge_get(a_tr, ring);
	(*ar_edge_count)++;

	/* Recurse into neighbor subtree. */
	tr_p_trti_update_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				 ar_edge_count);
    }
}

CW_P_INLINE void
tr_p_trti_update(cw_tr_t *a_tr, uint32_t a_nedges_prev)
{
    cw_assert(a_tr->modified == false);

    /* Make sure that the trti array is the right size. */
    if (a_tr->nedges > a_nedges_prev)
    {
	/* Grow the array. */
	if (a_tr->trti == NULL)
	{
	    cw_assert(a_nedges_prev == 0);

	    a_tr->trti = (cw_tr_edge_t *)
		cw_opaque_alloc(mema_alloc_get(a_tr->mema),
				mema_arg_get(a_tr->mema),
				sizeof(cw_tr_edge_t) * a_tr->nedges);
	}
	else
	{
	    a_tr->trti = (cw_tr_edge_t *)
		cw_opaque_realloc(mema_realloc_get(a_tr->mema),
				  mema_arg_get(a_tr->mema),
				  a_tr->trti,
				  sizeof(cw_tr_edge_t) * a_tr->nedges,
				  sizeof(cw_tr_edge_t) * a_nedges_prev);
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
	uint32_t edge_count;
	cw_trn_t *trn;
	cw_tr_ring_t ring;

	cw_assert(a_tr->base != CW_TR_NODE_NONE);

	edge_count = 0;
	trn = &a_tr->trns[a_tr->base];
	qli_foreach(ring, &trn->rings, a_tr->trrs, link)
	{
	    /* Record edge. */
	    a_tr->trti[edge_count] = tr_p_ring_edge_get(a_tr, ring);
	    edge_count++;

	    tr_p_trti_update_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				     &edge_count);
	}
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
	/* Reallocate trt.  There is never a need to deallocate trt here,
	 * since trt contains one extra element. */
	a_tr->trt = (cw_trt_t *) cw_opaque_realloc(mema_realloc_get(a_tr->mema),
						   mema_arg_get(a_tr->mema),
						   a_tr->trt,
						   sizeof(cw_trt_t)
						   * (a_tr->nedges + 1),
						   sizeof(cw_trt_t)
						   * (a_nedges_prev + 1));
    }

    /* Iteratively fill in trt. */
    for (i = j = offset = 0; i < a_tr->nedges; i++)
    {
	/* Record offset. */
	a_tr->trt[j].offset = offset;

	/* Record bisection edge. */
	a_tr->trt[j].bisect_edge = a_tr->trti[i];

	/* Update offset. */
	tr_p_bedges_gen(a_tr, a_tr->trti[i], NULL, NULL);
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

/* trt comparison function passed to bsearch(). */
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
	    = (cw_tr_edge_t *) cw_opaque_alloc(mema_alloc_get(a_tr->mema),
					       mema_arg_get(a_tr->mema),
					       sizeof(cw_tr_edge_t)
					       * a_tr->nedges);
    }
    else if (a_tr->nedges != a_nedges_prev)
    {
	if (a_tr->nedges > 0)
	{
	    /* Reallocate bedges. */
	    a_tr->bedges = (cw_tr_edge_t *)
		cw_opaque_realloc(mema_realloc_get(a_tr->mema),
				  mema_arg_get(a_tr->mema),
				  a_tr->bedges,
				  sizeof(cw_tr_edge_t) * a_tr->nedges,
				  sizeof(cw_tr_edge_t) * a_nedges_prev);
	}
	else
	{
	    /* Deallocate bedges. */
	    cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			      mema_arg_get(a_tr->mema),
			      a_tr->bedges,
			      sizeof(cw_tr_edge_t) * a_nedges_prev);
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

/* Used for canonizing trees. */
struct cw_tr_canonize_s
{
    cw_tr_ring_t ring;
    uint32_t min_taxon;
};

/* Comparison function that is passed to qsort(). */
static int
tr_p_canonize_compare(const void *a_a, const void *a_b)
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

/* Convert a tree node to canonical form by re-ordering the ring such that
 * subtrees are in increasing order of minimum taxon number contained. */
static uint32_t
tr_p_canonize(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    uint32_t retval, degree;
    cw_tr_node_t node;

    /* Get taxon number (an internal node has CW_TR_NODE_TAXON_NONE). */
    cw_dassert(tr_p_node_validate(a_tr, tr_p_ring_node_get(a_tr, a_ring)));
    node = tr_p_ring_node_get(a_tr, a_ring);
    retval = tr_node_taxon_num_get(a_tr, node);

    /* Get the degree of the node that this ring is a part of. */
    degree = tr_p_node_degree(a_tr, a_ring);

    if (degree > 1)
    {
	uint32_t i, min_taxon;
	cw_tr_ring_t ring;
	struct cw_tr_canonize_s *canonize;

	/* Allocate space for a temporary array that can be used to sort the
	 * ring. */
	canonize = (struct cw_tr_canonize_s *)
	    cw_opaque_alloc(mema_alloc_get(a_tr->mema),
			    mema_arg_get(a_tr->mema),
			    sizeof(struct cw_tr_canonize_s) * (degree - 1));

	/* Iteratively canonize subtrees, keeping track of the minimum taxon
	 * number seen overall, as well as for each subtree. */
	i = 0;
	retval = CW_TR_NODE_TAXON_NONE;
	qri_others_foreach(ring, a_tr->trrs, a_ring, link)
	{
	    min_taxon = tr_p_canonize(a_tr, tr_p_ring_other_get(a_tr, ring));
	    if (min_taxon < retval)
	    {
		retval = min_taxon;
	    }

	    canonize[i].ring = ring;
	    canonize[i].min_taxon = min_taxon;

	    i++;
	}
	cw_assert(i == degree - 1);

	/* Sort the subtrees. */
	qsort(canonize, degree - 1, sizeof(struct cw_tr_canonize_s),
	      tr_p_canonize_compare);

	/* Set the beginning of the ring to a_ring.  This makes it easier for
	 * external code to traverse a tree in canonical order. */
	qli_first(&a_tr->trns[node].rings) = a_ring;

	/* Re-arrange the ring.  The first element can be skipped, since the
	 * removal/re-insertion of all other elements eventually leaves the
	 * first element in the proper location. */
	for (i = 1; i < (degree - 1); i++)
	{
	    qri_remove(a_tr->trrs, canonize[i].ring, link);
	    qri_before_insert(a_tr->trrs, a_ring, canonize[i].ring, link);
	}

	/* Clean up. */
	cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			  mema_arg_get(a_tr->mema),
			  canonize,
			  sizeof(struct cw_tr_canonize_s) * (degree - 1));
    }

    return retval;
}

/* As part of TBR, extract a node that has only two neighbors.  Take care to
 * leave reconnection edges in the tree.  Return CW_TR_NODE_NONE, unless there
 * is only one node in the subtree; in that case, return the node so that it can
 * be used directly during reconnection. */
CW_P_INLINE cw_tr_node_t
tr_p_tbr_node_extract(cw_tr_t *a_tr, cw_tr_node_t a_node,
		      cw_tr_edge_t a_reconnect_a, cw_tr_edge_t a_reconnect_b,
		      cw_tr_edge_t *ar_tedges, uint32_t *ar_ntedges,
		      cw_tr_node_t *ar_tnodes, uint32_t *ar_ntnodes)
{
    cw_tr_node_t retval;

    switch (tr_node_degree(a_tr, a_node))
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
	    cw_tr_ring_t ring_a, ring_b;
	    cw_tr_ring_t ring_a_other, ring_b_other;
	    cw_tr_edge_t edge_a, edge_b;
	    cw_tr_node_t node_a, node_b;
	    cw_tr_ps_t *tps;

	    /* Get all variables that are necessary for careful extraction of
	     * a_node, and proper rematching of rings with nodes.  The
	     * rematching is critical to the maintenance of the character state
	     * sets in leaf nodes (which node_[ab] may or may not be). */
	    ring_a = qli_first(&a_tr->trns[a_node].rings);
	    edge_a = tr_p_ring_edge_get(a_tr, ring_a);
	    ring_a_other = tr_p_ring_other_get(a_tr, ring_a);
	    node_a = tr_p_ring_node_get(a_tr, ring_a_other);

	    ring_b = qri_next(a_tr->trrs, ring_a, link);
	    edge_b = tr_p_ring_edge_get(a_tr, ring_b);
	    ring_b_other = tr_p_ring_other_get(a_tr, ring_b);
	    node_b = tr_p_ring_node_get(a_tr, ring_b_other);

	    /* Detach. */
	    tr_edge_detach(a_tr, edge_a);
	    tr_edge_detach(a_tr, edge_b);

	    /* Store a_node as a spare. */
	    ar_tnodes[*ar_ntnodes] = a_node;
	    (*ar_ntnodes)++;

	    /* Be careful to preserve reconnection edges, which either edge_a or
	     * edge_b may be. */
	    if (edge_b != a_reconnect_a && edge_b != a_reconnect_b)
	    {
		/* Use edge_a when splicing node_[ab] back together. */

		/* Swap data in ring_a and ring_b_other. */
		tps = a_tr->trrs[ring_a].ps;
		a_tr->trrs[ring_a].ps = a_tr->trrs[ring_b_other].ps;
		a_tr->trrs[ring_b_other].ps = tps;

		/* Attach node_[ab].  Take care to keep the proper ends of
		 * edge_a associated with node_[ab]. */
		if (ring_a_other < ring_a)
		{
		    tr_edge_attach(a_tr, edge_a, node_a, node_b);
		}
		else
		{
		    tr_edge_attach(a_tr, edge_a, node_b, node_a);
		}

		/* Store edge_b as a spare. */
		ar_tedges[*ar_ntedges] = edge_b;
		(*ar_ntedges)++;
	    }
	    else
	    {
		/* Use edge_b when splicing node_[ab] back together. */
		cw_assert(edge_a != a_reconnect_a && edge_a != a_reconnect_b);

		/* Swap data in ring_b and ring_a_other. */
		tps = a_tr->trrs[ring_b].ps;
		a_tr->trrs[ring_b].ps = a_tr->trrs[ring_a_other].ps;
		a_tr->trrs[ring_a_other].ps = tps;

		/* Attach node_[ab].  Take care to keep the proper ends of
		 * edge_b associated with node_[ab]. */
		if (ring_b < ring_b_other)
		{
		    tr_edge_attach(a_tr, edge_b, node_a, node_b);
		}
		else
		{
		    tr_edge_attach(a_tr, edge_b, node_b, node_a);
		}

		/* Store edge_a as a spare. */
		ar_tedges[*ar_ntedges] = edge_a;
		(*ar_ntedges)++;
	    }

	    retval = CW_TR_NODE_NONE;
	    break;
	}
	default:
	{
	    /* Do nothing, since this node has enough neighbors to remain
	     * relevant (3 or more). */
	    retval = CW_TR_NODE_NONE;
	}
    }

    return retval;
}

/* Splice a node into the middle of a_edge, and return the node. */
CW_P_INLINE cw_tr_node_t
tr_p_tbr_node_splice(cw_tr_t *a_tr, cw_tr_edge_t a_edge, 
		     cw_tr_edge_t *ar_tedges, uint32_t *ar_ntedges,
		     cw_tr_node_t *ar_tnodes, uint32_t *ar_ntnodes,
		     cw_tr_edge_t (*a_edge_alloc_callback)(cw_tr_t *, void *),
		     cw_tr_node_t (*a_node_alloc_callback)(cw_tr_t *, void *),
		     void *a_arg)
{
    cw_tr_node_t retval, node_a, node_b;
    cw_tr_ring_t ring_a, ring_b, ring;
    cw_tr_edge_t edge;
    cw_tr_ps_t *tps;

    /* Get all variables that are necessary for careful splicing of a node into
     * a_edge, and proper rematching of rings with nodes.  The rematching is
     * critical to the maintenance of the character state sets in leaf nodes
     * (which node_[ab] may or may not be). */
    ring_a = tr_p_edge_ring_get(a_tr, a_edge, 0);
    node_a = tr_p_ring_node_get(a_tr, ring_a);

    ring_b = tr_p_edge_ring_get(a_tr, a_edge, 1);
    node_b = tr_p_ring_node_get(a_tr, ring_b);

    /* Get an edge. */
    if (*ar_ntedges > 0)
    {
	(*ar_ntedges)--;
	edge = ar_tedges[*ar_ntedges];
    }
    else
    {
	edge = a_edge_alloc_callback(a_tr, a_arg);
    }
    ring = tr_p_edge_ring_get(a_tr, edge, 0);

    /* Get a node. */
    if (*ar_ntnodes > 0)
    {
	(*ar_ntnodes)--;
	retval = ar_tnodes[*ar_ntnodes];
    }
    else
    {
	retval = a_node_alloc_callback(a_tr, a_arg);
    }

    /* Detach. */
    tr_edge_detach(a_tr, a_edge);

    /* Swap data in ring_b and ring. */
    tps = a_tr->trrs[ring_b].ps;
    a_tr->trrs[ring_b].ps = a_tr->trrs[ring].ps;
    a_tr->trrs[ring].ps = tps;

    /* Reattach. */
    tr_edge_attach(a_tr, a_edge, node_a, retval);
    tr_edge_attach(a_tr, edge, node_b, retval);

    return retval;
}

static void
tr_p_mp_ring_prepare(cw_tr_t *a_tr, cw_tr_ring_t a_ring, char *a_taxa[],
		     uint32_t a_ntaxa, uint32_t a_nchars,
		     bool *a_chars_mask, uint32_t a_ninformative)
{
    cw_trr_t *trr;
    uint32_t taxon_num;

    trr = &a_tr->trrs[a_ring];

    if (trr->ps == NULL)
    {
	trr->ps = tr_p_ps_new(a_tr);
    }
    tr_p_ps_prepare(a_tr, trr->ps, a_ninformative);

    /* If this is a leaf node, initialize the character state sets and
     * scores. */
    taxon_num = a_tr->trns[trr->node].taxon_num;
    if (taxon_num != CW_TR_NODE_TAXON_NONE)
    {
	uint32_t i, j;
	char *chars;

	trr->ps->score = 0;

	chars = a_taxa[taxon_num];
	for (i = j = 0; i < a_nchars; i++)
	{
	    /* Ignore uninformative characters. */
	    if (a_chars_mask[i] == false)
	    {
		continue;
	    }

	    switch (chars[i])
	    {
		case 'N':
		case 'n':
		case 'X':
		case 'x':
		/* Treat gaps as uncertainty.  This isn't the only way to do
		 * things, and may need to be made configurable. */
		case '-':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0xf, j);
		    break;
		}
		case 'V':
		case 'v':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0xe, j);
		    break;
		}
		case 'H':
		case 'h':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0xd, j);
		    break;
		}
		case 'M':
		case 'm':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0xc, j);
		    break;
		}
		case 'D':
		case 'd':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0xb, j);
		    break;
		}
		case 'R':
		case 'r':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0xa, j);
		    break;
		}
		case 'W':
		case 'w':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0x9, j);
		    break;
		}
		case 'A':
		case 'a':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0x8, j);
		    break;
		}
		case 'B':
		case 'b':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0x7, j);
		    break;
		}
		case 'S':
		case 's':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0x6, j);
		    break;
		}
		case 'Y':
		case 'y':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0x5, j);
		    break;
		}
		case 'C':
		case 'c':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0x4, j);
		    break;
		}
		case 'K':
		case 'k':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0x3, j);
		    break;
		}
		case 'G':
		case 'g':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0x2, j);
		    break;
		}
		case 'T':
		case 't':
		{
		    tr_p_ps_char_set(a_tr, trr->ps, 0x1, j);
		    break;
		}
		default:
		{
		    cw_not_reached();
		}
	    }
	    j++;
	}
    }
}

static void
tr_p_mp_prepare_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring,
			char *a_taxa[], uint32_t a_ntaxa, uint32_t a_nchars,
			bool *a_chars_mask, uint32_t a_ninformative)
{
    cw_tr_ring_t ring;

    /* Prepare a_ring. */
    tr_p_mp_ring_prepare(a_tr, a_ring, a_taxa, a_ntaxa, a_nchars,
			 a_chars_mask, a_ninformative);

    /* Recurse into subtrees. */
    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	/* Prepare ring. */
	tr_p_mp_ring_prepare(a_tr, ring, a_taxa, a_ntaxa, a_nchars,
			     a_chars_mask, a_ninformative);

	/* Recurse. */
	tr_p_mp_prepare_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				a_taxa, a_ntaxa, a_nchars,
				a_chars_mask, a_ninformative);
    }
}

static void
tr_p_mp_ring_finish(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    cw_trr_t *trr;

    trr = &a_tr->trrs[a_ring];

    if (trr->ps != NULL)
    {
	tr_p_ps_delete(a_tr, trr->ps);
	trr->ps = NULL;
    }
}

static void
tr_p_mp_finish_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring)
{
    cw_tr_ring_t ring;

    /* Clean up a_ring. */
    tr_p_mp_ring_finish(a_tr, a_ring);

    /* Recurse into subtrees. */
    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
    {
	/* Clean up ring. */
	tr_p_mp_ring_finish(a_tr, ring);

	/* Recurse. */
	tr_p_mp_finish_recurse(a_tr, tr_p_ring_other_get(a_tr, ring));
    }
}

#ifdef CW_CPU_IA32
CW_P_INLINE void
tr_p_mp_ia32_pscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
		    cw_tr_ps_t *a_b)
{
    uint32_t curlimit, i, nbytes, ns;
    cw_trc_t *chars_p, *chars_a, *chars_b;

    /* Calculate node score. */
    ns = 0;

    /* Calculate partial Fitch parsimony scores for each character. */
    chars_p = a_p->chars;
    chars_a = a_a->chars;
    chars_b = a_b->chars;

    nbytes = (a_p->nchars >> 1);

    /* Initialize SSE2 registers. */
    {
	static const unsigned char low[] __attribute__ ((aligned (16))) =
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f"
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f";

	asm volatile (
	    /* Clear pns. */
	    "pxor %%xmm4, %%xmm4;"

	    /* Fill xmm5 with masks for the least significant four bits of each
	     * byte. */
	    "movdqa %[low], %%xmm5;"

	    /* Fill xmm6 with masks for the most significant four bits of each
	     * byte. */
	    "pcmpeqb %%xmm6, %%xmm6;"
	    "pxor %%xmm5, %%xmm6;"

	    /* Fill xmm7 with 16 1's. */
	    "pxor %%xmm7, %%xmm7;"
	    "pcmpeqb %%xmm0, %%xmm0;"
	    "psubb %%xmm0, %%xmm7;"
	    :
	    : [low] "m" (*low)
	    : "%xmm4", "%xmm5", "%xmm6", "%xmm7"
	    );
    }

    /* The inner loop can be run a maximum of 127 times before the partial node
     * score results (stored in %xmm4) are added to ns (otherwise, overflow
     * could occur).  Therefore, the outer loop calculates the upper bound for
     * the inner loop, thereby avoiding extra computation in the inner loop. */
    curlimit = 127 * 16;
    if (curlimit > nbytes)
    {
	curlimit = nbytes;
    }
    for (i = 0;;)
    {
	/* Use SSE2 to evaluate as many of the characters as possible.  This
	 * loop handles 32 characters per iteration. */
	for (; i < curlimit; i += 16)
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

		/**************************************************************/
		/* Most significant bits. */

		/* Create bitmasks according to whether the character state sets
		 * are empty.
		 *
		 * c = p ? 0x00 : 0xff;
		 * e = (c & d);
		 * s = c ? 0 : 1;
		 * p = (p | e);
		 * ns += s;
		 */
		"pxor %%xmm2, %%xmm2;"
		"movdqa %%xmm1, %%xmm3;"
		"pand %%xmm6, %%xmm3;" /* Mask out unwanted bits of p. */
		"pcmpeqb %%xmm3, %%xmm2;" /* xmm2 contains c. */
		"movdqa %%xmm0, %%xmm3;"
		"pand %%xmm6, %%xmm3;" /* Mask out unwanted bits of d. */
		"pand %%xmm2, %%xmm3;" /* xmm3 contains e. */
		"por %%xmm3, %%xmm1;" /* Update p. */
		"pand %%xmm7, %%xmm2;" /* xmm2 contains s. */
		"paddusb %%xmm2, %%xmm4;" /* Update ns (add s). */

		/**************************************************************/
		/* Least significant bits. */

		/* Create bitmasks according to whether the character state sets
		 * are empty.
		 *
		 * c = p ? 0x00 : 0xff;
		 * e = (c & d);
		 * s = c ? 0 : 1;
		 * p = (p | e);
		 * ns += s;
		 */
		"pxor %%xmm2, %%xmm2;"
		"movdqa %%xmm1, %%xmm3;"
		"pand %%xmm5, %%xmm3;" /* Mask out unwanted bits of p. */
		"pcmpeqb %%xmm3, %%xmm2;" /* xmm2 contains c. */
		"pand %%xmm5, %%xmm0;" /* Mask out unwanted bits of d. */
		"pand %%xmm2, %%xmm0;" /* xmm0 contains e. */
		"por %%xmm0, %%xmm1;" /* Update p. */
		"pand %%xmm7, %%xmm2;" /* xmm2 contains s. */
		"paddusb %%xmm2, %%xmm4;" /* Update ns (add s). */

		/* Store results.
		 *
		 * *chars_p = p;
		 */
		"movdqa %%xmm1, %[p];"
		: [p] "=m" (chars_p[i])
		: [a] "m" (chars_a[i]), [b] "m" (chars_b[i])
		: "memory"
		);
	}

	/* Update ns and reset pns. */
	{
	    uint32_t pns;

	    asm volatile (
		/* Sum the upper 8 bytes and lower 8 bytes separately (that's
		 * what psadbw does). */
		"pxor %%xmm0, %%xmm0;"
		"psadbw %%xmm0, %%xmm4;"

		/* Combine the results of psadbw. */
		"movdqa %%xmm4, %%xmm3;"
		"punpckhqdq %%xmm0, %%xmm4;"
		"paddq %%xmm3, %%xmm4;"

		/* Store the result. */
		"movd %%xmm4, %[pns];"
		"pxor %%xmm4, %%xmm4;"
		: [pns] "=r" (pns)
		:
		: "memory"
		);

	    ns += pns;
	}

	/* Break out of the loop if the bound for the inner loop was the maximum
	 * possible. */
	if (curlimit == nbytes)
	{
	    break;
	}
	/* Update the bound for the inner loop, taking care not to exceed the
	 * maximum possible bound. */
	curlimit += 127 * 16;
	if (curlimit > nbytes)
	{
	    curlimit = nbytes;
	}
    }

    a_p->score += ns;
}
#endif

static void
tr_p_mp_c_pscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
		 cw_tr_ps_t *a_b)
{
    uint32_t i, nbytes, ns, a, b, p, r;
    cw_trc_t *chars_p, *chars_a, *chars_b;

    /* Calculate node score. */
    ns = 0;

    /* Calculate preliminary Fitch parsimony scores for each character. */
    chars_p = a_p->chars;
    chars_a = a_a->chars;
    chars_b = a_b->chars;
    for (i = 0, nbytes = (a_p->nchars >> 1); i < nbytes; i++)
    {
	a = chars_a[i];
	b = chars_b[i];

	if ((p = ((a & b) & 0xf0U)) != 0)
	{
	    r = p;
	}
	else
	{
	    r = ((a | b) & 0xf0U);
	    ns++;
	}

	if ((p = ((a & b) & 0x0fU)) != 0)
	{
	    r |= p;
	}
	else
	{
	    r |= ((a | b) & 0x0fU);
	    ns++;
	}

	chars_p[i] = r;
    }

    a_p->score += ns;
}

/* Unconditionally calculate the partial score for a_p, using a_a and a_b as
 * children. */
CW_P_INLINE void
tr_p_mp_pscore(cw_tr_t *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b)
{
    /* Calculate sum of subtree scores. */
    a_p->score = a_a->score + a_b->score;

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

#ifdef CW_CPU_IA32
CW_P_INLINE uint32_t
tr_p_mp_ia32_fscore(cw_tr_t *a_tr, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b,
		    uint32_t a_maxscore)
{
    uint32_t retval, i, nbytes, pns;
    cw_trc_t *chars_a, *chars_b;

    /* Calculate sum of subtree scores. */
    retval = a_a->score + a_b->score;

    /* Calculate partial Fitch parsimony scores for each character. */
    chars_a = a_a->chars;
    chars_b = a_b->chars;

    /* Initialize SSE2 registers. */
    {
	static const unsigned char low[] __attribute__ ((aligned (16))) =
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f"
	    "\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f";

	asm volatile (
	    /* Clear pns. */
	    "pxor %%xmm4, %%xmm4;"

	    /* Fill xmm5 with masks for the least significant four bits of each
	     * byte. */
	    "movdqa %[low], %%xmm5;"

	    /* Fill xmm6 with masks for the most significant four bits of each
	     * byte. */
	    "pcmpeqb %%xmm6, %%xmm6;"
	    "pxor %%xmm5, %%xmm6;"

	    /* Fill xmm7 with 16 1's. */
	    "pxor %%xmm7, %%xmm7;"
	    "pcmpeqb %%xmm0, %%xmm0;"
	    "psubb %%xmm0, %%xmm7;"
	    :
	    : [low] "m" (*low)
	    : "%xmm4", "%xmm5", "%xmm6", "%xmm7"
	    );
    }

    /* Use SSE2 to evaluate as many of the characters as possible.  This loop
     * handles 32 characters per iteration. */
    for (i = 0, nbytes = (a_a->nchars >> 1); i < nbytes; i += 16)
    {
	asm volatile (
	    /* Read character data, and'ing and or'ing them together.
	     *
	     * a = *chars_a;
	     * b = *chars_b;
	     * p = a & b;
	     */
	    "movdqa %[a], %%xmm1;"
	    "pand %[b], %%xmm1;" /* xmm1 contains p. */

	    /**************************************************************/
	    /* Most significant bits. */

	    /* Create bitmasks according to whether the character state sets are
	     * empty.
	     *
	     * c = p ? 0x00 : 0xff;
	     * s = c ? 0 : 1;
	     * retval += s;
	     */
	    "pxor %%xmm2, %%xmm2;"
	    "movdqa %%xmm1, %%xmm3;"
	    "pand %%xmm6, %%xmm3;" /* Mask out unwanted bits of p. */
	    "pcmpeqb %%xmm3, %%xmm2;" /* xmm2 contains c. */
	    "pand %%xmm7, %%xmm2;" /* xmm2 contains s. */
	    "paddusb %%xmm2, %%xmm4;" /* Update retval (add s). */

	    /**************************************************************/
	    /* Least significant bits. */

	    /* Create bitmasks according to whether the character state sets are
	     * empty.
	     *
	     * c = p ? 0x00 : 0xff;
	     * s = c ? 0 : 1;
	     * retval += s;
	     */
	    "pxor %%xmm2, %%xmm2;"
	    "pand %%xmm5, %%xmm1;" /* Mask out unwanted bits of p. */
	    "pcmpeqb %%xmm1, %%xmm2;" /* xmm2 contains c. */
	    "pand %%xmm7, %%xmm2;" /* xmm2 contains s. */
	    "paddusb %%xmm2, %%xmm4;" /* Update retval (add s). */

	    /* Sum the upper 8 bytes and lower 8 bytes separately (there's no
	     * choice in the matter -- that's what psadbw does). */
	    "pxor %%xmm0, %%xmm0;"
	    "psadbw %%xmm0, %%xmm4;"

	    /* Combine the results of psadbw. */
	    "movdqa %%xmm4, %%xmm3;"
	    "punpckhqdq %%xmm0, %%xmm4;"
	    "paddq %%xmm3, %%xmm4;"

	    /* Store the result. */
	    "movd %%xmm4, %[pns];"
	    "pxor %%xmm4, %%xmm4;"

	    : [pns] "=r" (pns)
	    : [a] "m" (chars_a[i]), [b] "m" (chars_b[i])
	    : "memory"
	    );

	/* Update retval and terminate if the max score was exceeded. */
	retval += pns;
	if (retval > a_maxscore)
	{
	    retval = UINT_MAX;
	    break;
	}
    }

    return retval;
}
#endif

static uint32_t
tr_p_mp_c_fscore(cw_tr_t *a_tr, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b,
		 uint32_t a_maxscore)
{
    uint32_t retval, i, nbytes, a, b;
    cw_trc_t *chars_a, *chars_b;

    /* Calculate sum of subtree scores. */
    retval = a_a->score + a_b->score;

    /* Calculate partial Fitch parsimony scores for each character. */
    chars_a = a_a->chars;
    chars_b = a_b->chars;
    for (i = 0, nbytes = (a_a->nchars >> 1); i < nbytes; i++)
    {
	a = chars_a[i];
	b = chars_b[i];

	if (((a & b) & 0xf0U) == 0)
	{
	    retval++;
	    if (retval > a_maxscore)
	    {
		retval = UINT_MAX;
		break;
	    }
	}

	if (((a & b) & 0x0fU) == 0)
	{
	    retval++;
	    if (retval > a_maxscore)
	    {
		retval = UINT_MAX;
		break;
	    }
	}
    }

    return retval;
}

/* Unconditionally calculate the final score of a tree, using a_a and a_b as
 * children. */
CW_P_INLINE uint32_t
tr_p_mp_fscore(cw_tr_t *a_tr, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b,
	       uint32_t a_maxscore)
{
    uint32_t retval;

#ifdef CW_CPU_IA32
    if (modcrux_ia32_use_sse2)
    {
	retval = tr_p_mp_ia32_fscore(a_tr, a_a, a_b, a_maxscore);
    }
    else
#endif
    {
	retval = tr_p_mp_c_fscore(a_tr, a_a, a_b, a_maxscore);
    }

    return retval;
}

static cw_tr_ps_t *
tr_p_mp_score_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring, cw_tr_edge_t a_bisect)
{
    cw_tr_ps_t *retval;
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
	    retval = a_tr->trrs[a_ring].ps;
	    break;
	}
	case 2:
	{
	    /* This is a trifurcating node that is adjacent to the bisection.
	     * Return the child node's ps, since this node's ps is
	     * irrelevant. */
	    cw_assert(adjacent);

	    /* Get the ring element that connects to the other portion of the
	     * subtree on this side of the bisection. */
	    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
	    {
		if (tr_p_ring_edge_get(a_tr, ring) != a_bisect)
		{
		    retval
			= tr_p_mp_score_recurse(a_tr,
						tr_p_ring_other_get(a_tr, ring),
						a_bisect);
		    break;
		}
	    }
	    break;
	}
	case 3:
	{
	    if (adjacent == false)
	    {
		cw_tr_ps_t *ps_a, *ps_b;

		/* This is a normal trifurcating node.  This is the common case,
		 * and is handled separately from the code below for performance
		 * reasons. */

		/* Recursively calculate partial scores for the subtrees. */
		ring = qri_next(a_tr->trrs, a_ring, link);
		ps_a = tr_p_mp_score_recurse(a_tr,
					     tr_p_ring_other_get(a_tr, ring),
					     a_bisect);

		ring = qri_next(a_tr->trrs, ring, link);
		ps_b = tr_p_mp_score_recurse(a_tr,
					     tr_p_ring_other_get(a_tr, ring),
					     a_bisect);

		/* Calculate the partial score for this node. */
		retval = a_tr->trrs[a_ring].ps;
		tr_p_mp_pscore(a_tr, retval, ps_a, ps_b);

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

    return retval;
}

static cw_tr_ps_t *
tr_p_mp_views_recurse(cw_tr_t *a_tr, cw_tr_ring_t a_ring, cw_tr_ps_t *a_ps,
		      cw_tr_edge_t a_bisect)
{
    cw_tr_ps_t *retval = NULL;
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
	    /* This is a trifurcating node that is adjacent to the bisection.
	     * Pass the parent's ps when recursing, since this node's ps is
	     * irrelevant. */
	    cw_assert(adjacent);

	    /* Get the ring element that connects to the other portion of the
	     * subtree on this side of the bisection. */
	    qri_others_foreach(ring, a_tr->trrs, a_ring, link)
	    {
		if (tr_p_ring_edge_get(a_tr, ring) != a_bisect)
		{
		    /* Recurse. */
		    tr_p_mp_views_recurse(a_tr,
					  tr_p_ring_other_get(a_tr, ring),
					  a_ps,
					  a_bisect);
		}
		else
		{
		    retval = a_tr->trrs[ring].ps;
		}
	    }
	    break;
	}
	case 3:
	{
	    if (adjacent == false)
	    {
		cw_tr_ring_t ring_a, ring_b;
		cw_tr_ring_t ring_a_other, ring_b_other;
		cw_tr_ps_t *ps_a, *ps_b;
		cw_tr_ps_t *ps_a_other, *ps_b_other;

		/* This is a normal trifurcating node.  This is the common case,
		 * and is handled separately from the code below for performance
		 * reasons. */

		/* Get all variables that are necessary for view calculation and
		 * recursion. */
		ring_a = qri_next(a_tr->trrs, a_ring, link);
		ps_a = a_tr->trrs[ring_a].ps;
		ring_a_other = tr_p_ring_other_get(a_tr, ring_a);
		ps_a_other = a_tr->trrs[ring_a_other].ps;

		ring_b = qri_next(a_tr->trrs, ring_a, link);
		ps_b = a_tr->trrs[ring_b].ps;
		ring_b_other = tr_p_ring_other_get(a_tr, ring_b);
		ps_b_other = a_tr->trrs[ring_b_other].ps;

		/* Calculate views and edges, and recurse. */
		tr_p_mp_pscore(a_tr, ps_a, a_ps, ps_b_other);
		tr_p_mp_pscore(a_tr, ps_b, a_ps, ps_a_other);

		tr_p_mp_views_recurse(a_tr, ring_a_other, ps_a, a_bisect);
		tr_p_mp_pscore(a_tr, ps_a, ps_a, ps_a_other);
		a_tr->tres[tr_p_ring_edge_get(a_tr, ring_a)].ps = ps_a;

		tr_p_mp_views_recurse(a_tr, ring_b_other, ps_b, a_bisect);
		tr_p_mp_pscore(a_tr, ps_b, ps_b, ps_b_other);
		a_tr->tres[tr_p_ring_edge_get(a_tr, ring_b)].ps = ps_b;

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

    return retval;
}

/* Calculate the partial score for each edge in a_edges.  a_edges[0] must either
 * be CW_TR_EDGE_NONE, or the edge connected to the node that is in turn
 * connected to the bisection edge. */
CW_P_INLINE bool
tr_p_bisection_edge_list_mp(cw_tr_t *a_tr, cw_tr_edge_t *a_edges,
			    uint32_t a_nedges, cw_tr_edge_t a_bisect,
			    uint32_t a_maxscore)
{
    bool retval;

    if (a_edges[0] != CW_TR_EDGE_NONE)
    {
	cw_tr_ring_t ring_a, ring_b;
	cw_tr_ps_t *ps_a, *ps_b, *ps_x, *ps_y;

	ring_a = tr_p_edge_ring_get(a_tr, a_edges[0], 0);
	ring_b = tr_p_edge_ring_get(a_tr, a_edges[0], 1);

	/* Recursively (post-order traversal) calculate the partial score at
	 * each node, as viewed from the first edge in a_edges.  This leaves one
	 * valid view at each node, which then makes it possible to calculate
	 * the rest of the views during a pre-order traversal of the tree. */
	ps_a = tr_p_mp_score_recurse(a_tr, ring_a, a_bisect);
	ps_b = tr_p_mp_score_recurse(a_tr, ring_b, a_bisect);

	/* The first edge must be calculated using ps_a and ps_b as children,
	 * rather than using the ps's at the ends of the edge.  This is because
	 * one of the connected nodes is in turn connected to the bisection
	 * edge, which means that the node does not have a useful ps.  The first
	 * edge is the only one for which this is an issue, so it is handled
	 * here. */
	// XXX This is extra work that need not be done, if things are
	// restructured a bit.
	if (tr_p_mp_fscore(a_tr, ps_a, ps_b, a_maxscore) > a_maxscore)
	{
	    /* Don't bother calculating other views or edge states, since this
	     * subtree exceeds the maximum score that's of interest. */
	    retval = true;
	    goto RETURN;
	}

	/* Perform the pre-order traversal, calculating the remaining views that
	 * were not calculated by the above post-order traversal, as well as
	 * calculating the state sets for the edges along the way.  Take care to
	 * pass the appropriate ps's. */
	ps_x = tr_p_mp_views_recurse(a_tr, ring_a, ps_b, a_bisect);
	ps_y = tr_p_mp_views_recurse(a_tr, ring_b, ps_a, a_bisect);

	if (ps_x == NULL)
	{
	    ps_x = ps_y;
	}
	cw_check_ptr(ps_x);
	tr_p_mp_pscore(a_tr, ps_x, ps_a, ps_b);
	a_tr->tres[a_edges[0]].ps = ps_x;

#ifdef CW_DBG
	/* Validate per-edge partial scores. */
	{
	    uint32_t i;

	    for (i = 1; i < a_nedges; i++)
	    {
		/* All edge partial scores should have the same value, since the
		 * location of the root is irrelevant to the score. */
		if (a_tr->tres[a_edges[i]].ps->score
		    != a_tr->tres[a_edges[0]].ps->score)
		{
		    fprintf(stderr,
			    "%s:%d:%s(): Expected %u, got %u\n",
			    __FILE__, __LINE__, __func__,
			    a_tr->tres[a_edges[0]].ps->score,
			    a_tr->tres[a_edges[i]].ps->score);
		    abort();
		}
	    }
	}
#endif
    }

    retval = false;
    RETURN:
    return retval;
}

/* Hold a tree.  If a_max_held is exceeded, the tree is not held.  This
 * introduces a bias in which trees are held.  There exist algorithms for making
 * this an unbiased process, but there is no need for that functionality at the
 * moment. */
CW_P_INLINE bool
tr_p_hold(cw_tr_t *a_tr, uint32_t a_max_hold, uint32_t a_neighbor,
	  uint32_t a_score)
{
    bool retval;

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

	retval = false;
    }
    else
    {
	retval = true;
    }

    return retval;
}

/* Calculate the Fitch parsimony scores for all TBR neighbors of a_tr, and hold
 * results according to the function parameters. */
CW_P_INLINE void
tr_p_tbr_neighbors_mp(cw_tr_t *a_tr, uint32_t a_max_hold,
		      uint32_t a_maxscore, cw_tr_hold_how_t a_how)
{
    uint32_t neighbor, i, j, k, curmax, score;
    cw_tr_edge_t bisect, edge_a, edge_b;
    cw_tr_node_t node_a, node_b;
    cw_tr_ps_t *ps_a, *ps_b;

    cw_dassert(tr_p_validate(a_tr));

    curmax = a_maxscore;

    /* Set up tree holding data structures. */
    a_tr->nheld = 0;

    /* Iteratively (logically) bisect at each edge in the tree. */
    neighbor = 0;
    for (i = 0; i < a_tr->nedges; i++)
    {
	bisect = a_tr->trti[i];

	/* Determine which edges are in each subtree. */
	tr_p_bedges_gen(a_tr, bisect, &node_a, &node_b);

	/* Calculate the partial score for each edge in the edge lists.  Don't
	 * bother scoring the trees if either subtree exceeds the max score. */
	if (tr_p_bisection_edge_list_mp(a_tr, a_tr->bedges,
					a_tr->nbedges_a, bisect, curmax)
	    || tr_p_bisection_edge_list_mp(a_tr, &a_tr->bedges[a_tr->nbedges_a],
					   a_tr->nbedges_b, bisect, curmax))
	{
	    neighbor += (a_tr->trt[i + 1].offset - a_tr->trt[i].offset);
	    continue;
	}

	/* Iteratively (logically) reconnect every legitimate pairing of edges
	 * between the two subtrees and calculate final parsimony scores. */
	for (j = 0; j < a_tr->nbedges_a; j++)
	{
	    edge_a = a_tr->bedges[j];
	    if (edge_a != CW_TR_EDGE_NONE)
	    {
		ps_a = a_tr->tres[edge_a].ps;
	    }
	    else
	    {
		ps_a = a_tr->trrs[qli_first(&a_tr->trns[node_a].rings)].ps;
	    }

	    for (k = 0; k < a_tr->nbedges_b; k++)
	    {
		/* Skip this iteration if the reconnection would result in
		 * reversing the bisection. */
		if (j == 0 && k == 0)
		{
		    continue;
		}

		edge_b = a_tr->bedges[a_tr->nbedges_a + k];
		if (edge_b != CW_TR_EDGE_NONE)
		{
		    ps_b = a_tr->tres[edge_b].ps;
		}
		else
		{
		    ps_b = a_tr->trrs[qli_first(&a_tr->trns[node_b].rings)].ps;
		}

		/* Calculate the final parsimony score for this reconnection. */
		score = tr_p_mp_fscore(a_tr, ps_a, ps_b, curmax);

		/* Hold the tree, if appropriate. */
		switch (a_how)
		{
		    case TR_HOLD_BEST:
		    {
			if (score < curmax)
			{
			    a_tr->nheld = 0;
			}

			if (score <= curmax || a_tr->nheld == 0)
			{
			    /* No trees held, or this tree is as good as those
			     * currently held. */
			    if (tr_p_hold(a_tr, a_max_hold, neighbor, score))
			    {
				/* No more room for trees. */
				curmax = score - 1;
			    }
			    else
			    {
				curmax = score;
			    }
			}
			break;
		    }
		    case TR_HOLD_BETTER:
		    {
			if (score <= curmax)
			{
			    /* No trees held, or this (neighboring) tree is
			     * better than the tree whose neighbors are being
			     * evaluated. */
			    tr_p_hold(a_tr, a_max_hold, neighbor, score);
			    curmax = score - 1;
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

    retval = (cw_tr_t *) cw_opaque_alloc(mema_alloc_get(a_mema),
					 mema_arg_get(a_mema),
					 sizeof(cw_tr_t));

    tr_p_new(retval, a_mema);

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

    retval->base = a_tr->base;
    retval->ntaxa = a_tr->ntaxa;
    retval->nedges = a_tr->nedges;

    /*
     * Copy trn-related data structures.
     */

    /* Allocate trns the same size as a_tr's, then copy. */
    if (a_tr->trns != NULL)
    {
	retval->trns
	    = (cw_trn_t *) cw_opaque_alloc(alloc, arg,
					   sizeof(cw_trn_t) * a_tr->ntrns);
	memcpy(retval->trns, a_tr->trns, sizeof(cw_trn_t) * a_tr->ntrns);
	retval->ntrns = a_tr->ntrns;

	/* Clean up the copied trn's. */
	for (i = 0; i < retval->ntrns; i++)
	{
	    a_tr->trns[i].u.aux = NULL;
	}

	/* The spare trns list is the same as for a_tr. */
	retval->sparetrns = a_tr->sparetrns;
    }

    /*
     * Copy tre-related data structures.
     */

    if (a_tr->tres != NULL)
    {
	/* Allocate tres the same size as a_tr's, then copy. */
	retval->tres
	    = (cw_tre_t *) cw_opaque_alloc(alloc, arg,
					   sizeof(cw_tre_t) * a_tr->ntres);
	memcpy(retval->tres, a_tr->tres, sizeof(cw_tre_t) * a_tr->ntres);
	retval->ntres = a_tr->ntres;

	/* Clean up the copied tre's. */
	for (i = 0; i < retval->ntres; i++)
	{
	    a_tr->tres[i].u.aux = NULL;
	    a_tr->tres[i].ps = NULL;
	}

	/* The spare tres list is the same as for a_tr. */
	retval->sparetres = a_tr->sparetres;

	/* Alocate trrs the same size as a_tr's, then copy. */
	retval->trrs
	    = (cw_trr_t *) cw_opaque_alloc(alloc, arg,
					   sizeof(cw_trr_t) * a_tr->ntres * 2);
	memcpy(retval->trrs, a_tr->trrs, sizeof(cw_trr_t *) * a_tr->ntres * 2);

	for (i = 0; i < retval->ntres * 2; i++)
	{
	    a_tr->trrs[i].ps = NULL;
	}
    }

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

    /* This assumes that all nodes are deallocated before tr_delete() is
     * called. */
    if (a_tr->trns != NULL)
    {
	cw_opaque_dealloc(dealloc, arg, a_tr->trns,
			  sizeof(cw_trn_t) * a_tr->ntrns);
    }

    if (a_tr->trt != NULL)
    {
	cw_opaque_dealloc(dealloc, arg, a_tr->trt,
			  sizeof(cw_trt_t) * (a_tr->nedges + 1));
	cw_opaque_dealloc(dealloc, arg, a_tr->trti,
			  sizeof(cw_tr_edge_t) * a_tr->nedges);
    }

    if (a_tr->bedges != NULL)
    {
	cw_opaque_dealloc(mema_dealloc_get(a_tr->mema),
			  mema_arg_get(a_tr->mema),
			  a_tr->bedges, sizeof(cw_tr_edge_t) * a_tr->nedges);
    }

    /* This assumes that all edges are deallocated before tr_delete() is
     * called. */
    if (a_tr->tres != NULL)
    {
	cw_opaque_dealloc(dealloc, arg, a_tr->tres,
			  sizeof(cw_tre_t) * a_tr->ntres);
	cw_opaque_dealloc(dealloc, arg, a_tr->trrs,
			  sizeof(cw_trr_t) * a_tr->ntres * 2);
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
    /* Update internal state, so that ntaxa and nedges are correct. */
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    if (a_tr->base != CW_TR_NODE_NONE)
    {
	uint32_t ntaxa, nedges;
	cw_tr_ring_t ring;

	/* Set base to be the lowest-numbered taxon. */
	ntaxa = 0;
	nedges = 0;
	a_tr->base = tr_p_lowest(a_tr, a_tr->base, &ntaxa, &nedges);

	/* Get base's ring. */
	ring = qli_first(&a_tr->trns[a_tr->base].rings);
	if (ring != CW_TR_RING_NONE)
	{
	    /* Canonize the tree. */
	    tr_p_canonize(a_tr, tr_p_ring_other_get(a_tr, ring));
	}
    }

    /* Re-update internal state. */
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));
}

void
tr_tbr(cw_tr_t *a_tr, cw_tr_edge_t a_bisect, cw_tr_edge_t a_reconnect_a,
       cw_tr_edge_t a_reconnect_b,
       cw_tr_edge_t (*a_edge_alloc_callback)(cw_tr_t *, void *),
       cw_tr_node_t (*a_node_alloc_callback)(cw_tr_t *, void *),
       void *a_arg)
{
    cw_tr_node_t node_a, node_b, nodes[4];
    cw_tr_edge_t tedges[3];
    uint32_t ntedges = 0;
    cw_tr_node_t tnodes[2];
    uint32_t ntnodes = 0;

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    /* Get the nodes to either side of the edge where the bisection will be
     * done. */
    node_a = tr_edge_node_get(a_tr, a_bisect, 0);
    node_b = tr_edge_node_get(a_tr, a_bisect, 1);

    /* Bisect.  a_bisect will be used below for reconnection. */
    tr_edge_detach(a_tr, a_bisect);
    
    /* For nodes_[ab], extract the node if it has only two neighbors.
     *
     * nodes_[0..1] are CW_TR_NODE_NONE, unless they refer to the only node in a
     * subtree. */
    nodes[0] = tr_p_tbr_node_extract(a_tr, node_a, a_reconnect_a, a_reconnect_b,
				     tedges, &ntedges, tnodes, &ntnodes);
    nodes[1] = tr_p_tbr_node_extract(a_tr, node_b, a_reconnect_a, a_reconnect_b,
				     tedges, &ntedges, tnodes, &ntnodes);

    /* For each reconnection edge, splice a node into the edge (if the subtree
     * has more than one node).
     *
     * nodes[2..3] are set to CW_TR_NODE_NONE if no reconnection edge is
     * specified. */
    if (a_reconnect_a != CW_TR_EDGE_NONE)
    {
	nodes[2] = tr_p_tbr_node_splice(a_tr, a_reconnect_a,
					tedges, &ntedges, tnodes, &ntnodes,
					a_edge_alloc_callback,
					a_node_alloc_callback,
					a_arg);
    }
    else
    {
	nodes[2] = CW_TR_NODE_NONE;
    }

    if (a_reconnect_b != CW_TR_EDGE_NONE)
    {
	nodes[3] = tr_p_tbr_node_splice(a_tr, a_reconnect_b,
					tedges, &ntedges, tnodes, &ntnodes,
					a_edge_alloc_callback,
					a_node_alloc_callback,
					a_arg);
    }
    else
    {
	nodes[3] = CW_TR_NODE_NONE;
    }

    /* If either subtree has only a single node, special care must be taken
     * during reconnection to re-associate the proper end of the bisection edge
     * with the single node.  This is because character state information for
     * leaf nodes is actually stored in the rings that attach them to the tree,
     * and breaking this association would require re-initializing the character
     * state set vectors. */
    if (nodes[0] != CW_TR_NODE_NONE)
    {
	/* nodes[0] (same as node_a) is a single node. */
	cw_assert(nodes[0] == node_a);
	cw_assert(nodes[1] == CW_TR_NODE_NONE);

	if (nodes[2] == CW_TR_EDGE_NONE)
	{
	    nodes[2] = nodes[3];
	}
	tr_edge_attach(a_tr, a_bisect, nodes[0], nodes[2]);
    }
    else if (nodes[1] != CW_TR_NODE_NONE)
    {
	/* nodes[1] (same as node_b) is a single node. */
	cw_assert(nodes[1] == node_b);

	if (nodes[2] == CW_TR_EDGE_NONE)
	{
	    nodes[2] = nodes[3];
	}
	tr_edge_attach(a_tr, a_bisect, nodes[2], nodes[1]);
    }
    else
    {
	/* Bisection was not done adjacent to a leaf node.  Attach the two
	 * spliced-in nodes. */
	tr_edge_attach(a_tr, a_bisect, nodes[2], nodes[3]);
    }

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
		    cw_tr_edge_t *r_bisect, cw_tr_edge_t *r_reconnect_a,
		    cw_tr_edge_t *r_reconnect_b)
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
    tr_p_bedges_gen(a_tr, trt->bisect_edge, NULL, NULL);

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
tr_mp_prepare(cw_tr_t *a_tr, bool a_uninformative_eliminate,
	      char *a_taxa[], uint32_t a_ntaxa, uint32_t a_nchars)
{
    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));
    if (a_tr->base != CW_TR_NODE_NONE)
    {
	cw_trn_t *trn;
	cw_tr_ring_t ring;
	uint32_t i, ninformative;
	bool chars_mask[a_nchars];

	if (a_uninformative_eliminate)
	{
	    uint32_t codes[15];
	    uint32_t j, k, x, y;

	    /* Preprocess the character data.  Eliminate uninformative
	     * characters, but keep track of their contribution to the parsimony
	     * score, were they to be left in. */
	    ninformative = 0;
	    for (i = 0; i < a_nchars; i++)
	    {
		for (k = 0; k < 15; k++)
		{
		    codes[k] = 0;
		}

		for (j = 0; j < a_ntaxa; j++)
		{
		    switch (a_taxa[j][i])
		    {
			case 'N':
			case 'n':
			case 'X':
			case 'x':
			/* Treat gaps as uncertainty.  This isn't the only way
			 * to do things, and may need to be made
			 * configurable. */
			case '-':
			{
			    break;
			}
			case 'V':
			case 'v':
			{
			    codes[14]++;
			    break;
			}
			case 'H':
			case 'h':
			{
			    codes[13]++;
			    break;
			}
			case 'M':
			case 'm':
			{
			    codes[12]++;
			    break;
			}
			case 'D':
			case 'd':
			{
			    codes[11]++;
			    break;
			}
			case 'R':
			case 'r':
			{
			    codes[10]++;
			    break;
			}
			case 'W':
			case 'w':
			{
			    codes[9]++;
			    break;
			}
			case 'A':
			case 'a':
			{
			    codes[8]++;
			    break;
			}
			case 'B':
			case 'b':
			{
			    codes[7]++;
			    break;
			}
			case 'S':
			case 's':
			{
			    codes[6]++;
			    break;
			}
			case 'Y':
			case 'y':
			{
			    codes[5]++;
			    break;
			}
			case 'C':
			case 'c':
			{
			    codes[4]++;
			    break;
			}
			case 'K':
			case 'k':
			{
			    codes[3]++;
			    break;
			}
			case 'G':
			case 'g':
			{
			    codes[2]++;
			    break;
			}
			case 'T':
			case 't':
			{
			    codes[1]++;
			    break;
			}
			default:
			{
			    cw_not_reached();
			}
		    }
		}

		/* Count the number of states in which two or more taxa
		 * exist. */
		chars_mask[i] = false;
		for (x = 1; x < 15; x++)
		{
		    for (y = 1; y < 15; y++)
		    {
			if ((x & y) == 0 && codes[x] >= 2 && codes[y] >= 2)
			{
			    if (chars_mask[i] == false)
			    {
				chars_mask[i] = true;
				ninformative++;
			    }
			}
		    }
		}
	    }
	}
	else
	{
	    /* Use all characters, regardless of whether they are
	     * informative. */
	    ninformative = a_nchars;

	    for (i = 0; i < a_nchars; i++)
	    {
		chars_mask[i] = true;
	    }
	}

	/* Prepare the tree. */
	trn = &a_tr->trns[a_tr->base];
	qli_foreach(ring, &trn->rings, a_tr->trrs, link)
	{
	    /* Prepare ring. */
	    tr_p_mp_ring_prepare(a_tr, ring, a_taxa, a_ntaxa, a_nchars,
				 chars_mask, ninformative);

	    /* Recurse. */
	    tr_p_mp_prepare_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				    a_taxa, a_ntaxa, a_nchars,
				    chars_mask, ninformative);
	}
    }
}

void
tr_mp_finish(cw_tr_t *a_tr)
{

    tr_p_update(a_tr);
    cw_dassert(tr_p_validate(a_tr));

    if (a_tr->base != CW_TR_NODE_NONE)
    {
	cw_trn_t *trn;
	cw_tr_ring_t ring;

	/* Clean up the tree. */
	trn = &a_tr->trns[a_tr->base];
	qli_foreach(ring, &trn->rings, a_tr->trrs, link)
	{
	    /* Clean up ring. */
	    tr_p_mp_ring_finish(a_tr, ring);

	    /* Recurse. */
	    tr_p_mp_finish_recurse(a_tr, tr_p_ring_other_get(a_tr, ring));
	}
    }
}

uint32_t
tr_mp_score(cw_tr_t *a_tr)
{
    uint32_t retval;
    cw_tr_ring_t ring;

    cw_dassert(tr_p_validate(a_tr));

    if (a_tr->base != CW_TR_NODE_NONE
	&& (ring = qli_first(&a_tr->trns[a_tr->base].rings)) != CW_TR_RING_NONE)
    {
	cw_tr_edge_t edge;
	cw_tr_ps_t *ps_a, *ps_b;

	edge = tr_p_ring_edge_get(a_tr, ring);

	/* Calculate partial scores for the subtrees on each end of edge. */
	ps_a = tr_p_mp_score_recurse(a_tr,
				     tr_p_edge_ring_get(a_tr, edge, 0),
				     CW_TR_EDGE_NONE);
	ps_b = tr_p_mp_score_recurse(a_tr,
				     tr_p_edge_ring_get(a_tr, edge, 1),
				     CW_TR_EDGE_NONE);

	/* Calculate the final score. */
	retval = tr_p_mp_fscore(a_tr, ps_a, ps_b, UINT_MAX);
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
    uint32_t score;

    cw_dassert(tr_p_validate(a_tr));

    score = tr_mp_score(a_tr);
    tr_p_tbr_neighbors_mp(a_tr, a_max_hold,
			  score > 0 ? score - 1 : 0,
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