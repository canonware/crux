/******************************************************************************
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
 * indices to be used for links between nodes and edges.
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

#include "../include/_cruxmodule.h"

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
#ifdef CxmDebug
    uint32_t magic;
#define CW_TRN_MAGIC 0x63329478
#endif

    union
    {
	/* Auxiliary opaque data pointer. */
	void *aux;

	/* Spares linkage. */
	CxtTrNode link;
    } u;

    /* If CxmTrNodeTaxonNone, then the node is not a leaf node. */
    uint32_t taxon_num;

    /* Ring of trr's, which are associated with tre's. */
    CxmQliHead rings;
};

/* Tree node edge ring element. */
struct cw_trr_s
{
    /* Ring linkage. */
    CxmQri link;

    /* Node whose ring this trr is a part of. */
    CxtTrNode node;

    /* Used for Fitch parsimony scoring. */
    cw_tr_ps_t *ps;
};

/* Tree edge information. */
struct cw_tre_s
{
#ifdef CxmDebug
    uint32_t magic;
#define CW_TRE_MAGIC 0xa683fa07
#endif

    union
    {
	/* Auxiliary opaque data pointer. */
	void *aux;

	/* Spares linkage. */
	CxtTrEdge link;
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
    CxtTrEdge bisect_edge;
};

/* Held neighbor tree. */
struct cw_trh_s
{
    /* Neighbor index for the tree.  This can be passed to CxTrTbrNeighborGet()
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

struct CxsTr
{
#ifdef CxmDebug
    uint32_t magic;
#define CW_TR_MAGIC 0x39886394
#endif

    /* Auxiliary opaque data pointer. */
    void *aux;

    /* True if this tree has been modified since the internal state (ntaxa,
     * nedges, bedges, trt) was updated, false otherwise. */
    bool modified;

    /* Base of the tree (may or may not be set). */
    CxtTrNode base;

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
    CxtTrEdge *bedges;
    uint32_t nbedges_a;
    uint32_t nbedges_b;

    /* trt is an array of elements that store per-edge information that is used
     * for TBR-related functions.  There is one more element in trt than there
     * are edges in the tree.  This is critical to the way binary searching on
     * the array is done, and it also makes it easy to get the total number of
     * TBR neighbors this tree has (trt[nedges].offset).
     *
     * Only the first trtused elements are valid, since not all bisection edges
     * necessarily result in neighbors. */
    cw_trt_t *trt;
    uint32_t trtused;

    /* Pointer to an array of trn's.  ntrns is the total number of trn's, not
     * all of which are necessarily in use.
     *
     * sparetrns is the index of the first spare trn in the spares stack. */
    cw_trn_t *trns;
    uint32_t ntrns;
    CxtTrNode sparetrns;

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
    CxtTrEdge sparetres;
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

CxmpInline void
tr_p_ring_init(CxtTr *a_tr, cw_tr_ring_t a_ring)
{
    cw_trr_t *trr;

    cxmAssert((a_ring >> 1) < a_tr->ntres);

    trr = &a_tr->trrs[a_ring];

    CxmQriNew(a_tr->trrs, a_ring, link);
    trr->node = CxmTrNodeNone;
    trr->ps = NULL;
}

CxmpInline CxtTrEdge
tr_p_ring_edge_get(CxtTr *a_tr, cw_tr_ring_t a_ring)
{
    cxmAssert((a_ring >> 1) < a_tr->ntres);

    return (a_ring >> 1);
}

CxmpInline CxtTrNode
tr_p_ring_node_get(CxtTr *a_tr, cw_tr_ring_t a_ring)
{
    cxmAssert((a_ring >> 1) < a_tr->ntres);

    return a_tr->trrs[a_ring].node;
}

CxmpInline cw_tr_ring_t
tr_p_ring_other_get(CxtTr *a_tr, cw_tr_ring_t a_ring)
{
    cxmAssert((a_ring >> 1) < a_tr->ntres);

    return (a_ring ^ 1);
}

/******************************************************************************/

/* Validation functions. */

#ifdef CxmDebug
/* Validate a ring. */
static bool
tr_p_ring_validate(CxtTr *a_tr, cw_tr_ring_t a_ring)
{
    cw_trr_t *trr;

    cxmCheckPtr(a_tr);
    cxmAssert(a_tr->magic == CW_TR_MAGIC);
    cxmAssert((a_ring >> 1) < a_tr->ntres);

    trr = &a_tr->trrs[a_ring];

    cxmAssert(trr->node < a_tr->ntrns || trr->node == CxmTrNodeNone);

    return true;
}

/* Validate an edge. */
static bool
tr_p_edge_validate(CxtTr *a_tr, CxtTrEdge a_edge)
{
    cw_tre_t *tre;

    cxmCheckPtr(a_tr);
    cxmAssert(a_tr->magic == CW_TR_MAGIC);
    cxmAssert(a_edge < a_tr->ntres);

    tre = &a_tr->tres[a_edge];

    cxmAssert(tre->magic == CW_TRE_MAGIC);

    tr_p_ring_validate(a_tr, (a_edge << 1));
    tr_p_ring_validate(a_tr, (a_edge << 1) + 1);

    return true;
}

/* Validate a node. */
static bool
tr_p_node_validate(CxtTr *a_tr, CxtTrNode a_node)
{
    cw_trn_t *trn;
    cw_tr_ring_t ring;
    uint32_t nneighbors;

    cxmCheckPtr(a_tr);
    cxmAssert(a_tr->magic == CW_TR_MAGIC);
    cxmAssert(a_node < a_tr->ntrns);

    trn = &a_tr->trns[a_node];

    cxmAssert(trn->magic == CW_TRN_MAGIC);

    nneighbors = 0;
    CxmQliForeach(ring, &trn->rings, a_tr->trrs, link)
    {
	/* Validate edge. */
	tr_p_edge_validate(a_tr, tr_p_ring_edge_get(a_tr, ring));

	nneighbors++;
    }

    if (trn->taxon_num != CxmTrNodeTaxonNone)
    {
	/* Only leaf nodes can have taxon numbers.  Leaf nodes have at most
	 * 1 neighbor. */
	cxmAssert(nneighbors <= 1);
    }

    return true;
}
#endif

/******************************************************************************/

/* tr_ps. */

CxmpInline cw_tr_ps_t *
tr_p_ps_new(CxtTr *a_tr)
{
    cw_tr_ps_t *retval;

    retval = (cw_tr_ps_t *) CxmMalloc(sizeof(cw_tr_ps_t));

    retval->parent = NULL;
    retval->chars = NULL;

    return retval;
}

CxmpInline void
tr_p_ps_delete(CxtTr *a_tr, cw_tr_ps_t *a_ps)
{
    if (a_ps->chars != NULL)
    {
	CxmFree(a_ps->achars);
    }

    CxmFree(a_ps);
}

#if (0) /* Unused (so far). */
CxmpInline cw_trc_t
tr_p_ps_char_get(CxtTr *a_tr, cw_tr_ps_t *a_ps, uint32_t a_offset)
{
    cw_trc_t retval;

    cxmCheckPtr(a_ps->chars);
    cxmAssert(a_offset < a_ps->nchars);

    retval = a_ps->chars[a_offset >> 1];
    retval >>= ((a_offset & 1) * 4);

    return retval;
}
#endif

CxmpInline void
tr_p_ps_char_set(CxtTr *a_tr, cw_tr_ps_t *a_ps, cw_trc_t a_char,
		 uint32_t a_offset)
{
    cxmCheckPtr(a_ps->chars);
    cxmAssert((a_char & 0xfU) == a_char);
    cxmAssert(a_offset
	      < a_ps->nchars + ((32 - (a_ps->nchars & 0x1fU)) & 0x1fU));

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

CxmpInline void
tr_p_ps_prepare(CxtTr *a_tr, cw_tr_ps_t *a_ps, uint32_t a_nchars)
{
    /* Clean up old character vector if it isn't the right size for a_nchars
     * characters. */
    if (a_ps->chars != NULL && a_ps->nchars != a_nchars)
    {
	CxmFree(a_ps->achars);
	a_ps->chars = NULL;
    }

    /* Allocate character vector if necessary. */
    if (a_ps->chars == NULL)
    {
	if (a_nchars != 0)
	{
	    uint32_t npad, i;

	    /* Calculate the number of pad bytes to append, such that the total
	     * number of bytes is a multiple of 16 (total number of taxonomical
	     * characters is a multiple of 32). */
	    npad = (32 - (a_nchars & 0x1fU)) & 0x1fU;
	    cxmAssert(((a_nchars + npad) & 0x1fU) == 0);

	    /* Tack on 8 bytes; all modern systems provide at least 8 byte
	     * alignment. */
	    a_ps->achars = (cw_trc_t *) CxmMalloc(sizeof(cw_trc_t)
						  * (((a_nchars + npad) >> 1))
						  + 8);

	    /* Make sure that chars is 16 byte-aligned.  Assume that achars is
	     * at least 8 byte-aligned. */
	    a_ps->chars = &a_ps->achars[((unsigned) a_ps->achars) & 0xfU];

	    a_ps->nchars = a_nchars + npad;

	    /* Set all pad characters to {ACGT}.  This allows the pad characters
	     * to be calculated along with the actual characters, without
	     * affecting the score. */
	    for (i = a_nchars; i < a_nchars + npad; i++)
	    {
		tr_p_ps_char_set(a_tr, a_ps, 0xfU, i);
	    }
	}
	else
	{
	    a_ps->achars = NULL;
	    a_ps->chars = NULL;
	    a_ps->nchars = 0;
	}
    }
}

/******************************************************************************/

/* tr_edge. */

CxmpInline void
tr_p_CxEdgeInit(CxtTr *a_tr, CxtTrEdge a_edge)
{
    cw_tre_t *tre;

    tre = &a_tr->tres[a_edge];

    tre->u.aux = NULL;
    tr_p_ring_init(a_tr, (a_edge << 1));
    tr_p_ring_init(a_tr, (a_edge << 1) + 1);
    tre->length = 0.0;
    tre->ps = NULL;

#ifdef CxmDebug
    tre->magic = CW_TRE_MAGIC;
#endif
}

CxmpInline CxtTrEdge
tr_p_edge_alloc(CxtTr *a_tr)
{
    CxtTrEdge retval;

    if (a_tr->sparetres == CxmTrEdgeNone)
    {
	uint32_t i, nspares;

	if (a_tr->tres == NULL)
	{
	    a_tr->tres = (cw_tre_t *) CxmMalloc(sizeof(cw_tre_t));
	    cxmAssert(a_tr->trrs == NULL);
	    a_tr->trrs = (cw_trr_t *) CxmMalloc(sizeof(cw_trr_t) * 2);
	    nspares = 1;
	    a_tr->ntres = 1;
	}
	else
	{
	    a_tr->tres = (cw_tre_t *) CxmRealloc(a_tr->tres,
						 sizeof(cw_tre_t)
						 * a_tr->ntres * 2);
	    cxmCheckPtr(a_tr->trrs);
	    a_tr->trrs = (cw_trr_t *) CxmRealloc(a_tr->trrs,
						 sizeof(cw_trr_t)
						 * a_tr->ntres * 4);
	    nspares = a_tr->ntres;
	    a_tr->ntres *= 2;
	}

	/* Initialize last spare. */
	a_tr->sparetres = a_tr->ntres - 1;
	a_tr->tres[a_tr->sparetres].u.link = CxmTrEdgeNone;

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
    tr_p_CxEdgeInit(a_tr, retval);

    return retval;
}

CxmpInline void
tr_p_edge_dealloc(CxtTr *a_tr, CxtTrEdge a_edge)
{
    cw_tre_t *tre;
    cw_trr_t *trr;

    tre = &a_tr->tres[a_edge];
    if (tre->ps != NULL)
    {
	tr_p_ps_delete(a_tr, tre->ps);
    }
#ifdef CxmDebug
    memset(tre, 0x5a, sizeof(cw_tre_t));
#endif

    trr = &a_tr->trrs[a_edge << 1];
    if (trr->ps != NULL)
    {
	tr_p_ps_delete(a_tr, trr->ps);
    }
#ifdef CxmDebug
    memset(trr, 0x5a, sizeof(cw_trr_t));
#endif

    trr = &a_tr->trrs[(a_edge << 1) + 1];
    if (trr->ps != NULL)
    {
	tr_p_ps_delete(a_tr, trr->ps);
    }
#ifdef CxmDebug
    memset(trr, 0x5a, sizeof(cw_trr_t));
#endif

    a_tr->tres[a_edge].u.link = a_tr->sparetres;
    a_tr->sparetres = a_edge;
}

CxmpInline cw_tr_ring_t
tr_p_edge_ring_get(CxtTr *a_tr, CxtTrEdge a_edge, uint32_t a_end)
{
    return ((a_edge << 1) + a_end);
}

CxtTrEdge
CxTrEdgeNew(CxtTr *a_tr)
{
    return tr_p_edge_alloc(a_tr);
}

void
CxTrEdgeDelete(CxtTr *a_tr, CxtTrEdge a_edge)
{
    cxmDassert(tr_p_edge_validate(a_tr, a_edge));
    cxmAssert(tr_p_ring_node_get(a_tr, tr_p_edge_ring_get(a_tr, a_edge, 0))
	      == CxmTrNodeNone);
    cxmAssert(tr_p_ring_node_get(a_tr, tr_p_edge_ring_get(a_tr, a_edge, 1))
	      == CxmTrNodeNone);

    tr_p_edge_dealloc(a_tr, a_edge);
}

CxtTrNode
CxTrEdgeNodeGet(CxtTr *a_tr, CxtTrEdge a_edge, uint32_t a_end)
{
    cxmDassert(tr_p_edge_validate(a_tr, a_edge));
    cxmAssert(a_end == 0 || a_end == 1);

    return tr_p_ring_node_get(a_tr, tr_p_edge_ring_get(a_tr, a_edge, a_end));
}

void
CxTrEdgeNextGet(CxtTr *a_tr, CxtTrEdge a_edge, uint32_t a_end,
		 CxtTrEdge *r_next, uint32_t *r_end)
{
    cw_tr_ring_t ringind;

    cxmDassert(tr_p_edge_validate(a_tr, a_edge));
    cxmAssert(a_end == 0 || a_end == 1);

    ringind = CxmQriNext(a_tr->trrs, tr_p_edge_ring_get(a_tr, a_edge, a_end),
		       link);
    *r_next = tr_p_ring_edge_get(a_tr, ringind);
    *r_end = (ringind & 1);
}

void
CxTrEdgePrevGet(CxtTr *a_tr, CxtTrEdge a_edge, uint32_t a_end,
		 CxtTrEdge *r_prev, uint32_t *r_end)
{
    cw_tr_ring_t ringind;

    cxmDassert(tr_p_edge_validate(a_tr, a_edge));
    cxmAssert(a_end == 0 || a_end == 1);

    ringind = CxmQriPrev(a_tr->trrs, tr_p_edge_ring_get(a_tr, a_edge, a_end),
		       link);
    *r_prev = tr_p_ring_edge_get(a_tr, ringind);
    *r_end = (ringind & 1);
}

double
CxTrEdgeLengthGet(CxtTr *a_tr, CxtTrEdge a_edge)
{
    cxmDassert(tr_p_edge_validate(a_tr, a_edge));

    return a_tr->tres[a_edge].length;
}

void
CxTrEdgeLengthSet(CxtTr *a_tr, CxtTrEdge a_edge, double a_length)
{
    cxmDassert(tr_p_edge_validate(a_tr, a_edge));

    a_tr->tres[a_edge].length = a_length;
}

void *
CxTrEdgeAuxGet(CxtTr *a_tr, CxtTrEdge a_edge)
{
    cxmDassert(tr_p_edge_validate(a_tr, a_edge));

    return a_tr->tres[a_edge].u.aux;
}

void
CxTrEdgeAuxSet(CxtTr *a_tr, CxtTrEdge a_edge, void *a_aux)
{
    cxmDassert(tr_p_edge_validate(a_tr, a_edge));

    a_tr->tres[a_edge].u.aux = a_aux;
}

void
CxTrEdgeAttach(CxtTr *a_tr, CxtTrEdge a_edge, CxtTrNode a_node_a,
	       CxtTrNode a_node_b)
{
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    cxmDassert(tr_p_edge_validate(a_tr, a_edge));
    cxmDassert(CxTrEdgeNodeGet(a_tr, a_edge, 0) == CxmTrNodeNone);
    cxmDassert(CxTrEdgeNodeGet(a_tr, a_edge, 1) == CxmTrNodeNone);
    cxmDassert(tr_p_node_validate(a_tr, a_node_a));
    cxmDassert(tr_p_node_validate(a_tr, a_node_b));
    cxmAssert(a_node_a != a_node_b);
    cxmAssert(CxTrNodeDistance(a_tr, a_node_a, a_node_b) == 0);

    /* First end. */
    ring = tr_p_edge_ring_get(a_tr, a_edge, 0);
    trn = &a_tr->trns[a_node_a];
    CxmQliTailInsert(&trn->rings, a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = a_node_a;

    /* Second end. */
    ring = tr_p_edge_ring_get(a_tr, a_edge, 1);
    trn = &a_tr->trns[a_node_b];
    CxmQliTailInsert(&trn->rings, a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = a_node_b;

    /* Mark tree as modified. */
    a_tr->modified = true;

    cxmDassert(tr_p_edge_validate(a_tr, a_edge));
    cxmDassert(tr_p_node_validate(a_tr, a_node_a));
    cxmDassert(tr_p_node_validate(a_tr, a_node_b));
    cxmAssert(CxTrNodeDistance(a_tr, a_node_a, a_node_b) == 1);
}

void
CxTrEdgeDetach(CxtTr *a_tr, CxtTrEdge a_edge)
{
    cw_tr_ring_t ring;
    cw_trn_t *trn;

    cxmDassert(tr_p_edge_validate(a_tr, a_edge));
    cxmDassert(CxTrEdgeNodeGet(a_tr, a_edge, 0) != CxmTrNodeNone);
    cxmDassert(CxTrEdgeNodeGet(a_tr, a_edge, 1) != CxmTrNodeNone);
    cxmAssert(CxTrNodeDistance(a_tr, CxTrEdgeNodeGet(a_tr, a_edge, 0),
			       CxTrEdgeNodeGet(a_tr, a_edge, 1)) == 1);

    /* Detach from neighboring nodes.  Use CxmQliRemove() to make sure that the
     * nodes still point to their rings. */
    ring = tr_p_edge_ring_get(a_tr, a_edge, 0);
    trn = &a_tr->trns[tr_p_ring_node_get(a_tr, ring)];
    CxmQliRemove(&trn->rings, a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = CxmTrNodeNone;

    ring = tr_p_edge_ring_get(a_tr, a_edge, 1);
    trn = &a_tr->trns[tr_p_ring_node_get(a_tr, ring)];
    CxmQliRemove(&trn->rings, a_tr->trrs, ring, link);
    a_tr->trrs[ring].node = CxmTrNodeNone;

    /* Mark tree as modified. */
    a_tr->modified = true;

    cxmDassert(tr_p_edge_validate(a_tr, a_edge));
}

/******************************************************************************/

/* tr_node. */

CxmpInline void
tr_p_CxNodeInit(CxtTr *a_tr, CxtTrNode a_node)
{
    cw_trn_t *trn;

    trn = &a_tr->trns[a_node];

    trn->u.aux = NULL;
    trn->taxon_num = CxmTrNodeTaxonNone;
    CxmQliNew(&trn->rings);

#ifdef CxmDebug
    trn->magic = CW_TRN_MAGIC;
#endif
}

CxmpInline CxtTrNode
tr_p_node_alloc(CxtTr *a_tr)
{
    CxtTrNode retval;

    if (a_tr->sparetrns == CxmTrNodeNone)
    {
	uint32_t i, nspares;

	/* Allocate spares. */
	if (a_tr->trns == NULL)
	{
	    a_tr->trns
		= (cw_trn_t *) CxmMalloc(sizeof(cw_trn_t));
	    nspares = 1;
	    a_tr->ntrns = 1;
	}
	else
	{
	    a_tr->trns = (cw_trn_t *) CxmRealloc(a_tr->trns,
						 sizeof(cw_trn_t)
						 * a_tr->ntrns * 2);
	    nspares = a_tr->ntrns;
	    a_tr->ntrns *= 2;
	}

	/* Initialize last spare. */
	a_tr->sparetrns = a_tr->ntrns - 1;
	a_tr->trns[a_tr->sparetrns].u.link = CxmTrNodeNone;

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
    tr_p_CxNodeInit(a_tr, retval);

    return retval;
}

CxmpInline void
tr_p_node_dealloc(CxtTr *a_tr, CxtTrNode a_node)
{
    cw_trn_t *trn;

    trn = &a_tr->trns[a_node];
    
#ifdef CxmDebug
    memset(trn, 0x5a, sizeof(cw_trn_t));
#endif

    a_tr->trns[a_node].u.link = a_tr->sparetrns;
    a_tr->sparetrns = a_node;
}

/* Calculate the number of edges connected to the node that a_ring is connected
 * to. */
CxmpInline uint32_t
tr_p_CxNodeDegree(CxtTr *a_tr, cw_tr_ring_t a_ring)
{
    uint32_t retval;
    cw_tr_ring_t ring;

    retval = 1;
    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
    {
	retval++;
    }

    return retval;
}

/* Calculate the number of edges between two nodes.  A distance of 0 means that
 * there is no path between the two nodes. */
static uint32_t
tr_p_node_distance(CxtTr *a_tr, cw_tr_ring_t a_ring, CxtTrNode a_other,
		   uint32_t a_distance)
{
    uint32_t retval;
    cw_tr_ring_t ring;

    if (tr_p_ring_node_get(a_tr, a_ring) == a_other)
    {
	retval = a_distance;
	goto RETURN;
    }

    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
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

CxtTrNode
CxTrNodeNew(CxtTr *a_tr)
{
    return tr_p_node_alloc(a_tr);
}

void
CxTrNodeDelete(CxtTr *a_tr, CxtTrNode a_node)
{
    cxmDassert(tr_p_node_validate(a_tr, a_node));
    cxmAssert(CxmQliFirst(&a_tr->trns[a_node].rings) == CW_TR_RING_NONE);

    tr_p_node_dealloc(a_tr, a_node);
}

uint32_t
CxTrNodeTaxonNumGet(CxtTr *a_tr, CxtTrNode a_node)
{
    cxmDassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].taxon_num;
}

void
CxTrNodeTaxonNumSet(CxtTr *a_tr, CxtTrNode a_node,
		      uint32_t a_taxon_num)
{
    cxmDassert(tr_p_node_validate(a_tr, a_node));

    a_tr->trns[a_node].taxon_num = a_taxon_num;

    a_tr->modified = true;
}

void
CxTrNodeEdgeGet(CxtTr *a_tr, CxtTrNode a_node,
		 CxtTrEdge *r_edge, uint32_t *r_end)
{
    cw_trn_t *trn;
    cw_tr_ring_t ringind;

    cxmDassert(tr_p_node_validate(a_tr, a_node));

    trn = &a_tr->trns[a_node];

    ringind = CxmQliFirst(&trn->rings);
    if (ringind != CxmTrEdgeNone)
    {
	*r_edge = tr_p_ring_edge_get(a_tr, ringind);
	if (r_end != NULL)
	{
	    *r_end = (ringind & 1);
	}
    }
    else
    {
	*r_edge = CxmTrEdgeNone;
    }
}

void *
CxTrNodeAuxGet(CxtTr *a_tr, CxtTrNode a_node)
{
    cxmDassert(tr_p_node_validate(a_tr, a_node));

    return a_tr->trns[a_node].u.aux;
}

void
CxTrNodeAuxSet(CxtTr *a_tr, CxtTrNode a_node, void *a_aux)
{
    cxmDassert(tr_p_node_validate(a_tr, a_node));

    a_tr->trns[a_node].u.aux = a_aux;
}

uint32_t
CxTrNodeDegree(CxtTr *a_tr, CxtTrNode a_node)
{
    uint32_t retval;
    cw_tr_ring_t ring;

    cxmDassert(tr_p_node_validate(a_tr, a_node));

    ring = CxmQliFirst(&a_tr->trns[a_node].rings);
    if (ring != CW_TR_RING_NONE)
    {
	retval = tr_p_CxNodeDegree(a_tr, ring);
    }
    else
    {
	retval = 0;
    }

    return retval;
}

uint32_t
CxTrNodeDistance(CxtTr *a_tr, CxtTrNode a_node, CxtTrNode a_other)
{
    uint32_t retval;
    cw_trn_t *trn;
    cw_tr_ring_t ring;

    cxmDassert(tr_p_node_validate(a_tr, a_node));
    cxmDassert(tr_p_node_validate(a_tr, a_other));
    cxmAssert(a_node != a_other);

    trn = &a_tr->trns[a_node];

    ring = CxmQliFirst(&trn->rings);
    if (ring != CW_TR_RING_NONE)
    {
	CxmQliForeach(ring, &trn->rings, a_tr->trrs, link)
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
CxmpInline void
tr_p_new(CxtTr *a_tr)
{
    a_tr->aux = NULL;
    a_tr->modified = false;
    a_tr->base = CxmTrNodeNone;
    a_tr->ntaxa = 0;
    a_tr->nedges = 0;
    a_tr->bedges = NULL;
    a_tr->nbedges_a = 0;
    a_tr->nbedges_b = 0;
    a_tr->trt = NULL;
    a_tr->trtused = 0;
    a_tr->trns = NULL;
    a_tr->ntrns = 0;
    a_tr->sparetrns = CxmTrNodeNone;
    a_tr->tres = NULL;
    a_tr->ntres = 0;
    a_tr->sparetres = CxmTrEdgeNone;
    a_tr->trrs = NULL;
    a_tr->held = NULL;
    a_tr->heldlen = 0;
    a_tr->nheld = 0;

#ifdef CxmDebug
    a_tr->magic = CW_TR_MAGIC;
#endif
}

/* Recursively traverse the tree, count the number of taxa, and find the lowest
 * numbered taxon. */
static CxtTrNode
tr_p_lowest_recurse(CxtTr *a_tr, cw_tr_ring_t a_ring, uint32_t *r_ntaxa,
		    uint32_t *r_nedges, CxtTrNode a_root)
{
    CxtTrNode retval, node, root, troot;
    cw_tr_ring_t ring;
    cw_trn_t *trn;

    node = tr_p_ring_node_get(a_tr, a_ring);
    trn = &a_tr->trns[node];

    if (trn->taxon_num != CxmTrNodeTaxonNone)
    {
	/* Leaf node. */
	(*r_ntaxa)++;
    }

    if (trn->taxon_num != CxmTrNodeTaxonNone
	&& (a_root == CxmTrNodeNone
	    || trn->taxon_num < a_tr->trns[a_root].taxon_num))
    {
	retval = node;
	root = node;
    }
    else
    {
	retval = CxmTrNodeNone;
	root = a_root;
    }

    /* Iterate over neighbors. */
    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
    {
	/* Count edge. */
	(*r_nedges)++;

	troot = tr_p_lowest_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				    r_ntaxa, r_nedges, root);
	if (troot != CxmTrNodeNone)
	{
	    retval = troot;
	    root = troot;
	}
    }

    return retval;
}

/* Recursively traverse the tree, count the number of taxa, and find the lowest
 * numbered taxon. */
static CxtTrNode
tr_p_lowest(CxtTr *a_tr, CxtTrNode a_node, uint32_t *r_ntaxa,
	    uint32_t *r_nedges)
{
    CxtTrNode retval, root, troot;
    cw_tr_ring_t ring;
    cw_trn_t *trn;

    cxmAssert(a_node != CxmTrNodeNone);

    trn = &a_tr->trns[a_node];

    if (trn->taxon_num != CxmTrNodeTaxonNone)
    {
	/* Leaf node. */
	(*r_ntaxa)++;
    }

    if (trn->taxon_num != CxmTrNodeTaxonNone)
    {
	retval = a_node;
	root = a_node;
    }
    else
    {
	retval = CxmTrNodeNone;
	root = CxmTrNodeNone;
    }

    /* Iterate over neighbors. */
    CxmQliForeach(ring, &trn->rings, a_tr->trrs, link)
    {
	/* Count edge. */
	(*r_nedges)++;
	
	troot = tr_p_lowest_recurse(a_tr, tr_p_ring_other_get(a_tr, ring),
				    r_ntaxa, r_nedges, root);
	if (troot != CxmTrNodeNone)
	{
	    retval = troot;
	    root = troot;
	}
    }

    return retval;
}

#ifdef CxmDebug
/* Validate a tree. */
static bool
tr_p_validate(CxtTr *a_tr)
{
    uint32_t i, ntaxa, nedges;

    cxmCheckPtr(a_tr);
    cxmAssert(a_tr->magic == CW_TR_MAGIC);
    cxmAssert(a_tr->modified == false);

    ntaxa = 0;
    nedges = 0;
    if (a_tr->base != CxmTrNodeNone)
    {
	tr_p_lowest(a_tr, a_tr->base, &ntaxa, &nedges);
    }
    cxmAssert(a_tr->ntaxa == ntaxa);
    cxmAssert(a_tr->nedges == nedges);

    /* Iterate over trns and do some basic sanity checks. */
    for (i = 0; i < a_tr->ntrns; i++)
    {
	if (a_tr->trns[i].magic == CW_TRN_MAGIC)
	{
	    tr_p_node_validate(a_tr, (CxtTrNode) i);
	}
	else
	{
	    /* Make sure there are no valid trn's in the free list. */
	    cxmAssert(a_tr->trns[i].u.link == CxmTrNodeNone
		      || a_tr->trns[a_tr->trns[i].u.link].magic
		      != CW_TRN_MAGIC);
	}
    }

    cxmAssert(a_tr->sparetrns == CxmTrNodeNone
	      || a_tr->trns[a_tr->sparetrns].magic != CW_TRN_MAGIC);

    return true;
}
#endif

static void
tr_p_ntaxa_nedges_update(CxtTr *a_tr)
{
    uint32_t ntaxa, nedges;

    /* Update ntaxa and nedges. */
    ntaxa = 0;
    nedges = 0;
    if (a_tr->base != CxmTrNodeNone)
    {
	tr_p_lowest(a_tr, a_tr->base, &ntaxa, &nedges);
    }

    a_tr->ntaxa = ntaxa;
    a_tr->nedges = nedges;
}

static void
tr_p_bisection_edge_list_gen_recurse(CxtTr *a_tr, cw_tr_ring_t a_ring,
				     CxtTrEdge *ar_edges,
				     uint32_t *ar_nedges)
{
    cw_tr_ring_t ring;

    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
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
CxmpInline CxtTrNode
tr_p_bisection_edge_list_gen(CxtTr *a_tr, cw_tr_ring_t a_ring,
			     CxtTrEdge *ar_edges, uint32_t *ar_nedges)
{
    CxtTrNode retval;
    cw_tr_ring_t ring;

    /* Initialize the length of the list before recursing. */
    *ar_nedges = 0;

    switch (tr_p_CxNodeDegree(a_tr, a_ring))
    {
	case 1:
	{
	    /* A subtree that is composed of a single node has no edges.  Add a
	     * single entry to the list, and return the node. */
	    ar_edges[0] = CxmTrEdgeNone;
	    (*ar_nedges)++;
	    retval = tr_p_ring_node_get(a_tr, a_ring);
	    break;
	}
	case 2:
	{
	    /* A tree should never have nodes of degree 2. */
	    cxmNotReached();
	}
	case 3:
	{
	    /* Take care to add only one of the edges that is connected to the
	     * node, since from the perspective of TBR, the node does not exist.
	     * (A node of degree 2 is a superfluous node.) */

	    /* First edge. */
	    ring = CxmQriNext(a_tr->trrs, a_ring, link);
	    ar_edges[0] = tr_p_ring_edge_get(a_tr, ring);
	    (*ar_nedges)++;

	    /* First subtree. */
	    tr_p_bisection_edge_list_gen_recurse(a_tr,
						 tr_p_ring_other_get(a_tr,
								     ring),
						 ar_edges, ar_nedges);

	    /* Second subtree. */
	    ring = CxmQriNext(a_tr->trrs, ring, link);
	    tr_p_bisection_edge_list_gen_recurse(a_tr,
						 tr_p_ring_other_get(a_tr,
								     ring),
						 ar_edges, ar_nedges);

	    retval = CxmTrNodeNone;
	    break;
	}
	default:
	{
	    /* Add all edges in the subtree.  Removing the bisection edge still
	     * leaves enough edges attached to the node for the node to have
	     * relevance. */
	    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
	    {
		/* Add edge to list. */
		ar_edges[*ar_nedges] = tr_p_ring_edge_get(a_tr, ring);
		(*ar_nedges)++;

		tr_p_bisection_edge_list_gen_recurse(a_tr,
						     tr_p_ring_other_get(a_tr,
									 ring),
						     ar_edges, ar_nedges);
	    }

	    retval = CxmTrNodeNone;
	    break;
	}
    }

    return retval;
}

/* Generate lists of edges in each half of a logical bisection at edge
 * a_bisect. */
CxmpInline void
tr_p_bedges_gen(CxtTr *a_tr, CxtTrEdge a_bisect, CxtTrNode *r_node_a,
		CxtTrNode *r_node_b)
{
    CxtTrNode node_a, node_b;

    cxmDassert(tr_p_edge_validate(a_tr, a_bisect));

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
tr_p_trt_bisect_edge_update_recurse(CxtTr *a_tr, cw_tr_ring_t a_ring,
				    uint32_t *ar_edge_count)
{
    cw_tr_ring_t ring;

    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
    {
	/* Record edge. */
	a_tr->trt[*ar_edge_count].bisect_edge = tr_p_ring_edge_get(a_tr, ring);
	(*ar_edge_count)++;

	/* Recurse into neighbor subtree. */
	tr_p_trt_bisect_edge_update_recurse(a_tr,
					    tr_p_ring_other_get(a_tr, ring),
					    ar_edge_count);
    }
}

static void
tr_p_trt_update(CxtTr *a_tr, uint32_t a_nedges_prev)
{
    uint32_t i, j, n, offset;

    cxmAssert(a_tr->modified == false);

    /* Allocate/reallocate/deallocate trt. */
    if (a_tr->trt == NULL)
    {
	/* Allocate trt. */
	a_tr->trt = (cw_trt_t *) CxmMalloc(sizeof(cw_trt_t)
					   * (a_tr->nedges + 1));
    }
    else if (a_tr->nedges != a_nedges_prev)
    {
	/* Reallocate trt.  There is never a need to deallocate trt here,
	 * since trt contains one extra element. */
	a_tr->trt = (cw_trt_t *) CxmRealloc(a_tr->trt,
					    sizeof(cw_trt_t)
					    * (a_tr->nedges + 1));
    }

    /* Recursively traverse the tree, and initialize trt->bisect_edge along the
     * way. */
    if (a_tr->nedges > 0)
    {
	uint32_t edge_count;
	cw_trn_t *trn;
	cw_tr_ring_t ring;

	cxmAssert(a_tr->base != CxmTrNodeNone);

	edge_count = 0;
	trn = &a_tr->trns[a_tr->base];
	CxmQliForeach(ring, &trn->rings, a_tr->trrs, link)
	{
	    /* Record edge. */
	    a_tr->trt[edge_count].bisect_edge = tr_p_ring_edge_get(a_tr, ring);
	    edge_count++;

	    tr_p_trt_bisect_edge_update_recurse(a_tr,
						tr_p_ring_other_get(a_tr, ring),
						&edge_count);
	}
	cxmAssert(edge_count == a_tr->nedges);
    }

    /* Iteratively fill in trt. */
    for (i = j = offset = 0; i < a_tr->nedges; i++)
    {
	/* Record offset. */
	a_tr->trt[j].offset = offset;

	/* Update offset. */
	tr_p_bedges_gen(a_tr, a_tr->trt[i].bisect_edge, NULL, NULL);
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
tr_p_bedges_update(CxtTr *a_tr, uint32_t a_nedges_prev)
{
    /* Allocate/reallocate/deallocate bedges.  To keep things simple, allocate
     * as big an array as there are edges, even though not quite that many are
     * ever used. */
    if (a_tr->bedges == NULL)
    {
	/* Allocate bedges. */
	a_tr->bedges
	    = (CxtTrEdge *) CxmMalloc(sizeof(CxtTrEdge) * a_tr->nedges);
    }
    else if (a_tr->nedges != a_nedges_prev)
    {
	if (a_tr->nedges > 0)
	{
	    /* Reallocate bedges. */
	    a_tr->bedges = (CxtTrEdge *) CxmRealloc(a_tr->bedges,
						       sizeof(CxtTrEdge)
						       * a_tr->nedges);
	}
	else
	{
	    /* Deallocate bedges. */
	    CxmFree(a_tr->bedges);
	    a_tr->bedges = NULL;
	}
    }

    /* Clear nbedges_[ab]. */
    a_tr->nbedges_a = 0;
    a_tr->nbedges_b = 0;
}

CxmpInline void
tr_p_update(CxtTr *a_tr)
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
struct cw_CxTrCanonize_s
{
    cw_tr_ring_t ring;
    uint32_t min_taxon;
};

/* Comparison function that is passed to qsort(). */
static int
tr_p_canonize_compare(const void *a_a, const void *a_b)
{
    const struct cw_CxTrCanonize_s *a = (const struct cw_CxTrCanonize_s *) a_a;
    const struct cw_CxTrCanonize_s *b = (const struct cw_CxTrCanonize_s *) a_b;

    if (a->min_taxon < b->min_taxon)
    {
	return -1;
    }
    else
    {
	cxmAssert(a->min_taxon > b->min_taxon);
	return 1;
    }
}

/* Convert a tree node to canonical form by re-ordering the ring such that
 * subtrees are in increasing order of minimum taxon number contained. */
static uint32_t
tr_p_canonize(CxtTr *a_tr, cw_tr_ring_t a_ring)
{
    uint32_t retval, degree;
    CxtTrNode node;

    /* Get taxon number (an internal node has CxmTrNodeTaxonNone). */
    cxmDassert(tr_p_node_validate(a_tr, tr_p_ring_node_get(a_tr, a_ring)));
    node = tr_p_ring_node_get(a_tr, a_ring);
    retval = CxTrNodeTaxonNumGet(a_tr, node);

    /* Get the degree of the node that this ring is a part of. */
    degree = tr_p_CxNodeDegree(a_tr, a_ring);

    if (degree > 1)
    {
	uint32_t i, min_taxon;
	cw_tr_ring_t ring;
	struct cw_CxTrCanonize_s *canonize;

	/* Allocate space for a temporary array that can be used to sort the
	 * ring. */
	canonize = (struct cw_CxTrCanonize_s *)
	    CxmMalloc(sizeof(struct cw_CxTrCanonize_s) * (degree - 1));

	/* Iteratively canonize subtrees, keeping track of the minimum taxon
	 * number seen overall, as well as for each subtree. */
	i = 0;
	retval = CxmTrNodeTaxonNone;
	CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
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
	cxmAssert(i == degree - 1);

	/* Sort the subtrees. */
	qsort(canonize, degree - 1, sizeof(struct cw_CxTrCanonize_s),
	      tr_p_canonize_compare);

	/* Set the beginning of the ring to a_ring.  This makes it easier for
	 * external code to traverse a tree in canonical order. */
	CxmQliFirst(&a_tr->trns[node].rings) = a_ring;

	/* Re-arrange the ring.  The first element can be skipped, since the
	 * removal/re-insertion of all other elements eventually leaves the
	 * first element in the proper location. */
	for (i = 1; i < (degree - 1); i++)
	{
	    CxmQriRemove(a_tr->trrs, canonize[i].ring, link);
	    CxmQriBeforeInsert(a_tr->trrs, a_ring, canonize[i].ring, link);
	}

	/* Clean up. */
	CxmFree(canonize);
    }

    return retval;
}

/* As part of TBR, extract a node that has only two neighbors.  Take care to
 * leave reconnection edges in the tree.  Return CxmTrNodeNone, unless there
 * is only one node in the subtree; in that case, return the node so that it can
 * be used directly during reconnection. */
CxmpInline CxtTrNode
tr_p_tbr_node_extract(CxtTr *a_tr, CxtTrNode a_node,
		      CxtTrEdge a_reconnect_a, CxtTrEdge a_reconnect_b,
		      CxtTrEdge *ar_tedges, uint32_t *ar_ntedges,
		      CxtTrNode *ar_tnodes, uint32_t *ar_ntnodes)
{
    CxtTrNode retval;

    switch (CxTrNodeDegree(a_tr, a_node))
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
	    cxmNotReached();
	}
	case 2:
	{
	    cw_tr_ring_t ring_a, ring_b;
	    cw_tr_ring_t ring_a_other, ring_b_other;
	    CxtTrEdge edge_a, edge_b;
	    CxtTrNode node_a, node_b;
	    cw_tr_ps_t *tps;

	    /* Get all variables that are necessary for careful extraction of
	     * a_node, and proper rematching of rings with nodes.  The
	     * rematching is critical to the maintenance of the character state
	     * sets in leaf nodes (which node_[ab] may or may not be). */
	    ring_a = CxmQliFirst(&a_tr->trns[a_node].rings);
	    edge_a = tr_p_ring_edge_get(a_tr, ring_a);
	    ring_a_other = tr_p_ring_other_get(a_tr, ring_a);
	    node_a = tr_p_ring_node_get(a_tr, ring_a_other);

	    ring_b = CxmQriNext(a_tr->trrs, ring_a, link);
	    edge_b = tr_p_ring_edge_get(a_tr, ring_b);
	    ring_b_other = tr_p_ring_other_get(a_tr, ring_b);
	    node_b = tr_p_ring_node_get(a_tr, ring_b_other);

	    /* Detach. */
	    CxTrEdgeDetach(a_tr, edge_a);
	    CxTrEdgeDetach(a_tr, edge_b);

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
		    CxTrEdgeAttach(a_tr, edge_a, node_a, node_b);
		}
		else
		{
		    CxTrEdgeAttach(a_tr, edge_a, node_b, node_a);
		}

		/* Store edge_b as a spare. */
		ar_tedges[*ar_ntedges] = edge_b;
		(*ar_ntedges)++;
	    }
	    else
	    {
		/* Use edge_b when splicing node_[ab] back together. */
		cxmAssert(edge_a != a_reconnect_a && edge_a != a_reconnect_b);

		/* Swap data in ring_b and ring_a_other. */
		tps = a_tr->trrs[ring_b].ps;
		a_tr->trrs[ring_b].ps = a_tr->trrs[ring_a_other].ps;
		a_tr->trrs[ring_a_other].ps = tps;

		/* Attach node_[ab].  Take care to keep the proper ends of
		 * edge_b associated with node_[ab]. */
		if (ring_b < ring_b_other)
		{
		    CxTrEdgeAttach(a_tr, edge_b, node_a, node_b);
		}
		else
		{
		    CxTrEdgeAttach(a_tr, edge_b, node_b, node_a);
		}

		/* Store edge_a as a spare. */
		ar_tedges[*ar_ntedges] = edge_a;
		(*ar_ntedges)++;
	    }

	    retval = CxmTrNodeNone;
	    break;
	}
	default:
	{
	    /* Do nothing, since this node has enough neighbors to remain
	     * relevant (3 or more). */
	    retval = CxmTrNodeNone;
	}
    }

    return retval;
}

/* Splice a node into the middle of a_edge, and return the node. */
CxmpInline CxtTrNode
tr_p_tbr_node_splice(CxtTr *a_tr, CxtTrEdge a_edge, 
		     CxtTrEdge *ar_tedges, uint32_t *ar_ntedges,
		     CxtTrNode *ar_tnodes, uint32_t *ar_ntnodes)
{
    CxtTrNode retval, node_a, node_b;
    cw_tr_ring_t ring_a, ring_b, ring;
    CxtTrEdge edge;
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
	// XXX edge = tr_p_edge_wrapped_new(a_tr);
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
	// XXX retval = tr_p_node_wrapped_new(a_tr);
    }

    /* Detach. */
    CxTrEdgeDetach(a_tr, a_edge);

    /* Swap data in ring_b and ring. */
    tps = a_tr->trrs[ring_b].ps;
    a_tr->trrs[ring_b].ps = a_tr->trrs[ring].ps;
    a_tr->trrs[ring].ps = tps;

    /* Reattach. */
    CxTrEdgeAttach(a_tr, a_edge, node_a, retval);
    CxTrEdgeAttach(a_tr, edge, node_b, retval);

    return retval;
}

static void
tr_p_mp_ring_prepare(CxtTr *a_tr, cw_tr_ring_t a_ring, char *a_taxa[],
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
    if (taxon_num != CxmTrNodeTaxonNone)
    {
	uint32_t i, j;
	char *chars;

	trr->ps->subtrees_score = 0;
	trr->ps->node_score = 0;

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
		    cxmNotReached();
		}
	    }
	    j++;
	}
    }
}

static void
tr_p_mp_prepare_recurse(CxtTr *a_tr, cw_tr_ring_t a_ring,
			char *a_taxa[], uint32_t a_ntaxa, uint32_t a_nchars,
			bool *a_chars_mask, uint32_t a_ninformative)
{
    cw_tr_ring_t ring;
    cw_tre_t *tre;

    /* Prepare a_ring. */
    tr_p_mp_ring_prepare(a_tr, a_ring, a_taxa, a_ntaxa, a_nchars,
			 a_chars_mask, a_ninformative);

    /* Recurse into subtrees. */
    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
    {
	/* Prepare edge before recursing. */
	tre = &a_tr->tres[tr_p_ring_edge_get(a_tr, ring)];
	if (tre->ps == NULL)
	{
	    tre->ps = tr_p_ps_new(a_tr);
	}
	tr_p_ps_prepare(a_tr, tre->ps, a_ninformative);

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
tr_p_mp_ring_finish(CxtTr *a_tr, cw_tr_ring_t a_ring)
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
tr_p_mp_finish_recurse(CxtTr *a_tr, cw_tr_ring_t a_ring)
{
    cw_tr_ring_t ring;
    cw_tre_t *tre;

    /* Clean up a_ring. */
    tr_p_mp_ring_finish(a_tr, a_ring);

    /* Recurse into subtrees. */
    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
    {
	/* Clean up edge before recursing. */
	tre = &a_tr->tres[tr_p_ring_edge_get(a_tr, ring)];
	if (tre->ps != NULL)
	{
	    tr_p_ps_delete(a_tr, tre->ps);
	    tre->ps = NULL;
	}

	/* Clean up ring. */
	tr_p_mp_ring_finish(a_tr, ring);

	/* Recurse. */
	tr_p_mp_finish_recurse(a_tr, tr_p_ring_other_get(a_tr, ring));
    }
}

#ifdef CW_CPU_IA32
CxmpInline void
tr_p_mp_ia32_pscore(CxtTr *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
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
	/* Use SSE2 to evaluate the characters.  This loop handles 32 characters
	 * per iteration. */
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

    a_p->node_score = ns;
}
#endif

static void
tr_p_mp_c_pscore(CxtTr *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
		 cw_tr_ps_t *a_b)
{
    uint32_t i, nwords, ns, a, b, m, r, un;
    uint32_t *chars_p, *chars_a, *chars_b;
    static const uint32_t bits_table[] =
	{
	    2, 1, -1, -1, -1, -1, -1, -1,
	    -1, -1, -1, -1, -1, -1, -1, -1,
	    1, 0
	};

    /* Calculate node score. */
    ns = 0;

#define MP_C_PSCORE_INNER()						\
    a = chars_a[i];							\
    b = chars_b[i];							\
									\
    /* Get 1's in the least significant bits of state sets that are	\
     * non-empty after the intersection operation. */			\
    r = m = a & b;							\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
									\
    /* Mask out garbage. */						\
    m &= 0x11111111;							\
									\
    /* Count up changes. */						\
    ns += bits_table[m & 0xff]						\
	+ bits_table[(m >> 8) & 0xff]					\
	+ bits_table[(m >> 16) & 0xff]					\
	+ bits_table[(m >> 24) & 0xff];					\
									\
    /* Propagate 1's to make a bit mask. */				\
    m |= (m << 1);							\
    m |= (m << 1);							\
    m |= (m << 1);							\
									\
    /* Combine results of intersection and union operations. */		\
    r &= m;								\
    un = a | b;								\
    m = (~m);								\
    un &= m;								\
    r |= un;								\
									\
    /* Store result. */							\
    chars_p[i] = r;							\
									\
    i++;

    /* Calculate preliminary Fitch parsimony scores for each character. */
    chars_p = (uint32_t *) a_p->chars;
    chars_a = (uint32_t *) a_a->chars;
    chars_b = (uint32_t *) a_b->chars;
    for (i = 0, nwords = (a_p->nchars >> 3); i < nwords;)
    {
	MP_C_PSCORE_INNER();
	MP_C_PSCORE_INNER();
	MP_C_PSCORE_INNER();
	MP_C_PSCORE_INNER();
    }
#undef MP_C_PSCORE_INNER

    a_p->node_score = ns;
}

/* Unconditionally calculate the partial score for a_p, using a_a and a_b as
 * children. */
CxmpInline void
tr_p_mp_pscore(CxtTr *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b)
{
    /* Reset this node's parent pointer, to keep the parent from using an
     * invalid cached value. */
    a_p->parent = NULL;

    /* Calculate sum of subtree scores. */
    a_p->subtrees_score
	= a_a->subtrees_score + a_a->node_score
	+ a_b->subtrees_score + a_b->node_score;

#ifdef CW_CPU_IA32
    if (CxgIa32UseSse2)
    {
	tr_p_mp_ia32_pscore(a_tr, a_p, a_a, a_b);
    }
    else
#endif
    {
	tr_p_mp_c_pscore(a_tr, a_p, a_a, a_b);
    }
}

/* The sole purpose of this function is to assure that the contents of
 * tr_p_mp_pscore() are not inlined in tr_p_mp_cache_pscore().  Most of the
 * time, the cache should be usable, so the actual scoring code doesn't usually
 * get called. */
static void
tr_p_no_inline_mp_pscore(CxtTr *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
			cw_tr_ps_t *a_b)
{
    tr_p_mp_pscore(a_tr, a_p, a_a, a_b);
}

/* Calculate the partial score for a_p, using a_a and a_b as children.  However,
 * do some extra bookkeeping in order to be able to cache the results, and later
 * recognize that precisely the same calculation was cached. */
CxmpInline void
tr_p_mp_cache_pscore(CxtTr *a_tr, cw_tr_ps_t *a_p, cw_tr_ps_t *a_a,
		     cw_tr_ps_t *a_b)
{
//#define CW_TR_MP_CACHE_PSCORE_VALIDATE
#ifdef CW_TR_MP_CACHE_PSCORE_VALIDATE
    bool cached;
    uint32_t cached_node_score;
#endif

    cxmCheckPtr(a_p);
    cxmCheckPtr(a_a);
    cxmCheckPtr(a_b);

    /* Only calculate the parent's node score if the cached value is invalid. */
    if (a_a->parent != a_p || a_b->parent != a_p)
#ifdef CW_TR_MP_CACHE_PSCORE_VALIDATE
    {
	cached = false;
    }
    else
    {
	cached = true;
	cached_node_score = a_p->node_score;

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
    {
	/* Set parent pointers, so that cached values may be used in future
	 * runs. */
	a_a->parent = a_p;
	a_b->parent = a_p;

	/* Calculate the partial score. */
	tr_p_no_inline_mp_pscore(a_tr, a_p, a_a, a_b);
    }

#ifdef CW_TR_MP_CACHE_PSCORE_VALIDATE
    if (cached)
    {
	if (cached_node_score != a_p->node_score)
	{
	    fprintf(stderr, "%s:%d:%s(): node_score %u (should be %u)\n",
		    __FILE__, __LINE__, __FUNCTION__,
		    cached_node_score, a_p->node_score);
	    abort();
	}
    }
#endif
}

CxmpInline void
tr_p_mp_cache_invalidate(CxtTr *a_tr, cw_tr_ps_t *a_ps)
{
    cxmCheckPtr(a_ps);

    /* Reset this node's parent pointer, to keep the old parent from using an
     * invalid cached value. */
    a_ps->parent = NULL;
}

#ifdef CW_CPU_IA32
CxmpInline uint32_t
tr_p_mp_ia32_fscore(CxtTr *a_tr, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b,
		    uint32_t a_maxscore)
{
    uint32_t retval, i, nbytes, pns;
    cw_trc_t *chars_a, *chars_b;

    /* Calculate sum of subtree scores. */
    retval
	= a_a->subtrees_score + a_a->node_score
	+ a_b->subtrees_score + a_b->node_score;

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

    /* Use SSE2 to evaluate the characters.  This loop handles 32 characters per
     * iteration. */
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
tr_p_mp_c_fscore(CxtTr *a_tr, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b,
		 uint32_t a_maxscore)
{
    uint32_t retval, i, nwords, a, b, m;
    uint32_t *chars_a, *chars_b;
    static const uint32_t bits_table[] =
	{
	    2, 1, -1, -1, -1, -1, -1, -1,
	    -1, -1, -1, -1, -1, -1, -1, -1,
	    1, 0
	};

    /* Calculate sum of subtree scores. */
    retval
	= a_a->subtrees_score + a_a->node_score
	+ a_b->subtrees_score + a_b->node_score;

#define MP_C_FSCORE_INNER()						\
    a = chars_a[i];							\
    b = chars_b[i];							\
									\
    /* Get 1's in the least significant bits of state sets that are	\
     * non-empty after the intersection operation. */			\
    m = a & b;								\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
    m |= (m >> 1);							\
									\
    /* Mask out garbage. */						\
    m &= 0x11111111;							\
									\
    /* Count up changes. */						\
    retval += bits_table[m & 0xff]					\
	+ bits_table[(m >> 8) & 0xff]					\
	+ bits_table[(m >> 16) & 0xff]					\
	+ bits_table[(m >> 24) & 0xff];					\
									\
    if (retval > a_maxscore)						\
    {									\
	retval = UINT_MAX;						\
	break;								\
    }									\
									\
    i++;

    /* Calculate partial Fitch parsimony scores for each character. */
    chars_a = (uint32_t *) a_a->chars;
    chars_b = (uint32_t *) a_b->chars;
    for (i = 0, nwords = (a_a->nchars >> 3); i < nwords;)
    {
	MP_C_FSCORE_INNER();
	MP_C_FSCORE_INNER();
	MP_C_FSCORE_INNER();
	MP_C_FSCORE_INNER();
    }
#undef MP_C_FSCORE_INNER

    return retval;
}

/* Unconditionally calculate the final score of a tree, using a_a and a_b as
 * children. */
CxmpInline uint32_t
tr_p_mp_fscore(CxtTr *a_tr, cw_tr_ps_t *a_a, cw_tr_ps_t *a_b,
	       uint32_t a_maxscore)
{
    uint32_t retval;

#ifdef CW_CPU_IA32
    if (CxgIa32UseSse2)
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
tr_p_mp_score_recurse(CxtTr *a_tr, cw_tr_ring_t a_ring, CxtTrEdge a_bisect)
{
    cw_tr_ps_t *retval
#ifdef CxmCcSilence
	= NULL
#endif
	;
    uint32_t degree;
    bool adjacent;
    cw_tr_ring_t ring;

    /* Get the degree of the node.  Don't count the bisection edge (only an
     * issue if this node is adjacent to the bisection). */
    degree = 1;
    adjacent = false;
    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
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
	    cxmAssert(adjacent);

	    /* Clear the cache for the view that is being bypassed.  This is
	     * critical to correctness of the caching machinery, since each view
	     * should never be claimed as the parent of more than two other
	     * views. */
	    tr_p_mp_cache_invalidate(a_tr, a_tr->trrs[a_ring].ps);

	    /* Get the ring element that connects to the other portion of the
	     * subtree on this side of the bisection. */
	    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
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
		ring = CxmQriNext(a_tr->trrs, a_ring, link);
		ps_a = tr_p_mp_score_recurse(a_tr,
					     tr_p_ring_other_get(a_tr, ring),
					     a_bisect);

		ring = CxmQriNext(a_tr->trrs, ring, link);
		ps_b = tr_p_mp_score_recurse(a_tr,
					     tr_p_ring_other_get(a_tr, ring),
					     a_bisect);

		/* Calculate the partial score for this node. */
		retval = a_tr->trrs[a_ring].ps;
		tr_p_mp_cache_pscore(a_tr, retval, ps_a, ps_b);

		break;
	    }
	    /* Fall through if this node is adjacent to the bisection. */
	}
	default:
	{
	    /* This is a multifurcating node. */
	    cxmError("XXX Not implemented");
	}
    }

    return retval;
}

static void
tr_p_mp_views_recurse(CxtTr *a_tr, cw_tr_ring_t a_ring, cw_tr_ps_t *a_ps,
		      CxtTrEdge a_bisect)
{
    uint32_t degree;
    bool adjacent;
    cw_tr_ring_t ring;

    /* Get the degree of the node.  Don't count the bisection edge (only an
     * issue if this node is adjacent to the bisection). */
    degree = 1;
    adjacent = false;
    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
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
	    cxmAssert(adjacent);

	    /* Get the ring element that connects to the other portion of the
	     * subtree on this side of the bisection. */
	    CxmQriOthersForeach(ring, a_tr->trrs, a_ring, link)
	    {
		if (tr_p_ring_edge_get(a_tr, ring) != a_bisect)
		{
		    /* Clear the cache for the view that is being bypassed.
		     * This is critical to correctness of the caching machinery,
		     * since each view should never be claimed as the parent of
		     * more than two other views. */
		    tr_p_mp_cache_invalidate(a_tr, a_tr->trrs[ring].ps);

		    /* Recurse. */
		    tr_p_mp_views_recurse(a_tr,
					  tr_p_ring_other_get(a_tr, ring),
					  a_ps,
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
		cw_tr_ring_t ring_a, ring_b;
		cw_tr_ring_t ring_a_other, ring_b_other;
		cw_tr_ps_t *ps_a, *ps_b;
		cw_tr_ps_t *ps_a_other, *ps_b_other;

		/* This is a normal trifurcating node.  This is the common case,
		 * and is handled separately from the code below for performance
		 * reasons. */

		/* Get all variables that are necessary for view calculation and
		 * recursion. */
		ring_a = CxmQriNext(a_tr->trrs, a_ring, link);
		ps_a = a_tr->trrs[ring_a].ps;
		ring_a_other = tr_p_ring_other_get(a_tr, ring_a);
		ps_a_other = a_tr->trrs[ring_a_other].ps;

		ring_b = CxmQriNext(a_tr->trrs, ring_a, link);
		ps_b = a_tr->trrs[ring_b].ps;
		ring_b_other = tr_p_ring_other_get(a_tr, ring_b);
		ps_b_other = a_tr->trrs[ring_b_other].ps;

		/* Calculate views and edges, and recurse. */
		tr_p_mp_pscore(a_tr, ps_a, a_ps, ps_b_other);
		tr_p_mp_pscore(a_tr,
			       a_tr->tres[tr_p_ring_edge_get(a_tr, ring_a)].ps,
			       ps_a,
			       a_tr->trrs[ring_a_other].ps);
		tr_p_mp_views_recurse(a_tr, ring_a_other, ps_a, a_bisect);

		tr_p_mp_pscore(a_tr, ps_b, a_ps, ps_a_other);
		tr_p_mp_pscore(a_tr,
			       a_tr->tres[tr_p_ring_edge_get(a_tr, ring_b)].ps,
			       ps_b,
			       a_tr->trrs[ring_b_other].ps);
		tr_p_mp_views_recurse(a_tr, ring_b_other, ps_b, a_bisect);

		break;
	    }
	    /* Fall through if this node is adjacent to the bisection. */
	}
	default:
	{
	    /* This is a multifurcating node. */
	    cxmError("XXX Not implemented");
	}
    }
}

/* Calculate the partial score for each edge in a_edges.  a_edges[0] must either
 * be CxmTrEdgeNone, or the edge connected to the node that is in turn
 * connected to the bisection edge. */
CxmpInline bool
tr_p_bisection_edge_list_mp(CxtTr *a_tr, CxtTrEdge *a_edges,
			    uint32_t a_nedges, CxtTrEdge a_bisect,
			    uint32_t a_maxscore)
{
    bool retval;

    if (a_edges[0] != CxmTrEdgeNone)
    {
	cw_tr_ring_t ring_a, ring_b;
	cw_tr_ps_t *ps, *ps_a, *ps_b;

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
	ps = a_tr->tres[a_edges[0]].ps;
	tr_p_mp_pscore(a_tr, ps, ps_a, ps_b);
	if (ps->subtrees_score + ps->node_score > a_maxscore)
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
	tr_p_mp_views_recurse(a_tr, ring_a, ps_b, a_bisect);
	tr_p_mp_views_recurse(a_tr, ring_b, ps_a, a_bisect);

#ifdef CxmDebug
	/* Validate per-edge partial scores. */
	{
	    uint32_t i;

	    for (i = 1; i < a_nedges; i++)
	    {
		/* All edge partial scores should have the same value, since the
		 * location of the root is irrelevant to the score. */
		if (a_tr->tres[a_edges[i]].ps->subtrees_score
		    + a_tr->tres[a_edges[i]].ps->node_score
		    != a_tr->tres[a_edges[0]].ps->subtrees_score
		    + a_tr->tres[a_edges[0]].ps->node_score)
		{
		    fprintf(stderr,
			    "%s:%d:%s(): Expected %u (%u + %u),"
			    " got %u (%u + %u)\n",
			    __FILE__, __LINE__, __func__,
			    a_tr->tres[a_edges[0]].ps->subtrees_score
			    + a_tr->tres[a_edges[0]].ps->node_score,
			    a_tr->tres[a_edges[0]].ps->subtrees_score,
			    + a_tr->tres[a_edges[0]].ps->node_score,
			    a_tr->tres[a_edges[i]].ps->subtrees_score
			    + a_tr->tres[a_edges[i]].ps->node_score,
			    a_tr->tres[a_edges[i]].ps->subtrees_score,
			    a_tr->tres[a_edges[i]].ps->node_score);
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
CxmpInline bool
tr_p_hold(CxtTr *a_tr, uint32_t a_max_hold, uint32_t a_neighbor,
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
	    a_tr->held = (cw_trh_t *) CxmMalloc(sizeof(cw_trh_t));
	    a_tr->heldlen = 1;
	}
	else if (a_tr->nheld == a_tr->heldlen)
	{
	    /* Reallocate. */
	    a_tr->held = (cw_trh_t *) CxmRealloc(a_tr->held,
						 sizeof(cw_trh_t)
						 * a_tr->heldlen * 2);
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
CxmpInline void
tr_p_tbr_neighbors_mp(CxtTr *a_tr, uint32_t a_max_hold,
		      uint32_t a_maxscore, cw_tr_hold_how_t a_how)
{
    uint32_t neighbor, i, j, k, curmax, score;
    CxtTrEdge bisect, edge_a, edge_b;
    CxtTrNode node_a, node_b;
    cw_tr_ps_t *ps_a, *ps_b;

    cxmDassert(tr_p_validate(a_tr));

    curmax = a_maxscore;

    /* Set up tree holding data structures. */
    a_tr->nheld = 0;

    /* Iteratively (logically) bisect at each edge in the tree. */
    neighbor = 0;
    for (i = 0; i < a_tr->nedges; i++)
    {
	bisect = a_tr->trt[i].bisect_edge;

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
	    if (edge_a != CxmTrEdgeNone)
	    {
		ps_a = a_tr->tres[edge_a].ps;
	    }
	    else
	    {
		ps_a = a_tr->trrs[CxmQliFirst(&a_tr->trns[node_a].rings)].ps;
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
		if (edge_b != CxmTrEdgeNone)
		{
		    ps_b = a_tr->tres[edge_b].ps;
		}
		else
		{
		    ps_b = a_tr->trrs[CxmQliFirst(&a_tr->trns[node_b].rings)].ps;
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
			cxmNotReached();
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

CxtTr *
CxTrNew(void)
{
    CxtTr *retval;

    retval = (CxtTr *) CxmMalloc(sizeof(CxtTr));
    tr_p_new(retval);

    return retval;
}

void
CxTrDelete(CxtTr *a_tr)
{
    cxmCheckPtr(a_tr);
    cxmAssert(a_tr->magic == CW_TR_MAGIC);

    if (a_tr->held != NULL)
    {
	CxmFree(a_tr->held);
    }

    /* This assumes that all nodes are deallocated before CxTrDelete() is
     * called. */
    if (a_tr->trns != NULL)
    {
	CxmFree(a_tr->trns);
    }

    if (a_tr->trt != NULL)
    {
	CxmFree(a_tr->trt);
    }

    if (a_tr->bedges != NULL)
    {
	CxmFree(a_tr->bedges);
    }

    /* This assumes that all edges are deallocated before CxTrDelete() is
     * called. */
    if (a_tr->tres != NULL)
    {
	CxmFree(a_tr->tres);
	CxmFree(a_tr->trrs);
    }

    CxmFree(a_tr);
}

uint32_t
CxTrNtaxaGet(CxtTr *a_tr)
{
    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));

    return a_tr->ntaxa;
}

uint32_t
CxTrNedgesGet(CxtTr *a_tr)
{
    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));

    return a_tr->nedges;
}

CxtTrNode
CxTrBaseGet(CxtTr *a_tr)
{
    cxmCheckPtr(a_tr);
    cxmAssert(a_tr->magic == CW_TR_MAGIC);

    return a_tr->base;
}

void
CxTrBaseSet(CxtTr *a_tr, CxtTrNode a_base)
{
    cxmCheckPtr(a_tr);
    cxmAssert(a_tr->magic == CW_TR_MAGIC);
#ifdef CxmDebug
    if (a_base != CxmTrNodeNone)
    {
	cxmDassert(tr_p_node_validate(a_tr, a_base));
    }
#endif

    a_tr->base = a_base;

    a_tr->modified = true;
}

void
CxTrCanonize(CxtTr *a_tr)
{
    /* Update internal state, so that ntaxa and nedges are correct. */
    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));

    if (a_tr->base != CxmTrNodeNone)
    {
	uint32_t ntaxa, nedges;
	cw_tr_ring_t ring;

	/* Set base to be the lowest-numbered taxon. */
	ntaxa = 0;
	nedges = 0;
	a_tr->base = tr_p_lowest(a_tr, a_tr->base, &ntaxa, &nedges);

	/* Get base's ring. */
	ring = CxmQliFirst(&a_tr->trns[a_tr->base].rings);
	if (ring != CW_TR_RING_NONE)
	{
	    /* Canonize the tree. */
	    tr_p_canonize(a_tr, tr_p_ring_other_get(a_tr, ring));
	}
    }

    /* Re-update internal state. */
    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));
}

void
CxTrTbr(CxtTr *a_tr, CxtTrEdge a_bisect, CxtTrEdge a_reconnect_a,
       CxtTrEdge a_reconnect_b)
{
    CxtTrNode node_a, node_b, nodes[4];
    CxtTrEdge tedges[3];
    uint32_t ntedges = 0;
    CxtTrNode tnodes[2];
    uint32_t ntnodes = 0;

    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));

    /* Get the nodes to either side of the edge where the bisection will be
     * done. */
    node_a = CxTrEdgeNodeGet(a_tr, a_bisect, 0);
    node_b = CxTrEdgeNodeGet(a_tr, a_bisect, 1);

    /* Bisect.  a_bisect will be used below for reconnection. */
    CxTrEdgeDetach(a_tr, a_bisect);
    
    /* For nodes_[ab], extract the node if it has only two neighbors.
     *
     * nodes_[0..1] are CxmTrNodeNone, unless they refer to the only node in a
     * subtree. */
    nodes[0] = tr_p_tbr_node_extract(a_tr, node_a, a_reconnect_a, a_reconnect_b,
				     tedges, &ntedges, tnodes, &ntnodes);
    nodes[1] = tr_p_tbr_node_extract(a_tr, node_b, a_reconnect_a, a_reconnect_b,
				     tedges, &ntedges, tnodes, &ntnodes);

    /* For each reconnection edge, splice a node into the edge (if the subtree
     * has more than one node).
     *
     * nodes[2..3] are set to CxmTrNodeNone if no reconnection edge is
     * specified. */
    if (a_reconnect_a != CxmTrEdgeNone)
    {
	nodes[2] = tr_p_tbr_node_splice(a_tr, a_reconnect_a,
					tedges, &ntedges, tnodes, &ntnodes);
    }
    else
    {
	nodes[2] = CxmTrNodeNone;
    }

    if (a_reconnect_b != CxmTrEdgeNone)
    {
	nodes[3] = tr_p_tbr_node_splice(a_tr, a_reconnect_b,
					tedges, &ntedges, tnodes, &ntnodes);
    }
    else
    {
	nodes[3] = CxmTrNodeNone;
    }

    /* If either subtree has only a single node, special care must be taken
     * during reconnection to re-associate the proper end of the bisection edge
     * with the single node.  This is because character state information for
     * leaf nodes is actually stored in the rings that attach them to the tree,
     * and breaking this association would require re-initializing the character
     * state set vectors. */
    if (nodes[0] != CxmTrNodeNone)
    {
	/* nodes[0] (same as node_a) is a single node. */
	cxmAssert(nodes[0] == node_a);
	cxmAssert(nodes[1] == CxmTrNodeNone);

	if (nodes[2] == CxmTrEdgeNone)
	{
	    nodes[2] = nodes[3];
	}
	CxTrEdgeAttach(a_tr, a_bisect, nodes[0], nodes[2]);
    }
    else if (nodes[1] != CxmTrNodeNone)
    {
	/* nodes[1] (same as node_b) is a single node. */
	cxmAssert(nodes[1] == node_b);

	if (nodes[2] == CxmTrEdgeNone)
	{
	    nodes[2] = nodes[3];
	}
	CxTrEdgeAttach(a_tr, a_bisect, nodes[2], nodes[1]);
    }
    else
    {
	/* Bisection was not done adjacent to a leaf node.  Attach the two
	 * spliced-in nodes. */
	CxTrEdgeAttach(a_tr, a_bisect, nodes[2], nodes[3]);
    }

    /* Update. */
    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));
}

uint32_t
CxTrTbrNneighborsGet(CxtTr *a_tr)
{
    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));

    return a_tr->trt[a_tr->trtused].offset;
}

void
CxTrTbrNeighborGet(CxtTr *a_tr, uint32_t a_neighbor,
		    CxtTrEdge *r_bisect, CxtTrEdge *r_reconnect_a,
		    CxtTrEdge *r_reconnect_b)
{
    cw_trt_t key, *trt;
    uint32_t rem;

    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));
    cxmAssert(a_neighbor < a_tr->trt[a_tr->trtused].offset);

    /* Get the bisection edge. */
    key.offset = a_neighbor;
    trt = bsearch(&key, a_tr->trt, a_tr->trtused, sizeof(cw_trt_t),
		  tr_p_trt_compare);
    cxmCheckPtr(trt);
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
CxTrAuxGet(CxtTr *a_tr)
{
    cxmCheckPtr(a_tr);
    cxmAssert(a_tr->magic == CW_TR_MAGIC);

    return a_tr->aux;
}

void
CxTrAuxSet(CxtTr *a_tr, void *a_aux)
{
    cxmCheckPtr(a_tr);
    cxmAssert(a_tr->magic == CW_TR_MAGIC);

    a_tr->aux = a_aux;
}

void
CxTrMpPrepare(CxtTr *a_tr, bool a_uninformative_eliminate,
	      char *a_taxa[], uint32_t a_ntaxa, uint32_t a_nchars)
{
    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));
    if (a_tr->base != CxmTrNodeNone)
    {
	cw_trn_t *trn;
	cw_tr_ring_t ring;
	cw_tre_t *tre;
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
			    cxmNotReached();
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
	CxmQliForeach(ring, &trn->rings, a_tr->trrs, link)
	{
	    /* Prepare edge before recursing. */
	    tre = &a_tr->tres[tr_p_ring_edge_get(a_tr, ring)];
	    if (tre->ps == NULL)
	    {
		tre->ps = tr_p_ps_new(a_tr);
	    }
	    tr_p_ps_prepare(a_tr, tre->ps, ninformative);

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
CxTrMpFinish(CxtTr *a_tr)
{

    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));

    if (a_tr->base != CxmTrNodeNone)
    {
	cw_trn_t *trn;
	cw_tr_ring_t ring;
	cw_tre_t *tre;

	/* Clean up the tree. */
	trn = &a_tr->trns[a_tr->base];
	CxmQliForeach(ring, &trn->rings, a_tr->trrs, link)
	{
	    /* Clean up edge before recursing. */
	    tre = &a_tr->tres[tr_p_ring_edge_get(a_tr, ring)];
	    if (tre->ps != NULL)
	    {
		tr_p_ps_delete(a_tr, tre->ps);
		tre->ps = NULL;
	    }

	    /* Clean up ring. */
	    tr_p_mp_ring_finish(a_tr, ring);

	    /* Recurse. */
	    tr_p_mp_finish_recurse(a_tr, tr_p_ring_other_get(a_tr, ring));
	}
    }
}

uint32_t
CxTrMpScore(CxtTr *a_tr)
{
    uint32_t retval;
    cw_tr_ring_t ring;

    cxmDassert(tr_p_validate(a_tr));

    if (a_tr->base != CxmTrNodeNone
	&& (ring = CxmQliFirst(&a_tr->trns[a_tr->base].rings)) != CW_TR_RING_NONE)
    {
	CxtTrEdge edge;
	cw_tr_ps_t *ps_a, *ps_b;

	edge = tr_p_ring_edge_get(a_tr, ring);

	/* Calculate partial scores for the subtrees on each end of edge. */
	ps_a = tr_p_mp_score_recurse(a_tr,
				     tr_p_edge_ring_get(a_tr, edge, 0),
				     CxmTrEdgeNone);
	ps_b = tr_p_mp_score_recurse(a_tr,
				     tr_p_edge_ring_get(a_tr, edge, 1),
				     CxmTrEdgeNone);

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
CxTrTbrBestNeighborsMp(CxtTr *a_tr, uint32_t a_max_hold)
{
    cxmDassert(tr_p_validate(a_tr));

    tr_p_tbr_neighbors_mp(a_tr, a_max_hold,
			  CxmTrMaxscoreNone,
			  TR_HOLD_BEST);
}

void
CxTrTbrBetterNeighborsMp(CxtTr *a_tr, uint32_t a_max_hold)
{
    uint32_t score;

    cxmDassert(tr_p_validate(a_tr));

    score = CxTrMpScore(a_tr);
    tr_p_tbr_neighbors_mp(a_tr, a_max_hold,
			  score > 0 ? score - 1 : 0,
			  TR_HOLD_BETTER);
}

void
CxTrTbrAllNeighborsMp(CxtTr *a_tr)
{
    cxmDassert(tr_p_validate(a_tr));

    tr_p_tbr_neighbors_mp(a_tr, CxmTrHoldAll,
			  CxmTrMaxscoreNone,
			  TR_HOLD_ALL);
}

void
CxTrHeldFinish(CxtTr *a_tr)
{
    cxmDassert(tr_p_validate(a_tr));

    if (a_tr->held != NULL)
    {
	CxmFree(a_tr->held);
	a_tr->held = NULL;
	a_tr->heldlen = 0;
	a_tr->nheld = 0;
    }
}

uint32_t
CxTrNheldGet(CxtTr *a_tr)
{
    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));

    return a_tr->nheld;
}

void
CxTrHeldGet(CxtTr *a_tr, uint32_t a_held, uint32_t *r_neighbor,
	    uint32_t *r_score)
{
    cw_trh_t *trh;

    tr_p_update(a_tr);
    cxmDassert(tr_p_validate(a_tr));
    cxmCheckPtr(a_tr->held);
    cxmAssert(a_held < a_tr->nheld);

    trh = &a_tr->held[a_held];
    *r_neighbor = trh->neighbor;
    *r_score = trh->score;
}
