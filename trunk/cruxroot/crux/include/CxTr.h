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

typedef struct CxsTr CxtTr;
#define CxmTrMaxscoreNone 0xffffffffU
#define CxmTrHoldAll 0xffffffffU

typedef uint32_t CxtTrNode;
#define CxmTrNodeNone 0xffffffffU
#define CxmTrNodeTaxonNone 0xffffffffU

typedef uint32_t CxtTrEdge;
#define CxmTrEdgeNone 0xffffffffU

typedef uint32_t CxtTrRing;
#define CxmTrRingNone UINT_MAX

typedef struct CxsTrPs CxtTrPs;
typedef struct CxsTrn CxtTrn;
typedef struct CxsTrr CxtTrr;
typedef struct CxsTre CxtTre;
typedef struct CxsTrt CxtTrt;
typedef struct CxsTrh CxtTrh;

/* Character (in the systematics sense of the word). */
typedef char CxtTrc;

/* Partial parsimony score information. */
struct CxsTrPs
{
    /* Parent which most recently used this node's partial score when caching
     * its results.  Both children must still point to the parent in order for
     * the cached results to be valid. */
    CxtTrPs *parent;

    /* Sum of the subtree scores, and this node's score, given particular
     * children.  In order for this to be useful, both childrens' parent
     * pointers must still point to this node. */
    uint32_t subtreesScore;
    uint32_t nodeScore;

    /* chars points to an array of Fitch parsimony state sets.  Each element in
     * the array contains a bitmap representation of a subset of {ACGT} in the 4
     * least significant bits.  T is the least significant bit.  1 means that a
     * nucleotide is in the set.
     *
     * There are nchars character state sets.
     *
     * achars is the actual allocation, which is padded in order to
     * be able to guarantee that chars is 16 byte-aligned. */
    CxtTrc *chars;
    uint32_t nchars;
    CxtTrc *achars;
};

/* Tree node for an unrooted bifurcating phylogenetic tree. */
struct CxsTrn
{
#ifdef CxmDebug
    uint32_t magic;
#define CxmTrnMagic 0x63329478
#endif

    union
    {
	/* Auxiliary opaque data pointer. */
	void *aux;

	/* Spares linkage. */
	CxtTrNode link;
    } u;

    /* If CxmTrNodeTaxonNone, then the node is not a leaf node. */
    uint32_t taxonNum;

    /* Ring of trr's, which are associated with tre's. */
    CxmQliHead rings;
};

/* Tree node edge ring element. */
struct CxsTrr
{
    /* Ring linkage. */
    CxmQri link;

    /* Auxiliary opaque data pointer. */
    void *aux;

    /* Node whose ring this trr is a part of. */
    CxtTrNode node;

    /* Used for Fitch parsimony scoring. */
    CxtTrPs *ps;
};

/* Tree edge information. */
struct CxsTre
{
#ifdef CxmDebug
    uint32_t magic;
#define CxmTreMagic 0xa683fa07
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
    CxtTrPs *ps;
};

/* TBR neighbor. */
struct CxsTrt
{
    /* Number of neighbors that can be reached by doing TBR at edges before this
     * one.  This is also the neighbor index of the first neighbor that can be
     * reached by doing TBR on this edge. */
    uint32_t offset;

    /* Bisection edge. */
    CxtTrEdge bisectEdge;
};

/* Held neighbor tree. */
struct CxsTrh
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
    CxeTrHoldBest,
    CxeTrHoldBetter,
    CxeTrHoldAll
} CxtTrHoldHow;

struct CxsTr
{
#ifdef CxmDebug
    uint32_t magic;
#define CxmTrMagic 0x39886394
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
     * The first list starts at offset 0 and has nbedgesA elements.  The second
     * list starts at offset nbedgesA and has nbedgesB elements. */
    CxtTrEdge *bedges;
    uint32_t nbedgesA;
    uint32_t nbedgesB;

    /* trt is an array of elements that store per-edge information that is used
     * for TBR-related functions.  There is one more element in trt than there
     * are edges in the tree.  This is critical to the way binary searching on
     * the array is done, and it also makes it easy to get the total number of
     * TBR neighbors this tree has (trt[nedges].offset).
     *
     * Only the first trtused elements are valid, since not all bisection edges
     * necessarily result in neighbors. */
    CxtTrt *trt;
    uint32_t trtused;

    /* Pointer to an array of trn's.  ntrns is the total number of trn's, not
     * all of which are necessarily in use.
     *
     * sparetrns is the index of the first spare trn in the spares stack. */
    CxtTrn *trns;
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
    CxtTre *tres;
    uint32_t ntres;
    CxtTrEdge sparetres;
    CxtTrr *trrs;

    /* held is an array of held TBR neighbors.  The array is iteratively doubled
     * as necessary.  heldlen is the actual length of the array, and nheld is
     * the number of elements in use. */
    CxtTrh *held;
    uint32_t heldlen;
    uint32_t nheld;
};

/******************************************************************************/

/* Validation function prototypes.  These are called by the debug versions of
 * some of the inline functions, so must be exposed here. */

#ifdef CxmDebug
bool
CxpTrRingValidate(CxtTr *aTr, CxtTrRing aRing);
bool
CxpTrEdgeValidate(CxtTr *aTr, CxtTrEdge aEdge);
bool
CxpTrNodeValidate(CxtTr *aTr, CxtTrNode aNode);
bool
CxpTrValidate(CxtTr *aTr);
#endif

/******************************************************************************/

/* CxTr. */

/* Constructor. */
CxtTr *
CxTrNew(void);

/* Destructor. */
void
CxTrDelete(CxtTr *aTr);

/* Get the number of taxa in the tree. */
uint32_t
CxTrNtaxaGet(CxtTr *aTr);

/* Get the number of edges in the tree. */
uint32_t
CxTrNedgesGet(CxtTr *aTr);

/* Get the base of the tree. */
CxtTrNode
CxTrBaseGet(CxtTr *aTr);

/* Set the base of the tree. */
void
CxTrBaseSet(CxtTr *aTr, CxtTrNode aBase);

/* Canonize the tree.  The base of the tree is set to the lowest numbered taxon,
 * and internal nodes are adjusted such that their edge rings are ordered
 * (subtrees with lower minimum taxon numbers come first), and the edge returned
 * by CxTrNodeEdgeGet() is the edge that leads back to the base. */
void
CxTrCanonize(CxtTr *aTr);

/* Perform TBR. */
void
CxTrTbr(CxtTr *aTr, CxtTrEdge aBisect, CxtTrEdge aReconnectA,
	CxtTrEdge aReconnectB);

/* Get the number of neighboring trees reachable via TBR. */
uint32_t
CxTrTbrNneighborsGet(CxtTr *aTr);

/* Get the parameters necessary to transorm this tree to neighbor aNeighbor. */
void
CxTrTbrNeighborGet(CxtTr *aTr, uint32_t aNeighbor,
		   CxtTrEdge *rBisect, CxtTrEdge *rReconnectA,
		   CxtTrEdge *rReconnectB);

/* Get the value of the auxiliary pointer associated with the tr. */
void *
CxTrAuxGet(CxtTr *aTr);

/* Set the value of the auxiliary pointer associated with the tr. */
void
CxTrAuxSet(CxtTr *aTr, void *aAux);

/* Prepare for calculating Fitch parsimony scores.  aTaxa points to an array of
 * character array pointers, where the index into aTaxa corresponds to taxon
 * number.  The character arrays need not be nil-terminated. */
void
CxTrMpPrepare(CxtTr *aTr, bool aUninformativeEliminate,
	      char *aTaxa[], uint32_t aNtaxa, uint32_t aNchars);

/* Clear the data structures used for calculating Fitch parsimony scores. */
void
CxTrMpFinish(CxtTr *aTr);

/* Calculate the Fitch parsimony score for this tree. */
uint32_t
CxTrMpScore(CxtTr *aTr);

/* Calculate the Fitch parsimony of all TBR neighbors, and keep track of up to
 * aMaxHold of the best neighbors (or all best neighbors, if aMaxHold is
 * CxmTrHoldAll). */
void
CxTrTbrBestNeighborsMp(CxtTr *aTr, uint32_t aMaxHold);

/* Calculate the Fitch parsimony of all TBR neighbors, and keep track of up to
 * aMaxHold of the better neighbors (or all better neighbors, if aMaxHold is
 * CxmTrHoldAll). */
void
CxTrTbrBetterNeighborsMp(CxtTr *aTr, uint32_t aMaxHold);

/* Calculate the Fitch parsimony of all TBR neighbors, and keep track of all
 * neighbors. */
void
CxTrTbrAllNeighborsMp(CxtTr *aTr);

/* Clear the data structures used to store held trees. */
void
CxTrHeldFinish(CxtTr *aTr);

/* Get the number of trees currently held. */
uint32_t
CxTrNheldGet(CxtTr *aTr);

/* Get the aHeld'th held tree, and its score.  *rNeighbor can be passed to
 * CxTrTbrNeighborGet() in order to get the TBR transformation parameters,
 * which can then be passed to CxTrTbr(). */
void
CxTrHeldGet(CxtTr *aTr, uint32_t aHeld, uint32_t *rNeighbor,
	    uint32_t *rScore);

/******************************************************************************/

/* CxTrNode. */

/* Constructor. */
CxtTrNode
CxTrNodeNew(CxtTr *aTr);

/* Destructor. */
void
CxTrNodeDelete(CxtTr *aTr, CxtTrNode aNode);

/* Get the taxon number associated with aNode.  Return CxmTrNodeTaxonNone if
 * no taxon number is set. */
uint32_t
CxTrNodeTaxonNumGet(CxtTr *aTr, CxtTrNode aNode);

/* Set the taxon number associated with aNode (use CxmTrNodeTaxonNone to
 * unset the taxon number. */
void
CxTrNodeTaxonNumSet(CxtTr *aTr, CxtTrNode aNode,
		    uint32_t aTaxonNum);

/* Get the first ring element in aNode's edge ring. */
CxtTrRing
CxTrNodeRingGet(CxtTr *aTr, CxtTrNode aNode);

/* Get the value of the auxiliary pointer associated with the node. */
void *
CxTrNodeAuxGet(CxtTr *aTr, CxtTrNode aNode);

/* Set the value of the auxiliary pointer associated with the node. */
void
CxTrNodeAuxSet(CxtTr *aTr, CxtTrNode aNode, void *aAux);

/* Get the degree (number of edges) of the node. */
uint32_t
CxTrNodeDegree(CxtTr *aTr, CxtTrNode aNode);

/* Get the number of edges between aNode and aOther.  0 means that there is
 * no path between the two nodes. */
uint32_t
CxTrNodeDistance(CxtTr *aTr, CxtTrNode aNode, CxtTrNode aOther);

/******************************************************************************/

/* CxTrEdge. */

/* Constructor. */
CxtTrEdge
CxTrEdgeNew(CxtTr *aTr);

/* Destructor. */
void
CxTrEdgeDelete(CxtTr *aTr, CxtTrEdge aEdge);

/* Get the edge length. */
double
CxTrEdgeLengthGet(CxtTr *aTr, CxtTrEdge aEdge);

/* Set the edge length. */
void
CxTrEdgeLengthSet(CxtTr *aTr, CxtTrEdge aEdge, double aLength);

/* Get the value of the auxiliary pointer associated with the edge. */
void *
CxTrEdgeAuxGet(CxtTr *aTr, CxtTrEdge aEdge);

/* Set the value of the auxiliary pointer associated with the edge. */
void
CxTrEdgeAuxSet(CxtTr *aTr, CxtTrEdge aEdge, void *aAux);

/* Attach edge to two nodes. */
void
CxTrEdgeAttach(CxtTr *aTr, CxtTrEdge aEdge, CxtTrNode aNodeA,
	       CxtTrNode aNodeB);

/* Detach edge from nodes. */
void
CxTrEdgeDetach(CxtTr *aTr, CxtTrEdge aEdge);

#if (!defined(CxmUseInlines))
CxtTrRing
CxTrEdgeRingGet(CxtTr *aTr, CxtTrEdge aEdge, uint32_t aEnd);
#endif
#if (defined(CxmUseInlines) || defined(CxmTr_c))
CxmInline CxtTrRing
CxTrEdgeRingGet(CxtTr *aTr, CxtTrEdge aEdge, uint32_t aEnd)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
    CxmAssert(aEnd == 0 || aEnd == 1);

    return ((aEdge << 1) + aEnd);
}
#endif

/******************************************************************************/

/* CxTrRing. */

#if (!defined(CxmUseInlines))
CxtTrNode
CxTrRingNodeGet(CxtTr *aTr, CxtTrRing aRing);
CxtTrEdge
CxTrRingEdgeGet(CxtTr *aTr, CxtTrRing aRing);
CxtTrRing
CxTrRingOtherGet(CxtTr *aTr, CxtTrRing aRing);
CxmInline CxtTrRing
CxTrRingNextGet(CxtTr *aTr, CxtTrRing aRing);
CxmInline CxtTrRing
CxTrRingPrevGet(CxtTr *aTr, CxtTrRing aRing);
#endif

#if (defined(CxmUseInlines) || defined(CxmTr_c))
/* Get the node associated with aRing. */
CxmInline CxtTrNode
CxTrRingNodeGet(CxtTr *aTr, CxtTrRing aRing)
{
    CxmDassert(CxpTrRingValidate(aTr, aRing));

    return aTr->trrs[aRing].node;
}

/* Get the edge that aRing is part of. */
CxmInline CxtTrEdge
CxTrRingEdgeGet(CxtTr *aTr, CxtTrRing aRing)
{
    CxmDassert(CxpTrRingValidate(aTr, aRing));

    return (aRing >> 1);
}

/* Get the ring element on the other end of the edge associated with aRing. */
CxmInline CxtTrRing
CxTrRingOtherGet(CxtTr *aTr, CxtTrRing aRing)
{
    CxmDassert(CxpTrRingValidate(aTr, aRing));

    return (aRing ^ 1);
}

/* Get the next element in the edge ring. */
CxmInline CxtTrRing
CxTrRingNextGet(CxtTr *aTr, CxtTrRing aRing)
{
    CxmDassert(CxpTrRingValidate(aTr, aRing));

    return CxmQriNext(aTr->trrs, aRing, link);
}

/* Get the next element in the edge ring. */
CxmInline CxtTrRing
CxTrRingPrevGet(CxtTr *aTr, CxtTrRing aRing)
{
    CxmDassert(CxpTrRingValidate(aTr, aRing));

    return CxmQriPrev(aTr->trrs, aRing, link);
}
#endif

void *
CxTrRingAuxGet(CxtTr *aTr, CxtTrRing aRing);
void
CxTrRingAuxSet(CxtTr *aTr, CxtTrRing aRing, void *aAux);
