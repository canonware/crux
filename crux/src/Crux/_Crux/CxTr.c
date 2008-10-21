//==============================================================================
//
// <Copyright = jasone>
// <License>
//
//==============================================================================
//
// Version: Crux <Version = crux>
//
//==============================================================================
//
// Tr implements multifurcating phylogenetic trees and various operations on
// them.  Nodes can be manipulated via the CxTrNode*() APIs, and edges can be
// manipulated via the CxTrEdge*() APIs.
//
// Nodes and edges are allocated from per-tr internal arrays.  This allows array
// indices to be used for links between nodes and edges.
//
// Each edge has two ring objects associated with it, one for each end of the
// edge.  These ring objects are also allocated from an internal array that has
// a direct correspondence to the edge array.  This provides for constant-time
// conversion between edges and rings.
//
// Each node optionally references a ring, which contains all the edges that
// connect the node to other nodes.  Each ring object in the ring refers to the
// node, which allows constant-time conversion from rings to nodes.
//
//==============================================================================

#define CxmTr_c
#include "Crux/_cruxmodule.h"

//==============================================================================

// CxTrRing.

CxmpInline void
CxpTrRingInit(CxtTr *aTr, CxtTrRing aRing)
{
    CxtTrr *trr;

    CxmAssert((aRing >> 1) < aTr->ntres);

    trr = &aTr->trrs[aRing];

    CxmQriNew(aTr->trrs, aRing, link);
    trr->aux = NULL;
    trr->node = CxmTrNodeNone;
    trr->ps = NULL;
}

void *
CxTrRingAuxGet(CxtTr *aTr, CxtTrRing aRing)
{
    CxmDassert(CxpTrRingValidate(aTr, aRing));

    return aTr->trrs[aRing].aux;
}

void
CxTrRingAuxSet(CxtTr *aTr, CxtTrRing aRing, void *aAux)
{
    CxmDassert(CxpTrRingValidate(aTr, aRing));

    aTr->trrs[aRing].aux = aAux;
}

//==============================================================================

// Validation functions.

#ifdef CxmDebug
// Validate a ring.
bool
CxpTrRingValidate(CxtTr *aTr, CxtTrRing aRing)
{
    CxtTrr *trr;

    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
    CxmAssert((aRing >> 1) < aTr->ntres);

    trr = &aTr->trrs[aRing];

    CxmAssert(trr->node < aTr->ntrns || trr->node == CxmTrNodeNone);

    return true;
}

// Validate an edge.
bool
CxpTrEdgeValidate(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxtTre *tre;

    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
    CxmAssert(aEdge < aTr->ntres);

    tre = &aTr->tres[aEdge];

    CxmAssert(tre->magic == CxmTreMagic);

    CxpTrRingValidate(aTr, (aEdge << 1));
    CxpTrRingValidate(aTr, (aEdge << 1) + 1);

    return true;
}

// Validate a node.
bool
CxpTrNodeValidate(CxtTr *aTr, CxtTrNode aNode)
{
    CxtTrn *trn;
    CxtTrRing ring;
    uint32_t nneighbors;

    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
    CxmAssert(aNode < aTr->ntrns);

    trn = &aTr->trns[aNode];

    CxmAssert(trn->magic == CxmTrnMagic);

    nneighbors = 0;
    CxmQliForeach(ring, &trn->rings, aTr->trrs, link)
	{
	    // Validate edge.
	    CxpTrEdgeValidate(aTr, CxTrRingEdgeGet(aTr, ring));

	    nneighbors++;
	}

    if (trn->taxonNum != CxmTrNodeTaxonNone)
    {
	// Only leaf nodes can have taxon numbers.  Leaf nodes have at most 1
	// neighbor.
	CxmAssert(nneighbors <= 1);
    }

    return true;
}
#endif

//==============================================================================

// CxTrPs.

CxmpInline void
CxpTrPsDelete(CxtTr *aTr, CxtTrPs *aPs)
{
    if (aPs->chars != NULL)
    {
	CxmFree(aPs->achars);
    }

    CxmFree(aPs);
}

//==============================================================================

// trEdge.

CxmpInline void
CxpTrEdgeInit(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxtTre *tre;

    tre = &aTr->tres[aEdge];

    tre->u.aux = NULL;
    CxpTrRingInit(aTr, (aEdge << 1));
    CxpTrRingInit(aTr, (aEdge << 1) + 1);
    tre->length = 0.0;
    tre->ps = NULL;

#ifdef CxmDebug
    tre->magic = CxmTreMagic;
#endif
}

CxmpInline CxtTrEdge
CxpTrEdgeAlloc(CxtTr *aTr)
{
    CxtTrEdge rVal;

    if (aTr->sparetres == CxmTrEdgeNone)
    {
	uint32_t i, nspares;

	if (aTr->tres == NULL)
	{
	    aTr->tres = (CxtTre *) CxmMalloc(sizeof(CxtTre));
	    CxmAssert(aTr->trrs == NULL);
	    aTr->trrs = (CxtTrr *) CxmMalloc(sizeof(CxtTrr) * 2);
	    nspares = 1;
	    aTr->ntres = 1;
	}
	else
	{
	    aTr->tres = (CxtTre *) CxmRealloc(aTr->tres,
					      sizeof(CxtTre)
					      * aTr->ntres * 2);
	    CxmCheckPtr(aTr->trrs);
	    aTr->trrs = (CxtTrr *) CxmRealloc(aTr->trrs,
					      sizeof(CxtTrr)
					      * aTr->ntres * 4);
	    nspares = aTr->ntres;
	    aTr->ntres *= 2;
	}

	// Initialize last spare.
	aTr->sparetres = aTr->ntres - 1;
#ifdef CxmDebug
	aTr->tres[aTr->sparetres].magic = 0xa5a5a5a5;
#endif
	aTr->tres[aTr->sparetres].u.link = CxmTrEdgeNone;

	// Insert other spares into spares stack.
	for (i = 1; i < nspares; i++)
	{
	    aTr->sparetres--;
#ifdef CxmDebug
	    aTr->tres[aTr->sparetres].magic = 0xa5a5a5a5;
#endif
	    aTr->tres[aTr->sparetres].u.link = aTr->sparetres + 1;
	}
    }

    // Remove a spare from the spares stack.
    rVal = aTr->sparetres;
    aTr->sparetres = aTr->tres[rVal].u.link;

    // Initialize rVal.
    CxpTrEdgeInit(aTr, rVal);

    return rVal;
}

CxmpInline void
CxpTrEdgeDealloc(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxtTre *tre;
    CxtTrr *trr;

    tre = &aTr->tres[aEdge];
    if (tre->ps != NULL)
    {
	CxpTrPsDelete(aTr, tre->ps);
    }
#ifdef CxmDebug
    memset(tre, 0x5a, sizeof(CxtTre));
#endif

    trr = &aTr->trrs[aEdge << 1];
    if (trr->ps != NULL)
    {
	CxpTrPsDelete(aTr, trr->ps);
    }
#ifdef CxmDebug
    memset(trr, 0x5a, sizeof(CxtTrr));
#endif

    trr = &aTr->trrs[(aEdge << 1) + 1];
    if (trr->ps != NULL)
    {
	CxpTrPsDelete(aTr, trr->ps);
    }
#ifdef CxmDebug
    memset(trr, 0x5a, sizeof(CxtTrr));
#endif

    aTr->tres[aEdge].u.link = aTr->sparetres;
    aTr->sparetres = aEdge;
}

CxtTrEdge
CxTrEdgeNew(CxtTr *aTr)
{
    return CxpTrEdgeAlloc(aTr);
}

void
CxTrEdgeDelete(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
    CxmAssert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 0))
	      == CxmTrNodeNone);
    CxmAssert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 1))
	      == CxmTrNodeNone);

    CxpTrEdgeDealloc(aTr, aEdge);
}

double
CxTrEdgeLengthGet(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));

    return aTr->tres[aEdge].length;
}

void
CxTrEdgeLengthSet(CxtTr *aTr, CxtTrEdge aEdge, double aLength)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));

    aTr->tres[aEdge].length = aLength;
}

void *
CxTrEdgeAuxGet(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));

    return aTr->tres[aEdge].u.aux;
}

void
CxTrEdgeAuxSet(CxtTr *aTr, CxtTrEdge aEdge, void *aAux)
{
    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));

    aTr->tres[aEdge].u.aux = aAux;
}

void
CxTrEdgeAttach(CxtTr *aTr, CxtTrEdge aEdge, CxtTrNode aNodeA,
	       CxtTrNode aNodeB)
{
    CxtTrn *trn;
    CxtTrRing ring;

    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
    CxmDassert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 0))
	       == CxmTrNodeNone);
    CxmDassert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 1))
	       == CxmTrNodeNone);
    CxmDassert(CxpTrNodeValidate(aTr, aNodeA));
    CxmDassert(CxpTrNodeValidate(aTr, aNodeB));
    CxmAssert(aNodeA != aNodeB);
    CxmAssert(CxTrNodeDistance(aTr, aNodeA, aNodeB) == 0);

    // First end.
    ring = CxTrEdgeRingGet(aTr, aEdge, 0);
    trn = &aTr->trns[aNodeA];
    CxmQliHeadInsert(&trn->rings, aTr->trrs, ring, link);
    aTr->trrs[ring].node = aNodeA;

    // Second end.
    ring = CxTrEdgeRingGet(aTr, aEdge, 1);
    trn = &aTr->trns[aNodeB];
    CxmQliHeadInsert(&trn->rings, aTr->trrs, ring, link);
    aTr->trrs[ring].node = aNodeB;

    // Mark tree as modified.
    aTr->modified = true;

    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
    CxmDassert(CxpTrNodeValidate(aTr, aNodeA));
    CxmDassert(CxpTrNodeValidate(aTr, aNodeB));
    CxmAssert(CxTrNodeDistance(aTr, aNodeA, aNodeB) == 1);
}

void
CxTrEdgeDetach(CxtTr *aTr, CxtTrEdge aEdge)
{
    CxtTrRing ring;
    CxtTrn *trn;

    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
    CxmDassert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 0))
	       != CxmTrNodeNone);
    CxmDassert(CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge, 1))
	       != CxmTrNodeNone);
    CxmAssert(CxTrNodeDistance(aTr,
			       CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge,
								    0)),
			       CxTrRingNodeGet(aTr, CxTrEdgeRingGet(aTr, aEdge,
								    1))) == 1);

    // Detach from neighboring nodes.  Use CxmQliRemove() to make sure that the
    // nodes still point to their rings.
    ring = CxTrEdgeRingGet(aTr, aEdge, 0);
    trn = &aTr->trns[CxTrRingNodeGet(aTr, ring)];
    CxmQliRemove(&trn->rings, aTr->trrs, ring, link);
    aTr->trrs[ring].node = CxmTrNodeNone;

    ring = CxTrEdgeRingGet(aTr, aEdge, 1);
    trn = &aTr->trns[CxTrRingNodeGet(aTr, ring)];
    CxmQliRemove(&trn->rings, aTr->trrs, ring, link);
    aTr->trrs[ring].node = CxmTrNodeNone;

    // Mark tree as modified.
    aTr->modified = true;

    CxmDassert(CxpTrEdgeValidate(aTr, aEdge));
}

//==============================================================================

// CxTrNode.

CxmpInline void
CxpTrNodeInit(CxtTr *aTr, CxtTrNode aNode)
{
    CxtTrn *trn;

    trn = &aTr->trns[aNode];

    trn->u.aux = NULL;
    trn->taxonNum = CxmTrNodeTaxonNone;
    CxmQliNew(&trn->rings);

#ifdef CxmDebug
    trn->magic = CxmTrnMagic;
#endif
}

CxmpInline CxtTrNode
CxpTrNodeAlloc(CxtTr *aTr)
{
    CxtTrNode rVal;

    if (aTr->sparetrns == CxmTrNodeNone)
    {
	uint32_t i, nspares;

	// Allocate spares.
	if (aTr->trns == NULL)
	{
	    aTr->trns
		= (CxtTrn *) CxmMalloc(sizeof(CxtTrn));
	    nspares = 1;
	    aTr->ntrns = 1;
	}
	else
	{
	    aTr->trns = (CxtTrn *) CxmRealloc(aTr->trns,
					      sizeof(CxtTrn)
					      * aTr->ntrns * 2);
	    nspares = aTr->ntrns;
	    aTr->ntrns *= 2;
	}

	// Initialize last spare.
	aTr->sparetrns = aTr->ntrns - 1;
#ifdef CxmDebug
	aTr->trns[aTr->sparetrns].magic = 0xa5a5a5a5;
#endif
	aTr->trns[aTr->sparetrns].u.link = CxmTrNodeNone;

	// Insert other spares into spares stack.
	for (i = 1; i < nspares; i++)
	{
	    aTr->sparetrns--;
#ifdef CxmDebug
	    aTr->trns[aTr->sparetrns].magic = 0xa5a5a5a5;
#endif
	    aTr->trns[aTr->sparetrns].u.link = aTr->sparetrns + 1;
	}
    }

    // Remove a spare from the spares stack.
    rVal = aTr->sparetrns;
    aTr->sparetrns = aTr->trns[rVal].u.link;

    // Initialize rVal.
    CxpTrNodeInit(aTr, rVal);

    return rVal;
}

CxmpInline void
CxpTrNodeDealloc(CxtTr *aTr, CxtTrNode aNode)
{
#ifdef CxmDebug
    CxtTrn *trn;

    trn = &aTr->trns[aNode];

    memset(trn, 0x5a, sizeof(CxtTrn));
#endif

    aTr->trns[aNode].u.link = aTr->sparetrns;
    aTr->sparetrns = aNode;
}

// Calculate the number of edges connected to the node that aRing is connected
// to.
CxmpInline uint32_t
CxpTrNodeDegree(CxtTr *aTr, CxtTrRing aRing)
{
    uint32_t rVal;
    CxtTrRing ring;

    rVal = 1;
    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    rVal++;
	}

    return rVal;
}

// Calculate the number of edges between two nodes.  A distance of 0 means that
// there is no path between the two nodes.
static uint32_t
CxpTrNodeDistance(CxtTr *aTr, CxtTrRing aRing, CxtTrNode aOther,
		  uint32_t aDistance)
{
    uint32_t rVal;
    CxtTrRing ring;

    if (CxTrRingNodeGet(aTr, aRing) == aOther)
    {
	rVal = aDistance;
	goto RETURN;
    }

    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    if ((rVal = CxpTrNodeDistance(aTr, CxTrRingOtherGet(aTr, ring),
					  aOther, aDistance + 1)) != 0)
	    {
		goto RETURN;
	    }
	}

    rVal = 0;
    RETURN:
    return rVal;
}

CxtTrNode
CxTrNodeNew(CxtTr *aTr)
{
    return CxpTrNodeAlloc(aTr);
}

void
CxTrNodeDelete(CxtTr *aTr, CxtTrNode aNode)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));
    CxmAssert(CxmQliFirst(&aTr->trns[aNode].rings) == CxmTrRingNone);

    CxpTrNodeDealloc(aTr, aNode);
}

uint32_t
CxTrNodeTaxonNumGet(CxtTr *aTr, CxtTrNode aNode)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    return aTr->trns[aNode].taxonNum;
}

void
CxTrNodeTaxonNumSet(CxtTr *aTr, CxtTrNode aNode,
		    uint32_t aTaxonNum)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    aTr->trns[aNode].taxonNum = aTaxonNum;

    aTr->modified = true;
}

CxtTrRing
CxTrNodeRingGet(CxtTr *aTr, CxtTrNode aNode)
{
    CxtTrRing rVal;
    CxtTrn *trn;

    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    trn = &aTr->trns[aNode];

    rVal = CxmQliFirst(&trn->rings);

    return rVal;
}

void *
CxTrNodeAuxGet(CxtTr *aTr, CxtTrNode aNode)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    return aTr->trns[aNode].u.aux;
}

void
CxTrNodeAuxSet(CxtTr *aTr, CxtTrNode aNode, void *aAux)
{
    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    aTr->trns[aNode].u.aux = aAux;
}

uint32_t
CxTrNodeDegree(CxtTr *aTr, CxtTrNode aNode)
{
    uint32_t rVal;
    CxtTrRing ring;

    CxmDassert(CxpTrNodeValidate(aTr, aNode));

    ring = CxmQliFirst(&aTr->trns[aNode].rings);
    if (ring != CxmTrRingNone)
    {
	rVal = CxpTrNodeDegree(aTr, ring);
    }
    else
    {
	rVal = 0;
    }

    return rVal;
}

// Set the head of the ring for aNode to aRing.
void
CxTrNodeRingSet(CxtTr *aTr, CxtTrNode aNode, CxtTrRing aRing)
{
    CxmAssert(CxTrRingNodeGet(aTr, aRing) == aNode);

    CxmQliFirst(&aTr->trns[aNode].rings) = aRing;
}

uint32_t
CxTrNodeDistance(CxtTr *aTr, CxtTrNode aNode, CxtTrNode aOther)
{
    uint32_t rVal;
    CxtTrn *trn;
    CxtTrRing ring;

    CxmDassert(CxpTrNodeValidate(aTr, aNode));
    CxmDassert(CxpTrNodeValidate(aTr, aOther));
    CxmAssert(aNode != aOther);

    trn = &aTr->trns[aNode];

    ring = CxmQliFirst(&trn->rings);
    if (ring != CxmTrRingNone)
    {
	CxmQliForeach(ring, &trn->rings, aTr->trrs, link)
	    {
		if ((rVal = CxpTrNodeDistance(aTr, CxTrRingOtherGet(aTr, ring),
					      aOther, 1)) != 0)
		{
		    break;
		}
	    }
    }
    else
    {
	rVal = 0;
    }

    return rVal;
}

//==============================================================================

// tr.

// Initialize everything except trns and sparetrns.
CxmpInline void
CxpTrNew(CxtTr *aTr)
{
    aTr->aux = NULL;
    aTr->modified = false;
    aTr->base = CxmTrNodeNone;
    aTr->ntaxa = 0;
    aTr->nedges = 0;
    aTr->trns = NULL;
    aTr->ntrns = 0;
    aTr->sparetrns = CxmTrNodeNone;
    aTr->tres = NULL;
    aTr->ntres = 0;
    aTr->sparetres = CxmTrEdgeNone;
    aTr->trrs = NULL;

#ifdef CxmDebug
    aTr->magic = CxmTrMagic;
#endif
}

// Recursively traverse the tree, count the number of taxa, and find the lowest
// numbered taxon.
static CxtTrNode
CxpTrLowestRecurse(CxtTr *aTr, CxtTrRing aRing, uint32_t *rNtaxa,
		   uint32_t *rNedges, CxtTrNode aRoot)
{
    CxtTrNode rVal, node, root, troot;
    CxtTrRing ring;
    CxtTrn *trn;

    node = CxTrRingNodeGet(aTr, aRing);
    trn = &aTr->trns[node];

    if (trn->taxonNum != CxmTrNodeTaxonNone)
    {
	// Leaf node.
	(*rNtaxa)++;
    }

    if (trn->taxonNum != CxmTrNodeTaxonNone
	&& (aRoot == CxmTrNodeNone
	    || trn->taxonNum < aTr->trns[aRoot].taxonNum))
    {
	rVal = node;
	root = node;
    }
    else
    {
	rVal = CxmTrNodeNone;
	root = aRoot;
    }

    // Iterate over neighbors.
    CxmQriOthersForeach(ring, aTr->trrs, aRing, link)
	{
	    // Count edge.
	    (*rNedges)++;

	    troot = CxpTrLowestRecurse(aTr, CxTrRingOtherGet(aTr, ring),
				       rNtaxa, rNedges, root);
	    if (troot != CxmTrNodeNone)
	    {
		rVal = troot;
		root = troot;
	    }
	}

    return rVal;
}

// Recursively traverse the tree, count the number of taxa, and find the lowest
// numbered taxon.
static CxtTrNode
CxpTrLowest(CxtTr *aTr, CxtTrNode aNode, uint32_t *rNtaxa,
	    uint32_t *rNedges)
{
    CxtTrNode rVal, root, troot;
    CxtTrRing ring;
    CxtTrn *trn;

    CxmAssert(aNode != CxmTrNodeNone);

    trn = &aTr->trns[aNode];

    if (trn->taxonNum != CxmTrNodeTaxonNone)
    {
	// Leaf node.
	(*rNtaxa)++;
    }

    if (trn->taxonNum != CxmTrNodeTaxonNone)
    {
	rVal = aNode;
	root = aNode;
    }
    else
    {
	rVal = CxmTrNodeNone;
	root = CxmTrNodeNone;
    }

    // Iterate over neighbors.
    CxmQliForeach(ring, &trn->rings, aTr->trrs, link)
	{
	    // Count edge.
	    (*rNedges)++;

	    troot = CxpTrLowestRecurse(aTr, CxTrRingOtherGet(aTr, ring),
				       rNtaxa, rNedges, root);
	    if (troot != CxmTrNodeNone)
	    {
		rVal = troot;
		root = troot;
	    }
	}

    return rVal;
}

#ifdef CxmDebug
// Validate a tree.
bool
CxpTrValidate(CxtTr *aTr)
{
    uint32_t i, ntaxa, nedges;

    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
    CxmAssert(aTr->modified == false);

    ntaxa = 0;
    nedges = 0;
    if (aTr->base != CxmTrNodeNone)
    {
	CxpTrLowest(aTr, aTr->base, &ntaxa, &nedges);
    }
    CxmAssert(aTr->ntaxa == ntaxa);
    CxmAssert(aTr->nedges == nedges);

    // Iterate over trns and do some basic sanity checks.
    for (i = 0; i < aTr->ntrns; i++)
    {
	if (aTr->trns[i].magic == CxmTrnMagic)
	{
	    CxpTrNodeValidate(aTr, (CxtTrNode) i);
	}
	else
	{
	    // Make sure there are no valid trn's in the free list.
	    CxmAssert(aTr->trns[i].u.link == CxmTrNodeNone
		      || aTr->trns[aTr->trns[i].u.link].magic
		      != CxmTrnMagic);
	}
    }

    CxmAssert(aTr->sparetrns == CxmTrNodeNone
	      || aTr->trns[aTr->sparetrns].magic != CxmTrnMagic);

    return true;
}
#endif

static void
CxpTrNtaxaNedgesUpdate(CxtTr *aTr)
{
    uint32_t ntaxa, nedges;

    // Update ntaxa and nedges.
    ntaxa = 0;
    nedges = 0;
    if (aTr->base != CxmTrNodeNone)
    {
	CxpTrLowest(aTr, aTr->base, &ntaxa, &nedges);
    }

    aTr->ntaxa = ntaxa;
    aTr->nedges = nedges;
}

CxmpInline void
CxpTrUpdate(CxtTr *aTr)
{
    if (aTr->modified)
    {
	uint32_t nedgesPrev;

	// Store nedges before updating.
	nedgesPrev = aTr->nedges;

	// Update ntaxa and nedges.
	CxpTrNtaxaNedgesUpdate(aTr);

	// Reset the modified flag.
	aTr->modified = false;
    }
}

CxtTr *
CxTrNew(void)
{
    CxtTr *rVal;

    rVal = (CxtTr *) CxmMalloc(sizeof(CxtTr));
    CxpTrNew(rVal);

    return rVal;
}

void
CxTrDelete(CxtTr *aTr)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);

    // This assumes that all nodes are deallocated before CxTrDelete() is
    // called.
    if (aTr->trns != NULL)
    {
	CxmFree(aTr->trns);
    }

    // This assumes that all edges are deallocated before CxTrDelete() is
    // called.
    if (aTr->tres != NULL)
    {
	CxmFree(aTr->tres);
	CxmFree(aTr->trrs);
    }

    CxmFree(aTr);
}

uint32_t
CxTrNtaxaGet(CxtTr *aTr)
{
    CxpTrUpdate(aTr);
    CxmDassert(CxpTrValidate(aTr));

    return aTr->ntaxa;
}

uint32_t
CxTrNedgesGet(CxtTr *aTr)
{
    CxpTrUpdate(aTr);
    CxmDassert(CxpTrValidate(aTr));

    return aTr->nedges;
}

CxtTrNode
CxTrBaseGet(CxtTr *aTr)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);

    return aTr->base;
}

void
CxTrBaseSet(CxtTr *aTr, CxtTrNode aBase)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);
#ifdef CxmDebug
    if (aBase != CxmTrNodeNone)
    {
	CxmDassert(CxpTrNodeValidate(aTr, aBase));
    }
#endif

    aTr->base = aBase;

    aTr->modified = true;
}

void *
CxTrAuxGet(CxtTr *aTr)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);

    return aTr->aux;
}

void
CxTrAuxSet(CxtTr *aTr, void *aAux)
{
    CxmCheckPtr(aTr);
    CxmAssert(aTr->magic == CxmTrMagic);

    aTr->aux = aAux;
}
