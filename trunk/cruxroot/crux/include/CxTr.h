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

/******************************************************************************/

/* tr. */

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

/* Get the first edge in aNode's edge ring. */
void
CxTrNodeEdgeGet(CxtTr *aTr, CxtTrNode aNode,
		CxtTrEdge *rEdge, uint32_t *rEnd);

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

/* trEdge. */

/* Constructor. */
CxtTrEdge
CxTrEdgeNew(CxtTr *aTr);

/* Destructor. */
void
CxTrEdgeDelete(CxtTr *aTr, CxtTrEdge aEdge);

/* Get node aI of the edge (aI must be 0 or 1). */
CxtTrNode
CxTrEdgeNodeGet(CxtTr *aTr, CxtTrEdge aEdge, uint32_t aI);

/* Get the next edge in the edge ring at the aI end of aEdge. */
void
CxTrEdgeNextGet(CxtTr *aTr, CxtTrEdge aEdge, uint32_t aI,
		CxtTrEdge *rNext, uint32_t *rEnd);

/* Get the previous edge in the edge ring at the aI end of aEdge. */
void
CxTrEdgePrevGet(CxtTr *aTr, CxtTrEdge aEdge, uint32_t aI,
		CxtTrEdge *rPrev, uint32_t *rEnd);

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
