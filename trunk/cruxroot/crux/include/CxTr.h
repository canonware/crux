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
CxTrDelete(CxtTr *a_tr);

/* Get the number of taxa in the tree. */
uint32_t
CxTrNtaxaGet(CxtTr *a_tr);

/* Get the number of edges in the tree. */
uint32_t
CxTrNedgesGet(CxtTr *a_tr);

/* Get the base of the tree. */
CxtTrNode
CxTrBaseGet(CxtTr *a_tr);

/* Set the base of the tree. */
void
CxTrBaseSet(CxtTr *a_tr, CxtTrNode a_base);

/* Canonize the tree.  The base of the tree is set to the lowest numbered taxon,
 * and internal nodes are adjusted such that their edge rings are ordered
 * (subtrees with lower minimum taxon numbers come first), and the edge returned
 * by CxTrNodeEdgeGet() is the edge that leads back to the base. */
void
CxTrCanonize(CxtTr *a_tr);

/* Perform TBR. */
void
CxTrTbr(CxtTr *a_tr, CxtTrEdge a_bisect, CxtTrEdge a_reconnect_a,
       CxtTrEdge a_reconnect_b);

/* Get the number of neighboring trees reachable via TBR. */
uint32_t
CxTrTbrNneighborsGet(CxtTr *a_tr);

/* Get the parameters necessary to transorm this tree to neighbor a_neighbor. */
void
CxTrTbrNeighborGet(CxtTr *a_tr, uint32_t a_neighbor,
		    CxtTrEdge *r_bisect, CxtTrEdge *r_reconnect_a,
		    CxtTrEdge *r_reconnect_b);

/* Get the value of the auxiliary pointer associated with the tr. */
void *
CxTrAuxGet(CxtTr *a_tr);

/* Set the value of the auxiliary pointer associated with the tr. */
void
CxTrAuxSet(CxtTr *a_tr, void *a_aux);

/* Prepare for calculating Fitch parsimony scores.  a_taxa points to an array of
 * character array pointers, where the index into a_taxa corresponds to taxon
 * number.  The character arrays need not be nil-terminated. */
void
CxTrMpPrepare(CxtTr *a_tr, bool a_uninformative_eliminate,
	      char *a_taxa[], uint32_t a_ntaxa, uint32_t a_nchars);

/* Clear the data structures used for calculating Fitch parsimony scores. */
void
CxTrMpFinish(CxtTr *a_tr);

/* Calculate the Fitch parsimony score for this tree. */
uint32_t
CxTrMpScore(CxtTr *a_tr);

/* Calculate the Fitch parsimony of all TBR neighbors, and keep track of up to
 * a_max_hold of the best neighbors (or all best neighbors, if a_max_hold is
 * CxmTrHoldAll). */
void
CxTrTbrBestNeighborsMp(CxtTr *a_tr, uint32_t a_max_hold);

/* Calculate the Fitch parsimony of all TBR neighbors, and keep track of up to
 * a_max_hold of the better neighbors (or all better neighbors, if a_max_hold is
 * CxmTrHoldAll). */
void
CxTrTbrBetterNeighborsMp(CxtTr *a_tr, uint32_t a_max_hold);

/* Calculate the Fitch parsimony of all TBR neighbors, and keep track of all
 * neighbors. */
void
CxTrTbrAllNeighborsMp(CxtTr *a_tr);

/* Clear the data structures used to store held trees. */
void
CxTrHeldFinish(CxtTr *a_tr);

/* Get the number of trees currently held. */
uint32_t
CxTrNheldGet(CxtTr *a_tr);

/* Get the a_held'th held tree, and its score.  *r_neighbor can be passed to
 * CxTrTbrNeighborGet() in order to get the TBR transformation parameters,
 * which can then be passed to CxTrTbr(). */
void
CxTrHeldGet(CxtTr *a_tr, uint32_t a_held, uint32_t *r_neighbor,
	    uint32_t *r_score);

/******************************************************************************/

/* tr_node. */

/* Constructor. */
CxtTrNode
CxTrNodeNew(CxtTr *a_tr);

/* Destructor. */
void
CxTrNodeDelete(CxtTr *a_tr, CxtTrNode a_node);

/* Get the taxon number associated with a_node.  Return CXTTRNODEAXON_NONE if
 * no taxon number is set. */
uint32_t
CxTrNodeTaxonNumGet(CxtTr *a_tr, CxtTrNode a_node);

/* Set the taxon number associated with a_node (use CXTTRNODEAXON_NONE to
 * unset the taxon number. */
void
CxTrNodeTaxonNumSet(CxtTr *a_tr, CxtTrNode a_node,
		      uint32_t a_taxon_num);

/* Get the first edge in a_node's edge ring. */
void
CxTrNodeEdgeGet(CxtTr *a_tr, CxtTrNode a_node,
		 CxtTrEdge *r_edge, uint32_t *r_end);

/* Get the value of the auxiliary pointer associated with the node. */
void *
CxTrNodeAuxGet(CxtTr *a_tr, CxtTrNode a_node);

/* Set the value of the auxiliary pointer associated with the node. */
void
CxTrNodeAuxSet(CxtTr *a_tr, CxtTrNode a_node, void *a_aux);

/* Get the degree (number of edges) of the node. */
uint32_t
CxTrNodeDegree(CxtTr *a_tr, CxtTrNode a_node);

/* Get the number of edges between a_node and a_other.  0 means that there is
 * no path between the two nodes. */
uint32_t
CxTrNodeDistance(CxtTr *a_tr, CxtTrNode a_node, CxtTrNode a_other);

/******************************************************************************/

/* tr_edge. */

/* Constructor. */
CxtTrEdge
CxTrEdgeNew(CxtTr *a_tr);

/* Destructor. */
void
CxTrEdgeDelete(CxtTr *a_tr, CxtTrEdge a_edge);

/* Get node a_i of the edge (a_i must be 0 or 1). */
CxtTrNode
CxTrEdgeNodeGet(CxtTr *a_tr, CxtTrEdge a_edge, uint32_t a_i);

/* Get the next edge in the edge ring at the a_i end of a_edge. */
void
CxTrEdgeNextGet(CxtTr *a_tr, CxtTrEdge a_edge, uint32_t a_i,
		 CxtTrEdge *r_next, uint32_t *r_end);

/* Get the previous edge in the edge ring at the a_i end of a_edge. */
void
CxTrEdgePrevGet(CxtTr *a_tr, CxtTrEdge a_edge, uint32_t a_i,
		 CxtTrEdge *r_prev, uint32_t *r_end);

/* Get the edge length. */
double
CxTrEdgeLengthGet(CxtTr *a_tr, CxtTrEdge a_edge);

/* Set the edge length. */
void
CxTrEdgeLengthSet(CxtTr *a_tr, CxtTrEdge a_edge, double a_length);

/* Get the value of the auxiliary pointer associated with the edge. */
void *
CxTrEdgeAuxGet(CxtTr *a_tr, CxtTrEdge a_edge);

/* Set the value of the auxiliary pointer associated with the edge. */
void
CxTrEdgeAuxSet(CxtTr *a_tr, CxtTrEdge a_edge, void *a_aux);

/* Attach edge to two nodes. */
void
CxTrEdgeAttach(CxtTr *a_tr, CxtTrEdge a_edge, CxtTrNode a_node_a,
	       CxtTrNode a_node_b);

/* Detach edge from nodes. */
void
CxTrEdgeDetach(CxtTr *a_tr, CxtTrEdge a_edge);
