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

typedef struct cw_tr_s cw_tr_t;
#define CW_TR_MAXSCORE_NONE 0xffffffffU
#define CW_TR_HOLD_ALL 0xffffffffU

typedef uint32_t cw_tr_node_t;
#define CW_TR_NODE_NONE 0xffffffffU
#define CW_TR_NODE_TAXON_NONE 0xffffffffU

typedef uint32_t cw_tr_edge_t;
#define CW_TR_EDGE_NONE 0xffffffffU

/******************************************************************************/

/* tr. */

/* Constructor. */
cw_tr_t *
tr_new(cw_mema_t *a_mema);

/* Destructor. */
void
tr_delete(cw_tr_t *a_tr);

/* Create a duplicate copy of a tree. */
cw_tr_t *
tr_dup(cw_tr_t *a_tr);

/* Get the number of taxa in the tree. */
uint32_t
tr_ntaxa_get(cw_tr_t *a_tr);

/* Get the number of edges in the tree. */
uint32_t
tr_nedges_get(cw_tr_t *a_tr);

/* Get the base of the tree. */
cw_tr_node_t
tr_base_get(cw_tr_t *a_tr);

/* Set the base of the tree. */
void
tr_base_set(cw_tr_t *a_tr, cw_tr_node_t a_base);

/* Canonize the tree.  The base of the tree is set to the lowest numbered taxon,
 * and internal nodes are adjusted such that their edge rings are ordered
 * (subtrees with lower minimum taxon numbers come first), and the edge returned
 * by tr_node_edge_get() is the edge that leads back to the base. */
void
tr_canonize(cw_tr_t *a_tr);

/* Perform TBR.  The callback functions are called when new edges or nodes are
 * needed; this never happens for strictly trifurcating trees. */
void
tr_tbr(cw_tr_t *a_tr, cw_tr_edge_t a_bisect, cw_tr_edge_t a_reconnect_a,
       cw_tr_edge_t a_reconnect_b,
       cw_tr_edge_t (*a_edge_alloc_callback)(cw_tr_t *, void *),
       cw_tr_node_t (*a_node_alloc_callback)(cw_tr_t *, void *),
       void *a_arg);

/* Get the number of neighboring trees reachable via TBR. */
uint32_t
tr_tbr_nneighbors_get(cw_tr_t *a_tr);

/* Get the parameters necessary to transorm this tree to neighbor a_neighbor. */
void
tr_tbr_neighbor_get(cw_tr_t *a_tr, uint32_t a_neighbor,
		    cw_tr_edge_t *r_bisect, cw_tr_edge_t *r_reconnect_a,
		    cw_tr_edge_t *r_reconnect_b);

/* Get the value of the auxiliary pointer associated with the tr. */
void *
tr_aux_get(cw_tr_t *a_tr);

/* Set the value of the auxiliary pointer associated with the tr. */
void
tr_aux_set(cw_tr_t *a_tr, void *a_aux);

/* Prepare for calculating Fitch parsimony scores.  a_taxa points to an array of
 * character array pointers, where the index into a_taxa corresponds to taxon
 * number.  The character arrays need not be nil-terminated. */
void
tr_mp_prepare(cw_tr_t *a_tr, bool a_uninformative_eliminate,
	      char *a_taxa[], uint32_t a_ntaxa, uint32_t a_nchars);

/* Clear the data structures used for calculating Fitch parsimony scores. */
void
tr_mp_finish(cw_tr_t *a_tr);

/* Calculate the Fitch parsimony score for this tree. */
uint32_t
tr_mp_score(cw_tr_t *a_tr);

/* Calculate the Fitch parsimony of all TBR neighbors, and keep track of up to
 * a_max_hold of the best neighbors (or all best neighbors, if a_max_hold is
 * CW_TR_HOLD_ALL). */
void
tr_tbr_best_neighbors_mp(cw_tr_t *a_tr, uint32_t a_max_hold);

/* Calculate the Fitch parsimony of all TBR neighbors, and keep track of up to
 * a_max_hold of the better neighbors (or all better neighbors, if a_max_hold is
 * CW_TR_HOLD_ALL). */
void
tr_tbr_better_neighbors_mp(cw_tr_t *a_tr, uint32_t a_max_hold);

/* Calculate the Fitch parsimony of all TBR neighbors, and keep track of all
 * neighbors. */
void
tr_tbr_all_neighbors_mp(cw_tr_t *a_tr);

/* Clear the data structures used to store held trees. */
void
tr_held_finish(cw_tr_t *a_tr);

/* Get the number of trees currently held. */
uint32_t
tr_nheld_get(cw_tr_t *a_tr);

/* Get the a_held'th held tree, and its score.  *r_neighbor can be passed to
 * tr_tbr_neighbor_get() in order to get the TBR transformation parameters,
 * which can then be passed to tr_tbr(). */
void
tr_held_get(cw_tr_t *a_tr, uint32_t a_held, uint32_t *r_neighbor,
	    uint32_t *r_score);

/******************************************************************************/

/* tr_node. */

/* Constructor. */
cw_tr_node_t
tr_node_new(cw_tr_t *a_tr);

/* Destructor. */
void
tr_node_delete(cw_tr_t *a_tr, cw_tr_node_t a_node);

/* Get the taxon number associated with a_node.  Return CW_TR_NODE_TAXON_NONE if
 * no taxon number is set. */
uint32_t
tr_node_taxon_num_get(cw_tr_t *a_tr, cw_tr_node_t a_node);

/* Set the taxon number associated with a_node (use CW_TR_NODE_TAXON_NONE to
 * unset the taxon number. */
void
tr_node_taxon_num_set(cw_tr_t *a_tr, cw_tr_node_t a_node,
		      uint32_t a_taxon_num);

/* Get the first edge in a_node's edge ring. */
void
tr_node_edge_get(cw_tr_t *a_tr, cw_tr_node_t a_node,
		 cw_tr_edge_t *r_edge, uint32_t *r_end);

/* Get the value of the auxiliary pointer associated with the node. */
void *
tr_node_aux_get(cw_tr_t *a_tr, cw_tr_node_t a_node);

/* Set the value of the auxiliary pointer associated with the node. */
void
tr_node_aux_set(cw_tr_t *a_tr, cw_tr_node_t a_node, void *a_aux);

/* Get the degree (number of edges) of the node. */
uint32_t
tr_node_degree(cw_tr_t *a_tr, cw_tr_node_t a_node);

/* Get the number of edges between a_node and a_other.  0 means that there is
 * no path between the two nodes. */
uint32_t
tr_node_distance(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_tr_node_t a_other);

/******************************************************************************/

/* tr_edge. */

/* Constructor. */
cw_tr_edge_t
tr_edge_new(cw_tr_t *a_tr);

/* Destructor. */
void
tr_edge_delete(cw_tr_t *a_tr, cw_tr_edge_t a_edge);

/* Get node a_i of the edge (a_i must be 0 or 1). */
cw_tr_node_t
tr_edge_node_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_i);

/* Get the next edge in the edge ring at the a_i end of a_edge. */
void
tr_edge_next_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_i,
		 cw_tr_edge_t *r_next, uint32_t *r_end);

/* Get the previous edge in the edge ring at the a_i end of a_edge. */
void
tr_edge_prev_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge, uint32_t a_i,
		 cw_tr_edge_t *r_prev, uint32_t *r_end);

/* Get the edge length. */
double
tr_edge_length_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge);

/* Set the edge length. */
void
tr_edge_length_set(cw_tr_t *a_tr, cw_tr_edge_t a_edge, double a_length);

/* Get the value of the auxiliary pointer associated with the edge. */
void *
tr_edge_aux_get(cw_tr_t *a_tr, cw_tr_edge_t a_edge);

/* Set the value of the auxiliary pointer associated with the edge. */
void
tr_edge_aux_set(cw_tr_t *a_tr, cw_tr_edge_t a_edge, void *a_aux);

/* Attach edge to two nodes. */
void
tr_edge_attach(cw_tr_t *a_tr, cw_tr_edge_t a_edge, cw_tr_node_t a_node_a,
	       cw_tr_node_t a_node_b);

/* Detach edge from nodes. */
void
tr_edge_detach(cw_tr_t *a_tr, cw_tr_edge_t a_edge);
