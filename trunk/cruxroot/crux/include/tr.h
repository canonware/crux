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
#define CW_TR_NODE_MAX_NEIGHBORS 3
#define CW_TR_NODE_EDGE_NONE 0xffffffffU

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

/* Get edge a_edge. */
void
tr_edge_get(cw_tr_t *a_tr, uint32_t a_edge, cw_tr_node_t *r_node_a,
	    cw_tr_node_t *r_node_b);

/* Get the edge index of the edge between two nodes. */
uint32_t
tr_edge_index_get(cw_tr_t *a_tr, cw_tr_node_t a_node_a, cw_tr_node_t a_node_b);

/* Get the base of the tree. */
cw_tr_node_t
tr_base_get(cw_tr_t *a_tr);

/* Set the base of the tree. */
void
tr_base_set(cw_tr_t *a_tr, cw_tr_node_t a_base);

/* Canonize the tree. */
void
tr_canonize(cw_tr_t *a_tr);

/* Perform TBR. */
void
tr_tbr(cw_tr_t *a_tr, uint32_t a_bisect, uint32_t a_reconnect_a,
       uint32_t a_reconnect_b);

/* Get the number of neighboring trees reachable via TBR. */
uint32_t
tr_tbr_nneighbors_get(cw_tr_t *a_tr);

/* Get the parameters necessary to transorm this tree to neighbor a_neighbor. */
void
tr_tbr_neighbor_get(cw_tr_t *a_tr, uint32_t a_neighbor,
		    uint32_t *r_bisect, uint32_t *r_reconnect_a,
		    uint32_t *r_reconnect_b);

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
tr_mp_prepare(cw_tr_t *a_tr, uint8_t *a_taxa[], uint32_t a_ntaxa,
	      uint32_t a_nchars);

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

/* Get neighbor a_i of the node. */
cw_tr_node_t
tr_node_neighbor_get(cw_tr_t *a_tr, cw_tr_node_t a_node, uint32_t a_i);

/* Swap two neighbors of a node. */
void
tr_node_neighbors_swap(cw_tr_t *a_tr, cw_tr_node_t a_node, uint32_t a_i,
		       uint32_t a_j);

/* Join two nodes. */
void
tr_node_join(cw_tr_t *a_tr, cw_tr_node_t a_a, cw_tr_node_t a_b);

/* Detatch two nodes. */
void
tr_node_detach(cw_tr_t *a_tr, cw_tr_node_t a_a, cw_tr_node_t a_b);

/* Get the value of the auxiliary pointer associated with the node. */
void *
tr_node_aux_get(cw_tr_t *a_tr, cw_tr_node_t a_node);

/* Set the value of the auxiliary pointer associated with the node. */
void
tr_node_aux_set(cw_tr_t *a_tr, cw_tr_node_t a_node, void *a_aux);
