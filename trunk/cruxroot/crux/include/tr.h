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

typedef cw_uint32_t cw_tr_node_t;
#define CW_TR_NODE_NONE 0xffffffffU
#define CW_TR_NODE_TAXON_NONE 0xffffffffU
#define CW_TR_NODE_MAX_NEIGHBORS 3
#define CW_TR_NODE_EDGE_NONE 0xffffffffU

/******************************************************************************/

/* tr. */

/* Constructor. */
cw_tr_t *
tr_new(void);

/* Destructor. */
void
tr_delete(cw_tr_t *a_tr);

/* Iterate over nodes.  Pass CW_TR_NODE_NONE the first time, and stop when
 * CW_TR_NODE_NONE is returned. */
cw_tr_node_t
tr_nodes_iterate(cw_tr_t *a_tr, cw_tr_node_t a_prev);

/* Get the number of taxa in the tree. */
cw_uint32_t
tr_ntaxa_get(cw_tr_t *a_tr);

/* Get the number of edges in the tree. */
cw_uint32_t
tr_nedges_get(cw_tr_t *a_tr);

/* Get edge a_edge. */
void
tr_edge_get(cw_tr_t *a_tr, cw_uint32_t a_edge, cw_tr_node_t *r_node,
	    cw_uint32_t *r_neighbor);

// XXX Remove?
/* Get the edge index of the edge between two nodes. */
cw_uint32_t
tr_edge_index_get(cw_tr_t *a_tr, cw_tr_node_t a_node_a, cw_tr_node_t a_node_b);

/* Get the root of the tree, were this a canonical tree (may or may not be).
 * This is always the lowest numbered taxon node. */
cw_tr_node_t
tr_croot_get(cw_tr_t *a_tr);

/* Canonize the tree. */
void
tr_canonize(cw_tr_t *a_tr);

/* Perform TBR. */
void
tr_tbr(cw_tr_t *a_tr, cw_uint32_t a_bisect, cw_uint32_t a_reconnect_a,
       cw_uint32_t a_reconnect_b, cw_uint32_t *r_bisect,
       cw_uint32_t *r_reconnect_a, cw_uint32_t *r_reconnect_b);

/* Get the number of neighboring trees reachable via TBR. */
cw_uint32_t
tr_tbr_nneighbors_get(cw_tr_t *a_tr);

/* Get the parameters necessary to transorm this tree to neighbor a_neighbor. */
void
tr_tbr_neighbor_get(cw_tr_t *a_tr, cw_uint32_t a_neighbor,
		    cw_uint32_t *r_bisect, cw_uint32_t *r_reconnect_a,
		    cw_uint32_t *r_reconnect_b);

/* Get the value of the auxiliary pointer associated with the tr. */
void *
tr_aux_get(cw_tr_t *a_tr);

/* Set the value of the auxiliary pointer associated with the tr. */
void
tr_aux_set(cw_tr_t *a_tr, void *a_aux);

/* Prepare for calculating Fitch parsimony.  a_taxa points to an array of
 * character array pointers, where the index into a_taxa corresponds to taxon
 * number.  The character arrays need not be nil-terminated. */
void
tr_mp_prepare(cw_tr_t *a_tr, cw_uint8_t *a_taxa[], cw_uint32_t a_ntaxa,
	      cw_uint32_t a_nchars);

/* Calculate the Fitch parsimony score for this tree.  If a_maxscore is not
 * CW_MAXSCORE_NONE, terminate scoring if a_maxscore is reached/exceeded and
 * return CW_TR_MAXSCORE_NONE. */
cw_uint32_t
tr_mp_score(cw_tr_t *a_tr, cw_uint32_t a_maxscore);

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
cw_uint32_t
tr_node_taxon_num_get(cw_tr_t *a_tr, cw_tr_node_t a_node);

/* Set the taxon number associated with a_node (use CW_TR_NODE_TAXON_NONE to
 * unset the taxon number. */
void
tr_node_taxon_num_set(cw_tr_t *a_tr, cw_tr_node_t a_node,
		      cw_uint32_t a_taxon_num);

/* Get neighbor a_i of the node. */
cw_tr_node_t
tr_node_neighbor_get(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_uint32_t a_i);

/* Swap two neighbors of a node. */
void
tr_node_neighbors_swap(cw_tr_t *a_tr, cw_tr_node_t a_node, cw_uint32_t a_i,
		       cw_uint32_t a_j);

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
