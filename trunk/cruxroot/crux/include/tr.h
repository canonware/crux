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

typedef struct cw_trn_s cw_trn_t;
typedef struct cw_trt_s cw_trt_t;
typedef struct cw_tr_s cw_tr_t;

/* Canonical unrooted bifurcating phylogenetic tree in compact form. */
typedef cw_uint8_t cw_trs_t;

/* Tree node for an unrooted bifurcating phylogenetic tree. */
struct cw_trn_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TRN_MAGIC 0x63329478
    /* If this trn has been inserted into a tree, it should no longer be
     * externally modified.  This variable aids in asserting this invariant. */
    cw_bool_t in_tr;
#endif
    /* If non-NULL, then the node was dynamically allocated. */
    cw_mema_t *mema;

    /* Auxiliary opaque data pointer.  This is used by the treenode wrapper code
     * for reference iteration. */
    void *aux;

    /* If 0xffffffff, then the node is not a leaf node. */
#define CW_TRN_TAXON_NONE 0xffffffffU
    cw_uint32_t taxon_num;

    /* Pointers to neighbors.  Only the first two elements are used if the node
     * is a leaf node. */
#define CW_TRN_MAX_NEIGHBORS 3
    cw_trn_t *neighbors[CW_TRN_MAX_NEIGHBORS];
#define CW_TRN_EDGE_NONE 0xffffffffU
};

/* TBR neighbor. */
struct cw_trt_s
{
    /* Number of neighbors that can be reached by doing TBR at edges before this
     * one.  This is also the neighbor number of the first neighbor that can be
     * reached by doing TBR on this edge. */
    cw_uint32_t offset;

    /* Number of edges in the two subtrees.  Note that 0 and 1 are different
     * logical cases, but the number of connections possible for those two cases
     * is the same. */
    cw_uint32_t nedges_a;
    cw_uint32_t nedges_b;
};

struct cw_tr_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TR_MAGIC 0x39886394
#endif
    /* If non-NULL, then the node was dynamically allocated. */
    cw_mema_t *mema;

    /* Auxiliary opaque data pointer.  This is used by the treenode wrapper code
     * for reference iteration. */
    void *aux;

    /* Number of taxa in tree. */
    cw_uint32_t ntaxa;

    /* TRUE if this is a rooted tree, false otherwise. */
    cw_bool_t rooted;

    /* Root node.  If this tree is unrooted, the root node has no neighbors. */
    cw_trn_t *root;

    /* Taxon 0. */
    cw_trn_t *croot;

    /* Undo information for the most recent TBR. */
    cw_bool_t tbr_undoable;
    cw_uint32_t tbr_undo_bisect;
    cw_uint32_t tbr_undo_edge_a;
    cw_uint32_t tbr_undo_edge_b;

    /* Array of triplets that store information that is used for TBR-related
     * functions. */
    cw_trt_t *trt;
};

/* trn. */

/* Constructor. */
cw_trn_t *
trn_new(cw_trn_t *a_trn, cw_mema_t *a_mema);

/* Destructor. */
void
trn_delete(cw_trn_t *a_trn);

/* Get the taxon number associated with a_trn.  Return CW_TRN_TAXON_NONE if no
 * taxon number is set. */
cw_uint32_t
trn_taxon_num_get(cw_trn_t *a_trn);

/* Set the taxon number associated with a_trn (use CW_TRN_TAXON_NONE to unset
 * the taxon number. */
void
trn_taxon_num_set(cw_trn_t *a_trn, cw_uint32_t a_taxon_num);

/* Get neighbor a_i of the trn. */
cw_trn_t *
trn_neighbor_get(cw_trn_t *a_trn, cw_uint32_t a_i);

/* Swap two neighbors of a trn. */
void
trn_neighbors_swap(cw_trn_t *a_trn, cw_uint32_t a_i, cw_uint32_t a_j);

/* Join two trn's. */
void
trn_join(cw_trn_t *a_a, cw_trn_t *a_b);

/* Detatch two trn's. */
void
trn_detach(cw_trn_t *a_a, cw_trn_t *a_b);

/* Get the value of the auxiliary pointer associated with the trn. */
void *
trn_aux_get(cw_trn_t *a_trn);

/* Set the value of the auxiliary pointer associated with the trn. */
void
trn_aux_set(cw_trn_t *a_trn, void *a_aux);

/* tr. */

/* Constructor. */
cw_tr_t *
tr_new(cw_tr_t *a_tr, cw_mema_t *a_mema, cw_trn_t *a_trn);

/* Destructor. */
void
tr_delete(cw_tr_t *a_tr);

/* Get the number of taxa in the tree. */
cw_uint32_t
tr_ntaxa_get(cw_tr_t *a_tr);

/* Get the number of edges in the tree. */
cw_uint32_t
tr_nedges_get(cw_tr_t *a_tr);

/* Get edge a_edge. */
void
tr_edge_get(cw_tr_t *a_tr, cw_uint32_t a_edge, cw_trn_t **r_trn,
	    cw_uint32_t *r_neighbor);

/* Get the edge index of the edge between two trn's. */
cw_uint32_t
tr_edge_index_get(cw_tr_t *a_tr, cw_trn_t *a_trn_a, cw_trn_t *a_trn_b);

/* Get the root of the tree, were this a canonical tree (may or may not be).
 * This is always the lowest numbered taxon node. */
cw_trn_t *
tr_croot_get(cw_tr_t *a_tr);

/* Get the root of the tree.  If the tree is rooted, this node will not have
 * any neighbors. */
cw_trn_t *
tr_root_get(cw_tr_t *a_tr);

/* Return whether the tree is rooted. */
cw_bool_t
tr_rooted(cw_tr_t *a_tr);

/* Root the tree at the edge between a_trn and neighbor a_neighbor. */
void
tr_root(cw_tr_t *a_tr, cw_trn_t *a_trn, cw_uint32_t a_neighbor);

/* Unroot the tree. */
void
tr_unroot(cw_tr_t *a_tr);

/* Canonize the tree.  The tree must be unrooted. */
void
tr_canonize(cw_tr_t *a_tr);

/* Perform TBR. */
void
tr_tbr(cw_tr_t *a_tr, cw_uint32_t a_bisect, cw_uint32_t a_reconnect_a,
       cw_uint32_t a_reconnect_b, cw_uint32_t *r_bisect,
       cw_uint32_t *r_reconnect_a, cw_uint32_t *r_reconnect_b);

/* Undo the previous TBR.  Using this rather than manually reversing the
 * previous TBR allows internal caching of neighbors to still be used. */
void
tr_tbr_undo(cw_tr_t *a_tr);

/* Get the number of neighboring trees reachable via TBR. */
cw_uint32_t
tr_tbr_nneighbors_get(cw_tr_t *a_tr);

/* Get the parameters necessary to transorm this tree to neighbor a_neighbor. */
void
tr_tbr_neighbor_get(cw_tr_t *a_tr, cw_uint32_t a_neighbor,
		    cw_uint32_t *r_bisect, cw_uint32_t *r_reconnect_a,
		    cw_uint32_t *r_reconnect_b);

/* trs. */

/* Constructor. */
cw_trs_t *
trs_new(cw_mema_t *a_mema, cw_tr_t *a_tr);

/* Destructor. */
void
trs_delete(cw_trs_t *a_trs, cw_mema_t *a_mema, cw_uint32_t a_ntaxa);

/* Get the number of taxa in the tree. */
cw_uint32_t
trs_ntaxa(const cw_trs_t *a_trs);

/* Get the size (in bytes) of a trs with a_ntaxa taxa. */
cw_uint32_t
trs_ntaxa2sizeof(cw_uint32_t a_ntaxa);

/* Copy a_trs to a string. */
void
trs_memcopy(cw_uint8_t *a_dest, cw_trs_t *a_trs);

/* Hash a trs. */
cw_uint32_t
trs_hash(const void *a_key);

/* Compare two trs's. */
cw_bool_t
trs_key_comp(const void *a_k1, const void *a_k2);
