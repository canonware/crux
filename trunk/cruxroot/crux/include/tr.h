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

/* Tree node for an unrooted bifurcating phylogenetic tree. */
struct cw_trn_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TRN_MAGIC 0x63329478
#endif

    /* Auxiliary opaque data pointer.  This is used by the treenode wrapper code
     * for reference iteration. */
    void *aux;

    /* If 0xffffffff, then the node is not a leaf node. */
#define CW_TRN_TAXON_NONE 0xffffffffU
    cw_uint32_t taxon_num;

    /* Pointers to neighbors.  Only the first element is used if the node is a
     * leaf node. */
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

    /* Bisection edge. */
    cw_uint32_t bisect_edge;

    /* Number of edges in the two subtrees.  Note that 0 and 1 are different
     * logical cases, but the number of connections possible for those two cases
     * is the same. */
    cw_uint32_t nedges_a;
    cw_uint32_t nedges_b;

    /* Edge indices for the edges that will reverse bisection.  This is used to
     * avoid enumerating reconnections that undo the bisections. */
    cw_uint32_t self_a;
    cw_uint32_t self_b;
};

struct cw_tr_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TR_MAGIC 0x39886394
#endif

    /* Auxiliary opaque data pointer.  This is used by the treenode wrapper code
     * for reference iteration. */
    void *aux;

    /* TRUE if this is a rooted tree, false otherwise. */
    cw_bool_t rooted;

    /* Root node.  If this tree is unrooted, the root node has no neighbors. */
    cw_trn_t *root;

    /* Taxon 0 (canonical tree root). */
    cw_trn_t *croot;

    /* Number of taxa in tree. */
    cw_uint32_t ntaxa;

    /* Number of edges in tree, assuming that the tree is unrooted.  This can be
     * derived from ntaxa, but is used often enough to make storing it
     * worthwhile. */
    cw_uint32_t nedges;

    /* Undo information for the most recent TBR. */
    cw_bool_t tbr_undoable;
    cw_uint32_t tbr_undo_bisect;
    cw_uint32_t tbr_undo_reconnect_a;
    cw_uint32_t tbr_undo_reconnect_b;
#ifdef CW_DBG
    cw_uint8_t *tbr_undo_canonical;
    cw_uint8_t *tbr_undone_canonical;
#endif

    /* Array of triplets that store per-edge information that is used for
     * TBR-related functions.  There is one more element in trt than there are
     * edges in the tree.  This is critical to the way binary searching on the
     * array is done, and it also makes it easy to get the total number of
     * TBR neighbors this tree has (trt[nedges].offset).
     *
     * Only the first trtlen elements are valid, since not all bisection edges
     * necessarily result in neighbors. */
    cw_trt_t *trt;
    cw_uint32_t trtlen;

    /* Array of integers that is used for randomly iterating over all TBR
     * neighbors (which are enumerated by trt).  The following algorithm is used
     * for the random iteration:
     *
     *   1) memset(trr, 0xff, trt[nedges].offset * sizeof(cw_uint32_t)).
     *
     *   2) Initialize trri <-- 0.
     * 
     *   3) Choose a random slot r in [trri,trt[nedges].offset).
     *
     *   4) If trr[r] == 0xffffffff, set trr[r] <-- r.
     *
     *   5) If trr[trri] == 0xffffffff, set trr[trri] <-- trri.
     *
     *   6) Swap trr[trri] <--> trr[r].  trr[trri] contains the choice.
     *
     *   7) Increment trri.
     *
     *   8) If trri < trt[nedges].offset, go to step 3.
     *
     * trr can be re-initialized by iterating over the first trri elements and
     * clearing the elements associated with the values stored, then clearing
     * the first trri elements. */
    cw_uint32_t *trr;
    /* Number of elements in trr.  This is not necessarily the number of
     * elements that are actually being used. */
    cw_uint32_t trrlen;
    /* Number of integers that have been chosen from trr. */
    cw_uint32_t trri;
};

/* trn. */

/* Constructor. */
void
trn_new(cw_trn_t *a_trn);

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
void
tr_new(cw_tr_t *a_tr, cw_trn_t *a_root);

/* Destructor. */
void
tr_delete(cw_tr_t *a_tr, cw_bool_t a_delete_trns);

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

/* Root the tree at the edge between a_trn_a and a_trn_b. */
void
tr_root(cw_tr_t *a_tr, cw_trn_t *a_trn_a, cw_trn_t *a_trn_b);

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

/* Get the number of neighbors that have been randomly chosen via
 * tr_tbr_rneighbor_get(). */
cw_uint32_t
tr_tbr_rneighbor_nchosen_get(cw_tr_t *a_tr);

/* Get the next randomly chosen neighbor index, which can in turn be used as the
 * a_neighbor argument to tr_tbr_neighbor_get(). */
cw_uint32_t
tr_tbr_rneighbor_get(cw_tr_t *a_tr, cw_mt_t *a_mt);

/* Create a canonical string representation of a_tr (which must be unrooted).
 * ar_string must point to a_len bytes of storage, which must in turn be
 * tr_string_ntaxa2sizeof(tr_ntaxa_get(a_tr)). */
void
tr_string(cw_tr_t *a_tr, cw_uint8_t *ar_string, cw_uint32_t a_len);

/* Get the number of taxa in the canonical tree represented by a_string. */
cw_uint32_t
tr_string_ntaxa(const cw_uint8_t *a_string);

/* Return the number of bytes needed to store a canonical tree with a_ntaxa taxa
 * in string format. */
cw_uint32_t
tr_string_ntaxa2sizeof(cw_uint32_t a_ntaxa);

/* Hash a string that represents a canonical tree. */
cw_uint32_t
tr_string_hash(const void *a_key);

/* Compare two strings that represent canonical trees. */
cw_bool_t
tr_string_key_comp(const void *a_k1, const void *a_k2);

/* Get the value of the auxiliary pointer associated with the tr. */
void *
tr_aux_get(cw_tr_t *a_tr);

/* Set the value of the auxiliary pointer associated with the tr. */
void
tr_aux_set(cw_tr_t *a_tr, void *a_aux);
