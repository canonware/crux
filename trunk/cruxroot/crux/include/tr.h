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

/* Canonical unrooted bifurcating phylogenetic tree in compact form. */
typedef cw_uint8_t cw_tr_t;

/* Tree node for an unrooted bifurcating phylogenetic tree. */
struct cw_trn_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TRN_MAGIC 0x63329478
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
};

/* tr. */
cw_tr_t *
tr_new(cw_mema_t *a_mema, cw_trn_t *a_trn, cw_uint32_t a_ntaxa);

void
tr_delete(cw_tr_t *a_tr, cw_mema_t *a_mema, cw_uint32_t a_ntaxa);

cw_uint32_t
tr_ntaxa(const cw_tr_t *a_tr);

cw_uint32_t
tr_ntaxa2sizeof(cw_uint32_t a_ntaxa);

void
tr_memcpy(cw_uint8_t *a_dest, cw_tr_t *a_tr);

cw_trn_t *
tr_trn(cw_tr_t *a_tr, cw_mema_t *a_mema, cw_uint32_t a_ntaxa);

cw_uint32_t
tr_hash(const void *a_key);

cw_bool_t
tr_key_comp(const void *a_k1, const void *a_k2);

/* trn. */
cw_trn_t *
trn_new(cw_trn_t *a_trn, cw_mema_t *a_mema);

void
trn_delete(cw_trn_t *a_trn);

cw_uint32_t
trn_taxon_num_get(cw_trn_t *a_trn);

void
trn_taxon_num_set(cw_trn_t *a_trn, cw_uint32_t a_taxon_num);

cw_trn_t *
trn_neighbor_get(cw_trn_t *a_trn, cw_uint32_t a_i);

void
trn_neighbors_swap(cw_trn_t *a_trn, cw_uint32_t a_i, cw_uint32_t a_j);

void
trn_join(cw_trn_t *a_a, cw_trn_t *a_b);

void
trn_detach(cw_trn_t *a_a, cw_trn_t *a_b);

void *
trn_aux_get(cw_trn_t *a_trn);

void
trn_aux_set(cw_trn_t *a_trn, void *a_aux);

/* Functions that operate on trn-based trees. */
cw_trn_t *
trn_tree_random(cw_mema_t *a_mema, cw_uint32_t a_ntaxa);

void
trn_tree_delete(cw_trn_t *a_trn);

cw_uint32_t
trn_tree_ntaxa_get(cw_trn_t *a_trn);

cw_trn_t *
trn_tree_root_get(cw_trn_t *a_trn);
