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
typedef struct cw_trn_s cw_trn_t;

/* Canonical unrooted bifurcating phylogenetic tree in compact form. */
typedef cw_uint8_t * cw_trc_t;

/* Tree node, used by cw_tr_t. */
struct cw_trn_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TRN_MAGIC 0x63329478
#endif
    cw_bool_t is_malloced;

    /* If this is NULL, then the node is not a leaf node. */
    cw_tx_t *taxon;

    /* Pointers to neighbors.  Only the first two elements are used if the node
     * is a leaf node. */
    cw_trn_t *neighbors[3];
};

/* Unrooted bifurcating phylogenetic tree. */
struct cw_tr_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TR_MAGIC 0x37d478a3
#endif

    /* Although this is an unrooted tree, we always start traversal from the
     * node connected to the first taxon, which effectively makes this the
     * root. */
    cw_trn_t *start;
};

/* trc. */
cw_trc_t *
trc_new(cw_trc_t *a_trc, cw_tr_t *a_tr);

void
trc_delete(cw_trc_t *a_trc);

void
trc_tr(cw_trc_t *a_trc, cw_tr_t *a_tr);

/* trn. */
cw_trn_t *
trn_new(cw_trn_t *a_trn);

void
trn_delete(cw_trn_t *a_trn);

void
trn_leaf(cw_trn_t *a_trn);

void
trn_internal(cw_trn_t *a_trn);

cw_tx_t *
trn_tx_get(cw_trn_t *a_trn);

void
trn_tx_set(cw_trn_t *a_trn, cw_tx_t *a_tx);

cw_trn_t *
trn_neighbor_get(cw_trn_t *a_trn, cw_uint32_t a_neighbor);

void
trn_neighbors_swap(cw_trn_t *a_trn, cw_uint32_t a_i, cw_uint32_t a_j);

void
trn_join(cw_trn_t *a_a, cw_uint32_t a_a_i, cw_trn_t *a_b, cw_uint32_t a_b_i);

void
trn_detach(cw_trn_t *a_a, cw_trn_t *a_b);

/* tr. */
cw_tr_t *
tr_new(cw_tr_t *a_tr);

void
tr_delete(cw_tr_t *a_tr);

cw_bool_t
tr_is_complete(cw_tr_t *a_tr, cw_txl_t *a_txl);

cw_trc_t *
tr_trc(cw_tr_t *a_tr, cw_trc_t *a_trc);

cw_uint32_t
tr_parsimony(cw_tr_t *a_tr);
