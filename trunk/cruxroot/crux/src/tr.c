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
 ******************************************************************************
 *
 * cw_tr_t contains a space-efficient representation of an unrooted bifurcating
 * phylogenetic tree.  The feasible duration of each iteration of tabu search is
 * bounded by the number of trees (and assoicated tabu lists) that can be
 * stored.  As such, the internal tree representation is as compact as
 * reasonably possible.
 *
 * Since tabu search depends on being able to determine if a tree has been
 * previously visited, tree comparison must be fast.  Therefore, a canonical
 * form is defined for trees.  Since the internal representation is effectively
 * string, it can be used as the input to a string hashing function, which makes
 * searching for a specific tree among previously visited trees a constant time
 * operation.
 *
 * A tree is stored in parenthetical form, generated by doing an in-order
 * traversal, starting at the first taxon, and ordering subtrees (left vs right)
 * according to which subtree has the lowest-numbered taxon.  The leaf node to
 * taxon mapping is stored as an ordered list that corresponds to the order in
 * which leaves are visited during the aforementioned in-order traversal.  The
 * first taxon (effectively the root) is left out of the list, since its
 * position is always the same.
 *
 * The parentheses are represented as a string of bits, where 0 is `(', and 1 is
 * `)'.
 *
 * The taxon list can be thought of as a series of permutations, where element i
 *             / n - i \
 * represents: |       |.  As the number of choices decreases, so does the
 *             \   1   /
 * number of bits that are necessary to distinguish among the choices.
 *
 * The following example is in canonical form:
 *
 *            A           C
 *             \         /
 *              \       /
 *               \     /
 *                \   /
 *                 \ /
 *                  V
 *                  |
 *                  |
 *                  |
 *                  |
 *                  |
 *                 / \
 *                /   \
 *               /     \
 *              /       \
 *   F---------/         \---------B
 *             |         |
 *             |         |
 *             |         |
 *             |         |
 *             |         |
 *             D         E
 *
 * Doing an in-order traversal of the tree results in the following expression:
 *
 *   (((BE)(DF))C)
 *
 * The topology for this tree can be represented by the following bit string:
 *
 *   ((()()))
 *   00010111
 *
 * The taxon visitation order can be represented by the following bit string
 * (keep in mind that only as many bits as are necessary are used for each
 * element):
 *
 *   A  B  E D F C
 *   - 00 10 1 1 -
 *
 * The `-' characters represent implicit information (A always comes first, and
 * 1 choose 1 is always the same).
 *
 * In general, the number of bytes that are needed to store a tree can be
 * calculated via the following formula:
 *
 *   __                          __
 *   |           n-2              |
 *   |          ____  __       __ |
 *   |          \     |         | |
 *   | 2(n-3) +  >    |log (n-i)| |
 *   |          /___  |   2     | |
 *   |           i=0              |
 *   | -------------------------- |
 *   |              8             |
 *   |                            |
 *
 * The first term in the numerator corresponds to the parenthetical expression,
 * and the summation term corresponds to the taxon permutation.
 *
 * In summary, a tree is represented by a parenthetical expression, immediately
 * followed by the taxon visitation order permutation.  The number of taxa in
 * the tree is implied by the parenthetical expression, and a well known taxon
 * ordering is assumed.
 *
 ******************************************************************************
 *
 * cw_trn_t contains a time-efficient representation of a node of an unrooted
 * bifurcating phlyogenetic tree.  This data structure takes much more space
 * than cw_tr_t does, but tree operations are much more efficient than they
 * would be for cw_tr_t.
 *
 ******************************************************************************/

#include "../include/modcrux.h"

/* tr. */
static cw_uint32_t
tr_p_len(const cw_tr_t *a_tr)
{
    cw_uint32_t retval;

    return retval;
}

cw_tr_t *
tr_new(cw_tr_t *a_tr, cw_mema_t *a_mema, cw_trn_t *a_trn, cw_uint32_t a_ntaxa)
{
    return NULL; /* XXX */
}

void
tr_delete(cw_tr_t *a_tr, cw_mema_t *a_mema, cw_uint32_t a_ntaxa)
{
    /* XXX */
}

cw_trn_t *
tr_trn(cw_tr_t *a_tr, cw_mema_t *a_mema, cw_uint32_t a_ntaxa)
{
    return NULL; /* XXX */
}

cw_uint32_t
tr_hash(const void *a_key)
{
    cw_uint32_t retval;
    const cw_tr_t *tr;
    cw_uint32_t len, i;

    cw_check_ptr(a_key);

    tr = (const cw_tr_t *) a_key;;

    len = tr_p_len(tr);
    for (i = retval = 0; i < len; i++)
    {
	retval = retval * 33 + ((cw_uint8_t *) tr)[i];
    }

    return retval;
}

cw_bool_t
tr_key_comp(const void *a_k1, const void *a_k2)
{
    cw_bool_t retval;
    const cw_tr_t *tr1, *tr2;
    cw_uint32_t len1, len2;

    cw_check_ptr(a_k1);
    cw_check_ptr(a_k2);

    tr1 = (const cw_tr_t *) a_k1;
    tr2 = (const cw_tr_t *) a_k2;

    len1 = tr_p_len(tr1);
    len2 = tr_p_len(tr2);
    if (len1 != len2)
    {
	retval = FALSE;
    }
    else
    {
	retval = memcmp(tr1, tr2, len1) ? FALSE : TRUE;
    }

    return retval;
}

/* trn. */
cw_trn_t *
trn_new(cw_trn_t *a_trn, cw_mema_t *a_mema)
{
    return NULL; /* XXX */
}

void
trn_delete(cw_trn_t *a_trn)
{
    /* XXX */
}

cw_uint32_t
trn_taxon_num_get(cw_trn_t *a_trn)
{
    cw_check_ptr(a_trn);
    cw_dassert(a_trn->magic == CW_TRN_MAGIC);

    return a_trn->taxon_num;
}

void
trn_taxon_num_set(cw_trn_t *a_trn, cw_uint32_t a_taxon_num)
{
    cw_check_ptr(a_trn);
    cw_dassert(a_trn->magic == CW_TRN_MAGIC);

    a_trn->taxon_num = a_taxon_num;
}

cw_trn_t *
trn_neighbor_get(cw_trn_t *a_trn, cw_uint32_t a_i)
{
    return NULL; /* XXX */
}

void
trn_neighbors_swap(cw_trn_t *a_trn, cw_uint32_t a_i, cw_uint32_t a_j)
{
    /* XXX */
}

void
trn_join(cw_trn_t *a_a, cw_uint32_t a_a_i, cw_trn_t *a_b, cw_uint32_t a_b_i)
{
    /* XXX */
}

void
trn_detach(cw_trn_t *a_a, cw_trn_t *a_b)
{
    /* XXX */
}

cw_tr_t *
trn_tr(cw_trn_t *a_tr, cw_mema_t *a_mema, cw_uint32_t a_ntaxa)
{
    return NULL; /* XXX */
}
