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
 * Consider the following example tree:
 *
 *
 *          E           B
 *           \         /
 *            \       /
 *             \     /
 *              \   /
 *               \ /
 *                V         G
 *                |         |
 *                |         |
 *                |         |
 *                |         |
 *                |         |
 *               / \        /---------D
 *              /   \      /
 *             /     \    /
 *            /       \  /
 *   A-------/         \/
 *           |         |
 *           |         |
 *           |         |
 *           |         |
 *           |         |
 *           C         F
 *
 * Assume the following ordering of taxa: A B C D E F G.  The above tree is
 * converted to canonical form by placing the lowest taxon (A) above the root
 * node, and re-ordering subtrees such that for every internal node, the left
 * subtree contains a lower taxon than any of the taxa in the right subtree.
 * Following is the canonical form of the above tree:
 *
 *
 *
 *
 *             A
 *             |
 *             |
 *             |
 *            /^\
 *           /   \
 *          /     \
 *         /\      C
 *        /  \
 *       /    \
 *      /\     \
 *     /  \     \
 *    /    \     \
 *   B      E    /\
 *              /  \
 *             /    \
 *            /\     F
 *           /  \
 *          /    \
 *         D      G
 *
 * Doing an in-order traversal of the tree results in the following expression:
 *
 *   A(((BE)((DG)F))C)
 *
 * The topology for this tree can be represented by the following bit string:
 *
 *   ((()(())))
 *   -001001111
 *
 * The `-' character represents implicit information.
 *
 * The taxon visitation order can be represented by the following bit string
 * (keep in mind that only as many bits as are necessary are used for each
 * element):
 *
 *   A B   E  D  G F C
 *   - - 010 01 10 1 -
 *
 * The `-' characters represent implicit information (A always comes first, B
 * always comes second, and 1 choose 1 is always the same).
 *
 * Finally, padding bits are added to round the total number of bit up to a
 * multiple of 8, in order to work well on byte-addressable machinery.
 *
 * 
 *   001001111010011010000000
 *
 * In general, the number of bytes that are needed to store an n-taxon tree can
 * be calculated via the following formula:
 *
 *   __                                __
 *   |                 n-3              |
 *   |                ____  __       __ |
 *   |                \     |         | |
 *   | (2(n-3) + 1) +  >    |log (n-i)| |
 *   |                /___  |   2     | |
 *   |                 i=0              |
 *   | -------------------------------- |
 *   |                    8             |
 *   |                                  |
 *
 * The first term in the numerator corresponds to the parenthetical tree
 * expression, and the summation term corresponds to the taxa permutation.
 *
 * In summary, a tree is represented by a parenthetical expression, immediately
 * followed by the taxon visitation order permutation, and padded to the nearest
 * byte boundary.  The number of taxa in the tree is implied by the
 * parenthetical expression, and a well known taxon ordering is assumed.
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

/* Prototypes. */
/* tr. */
static cw_uint32_t
tr_p_log2ceil(cw_uint32_t a_val);
static cw_uint32_t
tr_p_ntaxa2nbits(cw_uint32_t a_ntaxa);
static cw_uint32_t
tr_p_sizeof(const cw_tr_t *a_tr);
static int
tr_p_taxon_num_compar(const void *a_a, const void *a_b);
#ifdef CW_DBG
static cw_bool_t
tr_p_validate(const cw_tr_t *a_tr);
#endif
static void
tr_p_new_recurse(cw_tr_t *a_tr, cw_uint32_t *a_bitind_paren,
		 cw_uint32_t *a_bitind_perm, cw_trn_t *a_trn, cw_trn_t *a_prev,
		 cw_uint32_t *a_unchosen, cw_uint32_t *a_nunchosen);

/* trn. */
#ifdef CW_DBG
static cw_bool_t
trn_p_validate(cw_trn_t *a_trn);
static cw_uint32_t
trn_p_tree_validate_recurse(cw_trn_t *a_trn, cw_trn_t *a_prev,
			    cw_uint32_t a_taxon_num);
static cw_bool_t
trn_p_tree_validate(cw_trn_t *a_trn);
#endif
static void
trn_p_tree_delete(cw_trn_t *a_trn);
static cw_uint32_t
trn_p_tree_ntaxa_get(cw_trn_t *a_trn, cw_trn_t *a_prev);
static cw_trn_t *
trn_p_tree_root_get(cw_trn_t *a_trn, cw_trn_t *a_prev, cw_trn_t *a_root);
static cw_uint32_t
trn_p_tree_canonize(cw_trn_t *a_trn, cw_trn_t *a_prev);
static cw_bool_t
trn_p_tree_edge_get_recurse(cw_trn_t *a_trn, cw_uint32_t a_edge,
			    cw_trn_t *a_prev, cw_uint32_t *r_edge_count,
			    cw_trn_t **r_trn, cw_uint32_t *r_neighbor);
static void
trn_p_tree_edge_get(cw_trn_t *a_trn, cw_uint32_t a_edge, cw_trn_t **r_trn,
		    cw_uint32_t *r_neighbor);

/* Get bit a_i in a_vec (cw_uint8_t *). */
#define TR_BIT_GET(a_vec, a_i)						\
    (( ((a_vec)[(a_i) >> 3]) >> (7 - ((a_i) & 0x7)) ) & 0x1)

/* Set bit a_i in a_vec (cw_uint8_t *) to a_val (0 or 1). */
#define TR_BIT_SET(a_vec, a_i, a_val)					\
    (a_vec)[(a_i) >> 3] =						\
	((a_vec)[(a_i) >> 3] & (~(0x1 << (7 - ((a_i) & 0x7)))))		\
	| ((a_val) << ((7 - ((a_i) & 0x7))))

/* tr. */

/*           __         __
 *           |           |
 * Calculate |log (a_val)|.  This is accomplished by counting the number of bits
 *           |   2       |
 * set to 1 while iterating to find the most significant bit set.  If there
 * was more than one 1 bit, the result is incremented (ceiling).
 *
 * Note that the result for a_val == 1 should be 1, but this function would
 * return 0.  However, we never care about that case in this code, so it doesn't
 * matter.
 */
static cw_uint32_t
tr_p_log2ceil(cw_uint32_t a_val)
{
    cw_uint32_t retval;
    cw_uint32_t ones;

    cw_assert(a_val > 1);

    /* Find the most significant 1 bit. */
    for (retval = ones = 0; a_val != 0; retval++, a_val >>= 1)
    {
	ones += (a_val & 0x1);
    }

    /* Decrement unless ceiling needs to be taken. */
    if (ones == 1)
    {
	retval--;
    }

    return retval;
}

/* Calculatet bit vector bit count, given the number of taxa in the tree. */
static cw_uint32_t
tr_p_ntaxa2nbits(cw_uint32_t a_ntaxa)
{
    cw_uint32_t retval;
    cw_uint32_t i;

    /* Use the numerator of the formula at the beginning of this file to
     * calculate the number of bits required to store a tree with a_ntaxa. */

    /* First term of numerator (parenthetical tree). */
    retval = 2 * (a_ntaxa - 3) + 1;

    /* Second term of numerator (summation, taxa permutation). */
    for (i = 0; i < a_ntaxa - 2; i++)
    {
	retval += tr_p_log2ceil(a_ntaxa - i);
    }

    return retval;
}

/* Calculate bit vector byte count, using the bit representation of the encoded
 * tree. */
static cw_uint32_t
tr_p_sizeof(const cw_tr_t *a_tr)
{
    return tr_ntaxa2sizeof(tr_ntaxa(a_tr));
}

/* Comparison function passed to bsearch(3) when searching for a taxon to build
 * a taxa permutation. */
static int
tr_p_taxon_num_compar(const void *a_a, const void *a_b)
{
    int retval;
    cw_uint32_t *a = (cw_uint32_t *) a_a;
    cw_uint32_t *b = (cw_uint32_t *) a_b;

    if (*a < *b)
    {
	retval = -1;
    }
    else if (*a > *b)
    {
	retval = 1;
    }
    else
    {
	retval = 0;
    }

    return retval;
}

#ifdef CW_DBG
static cw_bool_t
tr_p_validate(const cw_tr_t *a_tr)
{
    cw_uint32_t ntaxa, nbytes, nbits;
    cw_uint32_t i, j, bitind, npbits, pbits;

    cw_check_ptr(a_tr);

    ntaxa = tr_ntaxa(a_tr);
    nbytes = tr_ntaxa2sizeof(ntaxa);
    nbits = tr_p_ntaxa2nbits(ntaxa);

    /* Assert that trailing bits are all 0. */
    cw_assert((a_tr[nbytes - 1] & (0xff >> (nbits & 0x7))) == 0);

    /* Assert that no invalid choices exist in the taxon permutation.  The taxa
     * permutation is stored in a format that makes it impossible to encode the
     * same taxon number twice.  However, there can still be invalid numbers in
     * the encoding (for example 7 cannot be used for a 5 choose 1 choice). */
    for (i = ntaxa - 2, bitind = 2 * (ntaxa - 3) + 1;
	 i > 1;
	 i--)
    {
	for (j = 0, npbits = tr_p_log2ceil(i), pbits = 0;
	     j < npbits;
	     j++, bitind++)
	{
	    pbits <<= 1;
	    pbits |= TR_BIT_GET(a_tr, bitind);

	    cw_assert(pbits < i);
	}
    }

    return TRUE;
}
#endif

static void
tr_p_new_recurse(cw_tr_t *a_tr, cw_uint32_t *a_bitind_paren,
		 cw_uint32_t *a_bitind_perm, cw_trn_t *a_trn, cw_trn_t *a_prev,
		 cw_uint32_t *a_unchosen, cw_uint32_t *a_nunchosen)
{
    cw_uint32_t i;

    if (a_trn->taxon_num != CW_TRN_TAXON_NONE)
    {
	if (*a_nunchosen > 1)
	{
	    cw_uint32_t *taxon, offset, nbits;

	    /* Leaf node. */

	    /* Avoid encoding taxon 1, since it is implicit. */
	    if (a_trn->taxon_num != 1)
	    {
		/* Get the offset of this taxon within the array of unchosen
		 * taxa. */
		taxon = (cw_uint32_t *) bsearch(&a_trn->taxon_num,
						a_unchosen, *a_nunchosen,
						sizeof(cw_uint32_t),
						tr_p_taxon_num_compar);
		cw_check_ptr(taxon);
		offset = (cw_uint32_t) (taxon - a_unchosen);

		/* Determine how many bits to use in storing this choice. */
		nbits = tr_p_log2ceil(*a_nunchosen);

		/* Remove the taxon from the array of unchosen taxa. */
		if (offset < *a_nunchosen - 1)
		{
		    memmove(&a_unchosen[offset], &a_unchosen[offset + 1],
			    (*a_nunchosen - offset - 1) * sizeof(cw_uint32_t));
		}
		(*a_nunchosen)--;

		/* Store this choice. */
		for (i = 0; i < nbits; i++)
		{
		    TR_BIT_SET(a_tr, *a_bitind_perm,
			       ((offset >> (nbits - i - 1) & 0x1)));
//		    fprintf(stderr, "%c",
//			    ((offset >> (nbits - i - 1) & 0x1)) ? '1' : '0');
		    (*a_bitind_perm)++;
		}
//		fprintf(stderr, " ");
	    }
	}
    }
    else
    {
	cw_bool_t did_paren = FALSE;

	/* Internal node. */

	for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	    {
		if (did_paren == FALSE)
		{
		    did_paren = TRUE;

		    /* Insert open paren. */
		    TR_BIT_SET(a_tr, *a_bitind_paren, 0);
		    (*a_bitind_paren)++;
//		    fprintf(stderr, "(");
		}

		/* Recurse. */
		tr_p_new_recurse(a_tr, a_bitind_paren, a_bitind_perm,
				 a_trn->neighbors[i], a_trn,
				 a_unchosen, a_nunchosen);
	    }
	}

	if (did_paren)
	{
	    /* Insert close paren. */
	    TR_BIT_SET(a_tr, *a_bitind_paren, 1);
	    (*a_bitind_paren)++;
//	    fprintf(stderr, ")");
	}
    }
}

cw_tr_t *
tr_new(cw_mema_t *a_mema, cw_trn_t *a_trn, cw_uint32_t a_ntaxa)
{
    cw_tr_t *retval;
    cw_trn_t *root, *node;
    cw_uint32_t tr_sizeof, bitind_paren, bitind_perm, i, *unchosen, nunchosen;

    cw_check_ptr(a_mema);
    cw_check_ptr(mema_alloc_get(a_mema));
    cw_dassert(trn_p_tree_validate(a_trn));
    cw_assert(trn_tree_ntaxa_get(a_trn) == a_ntaxa);

    tr_sizeof = tr_ntaxa2sizeof(a_ntaxa);
    retval = cw_opaque_alloc(mema_alloc_get(a_mema),
			     mema_arg_get(a_mema), tr_sizeof);

    /* Set the last byte to 0 to assure that there is no trailing garbage. */
    retval[tr_sizeof - 1] = 0;

    /* Canonize the tree. */
    root = trn_tree_root_get(a_trn);
    trn_p_tree_canonize(root, NULL);

    /* The tree can now be converted to parenthetical form and taxa permutation
     * using an in-order traversal. */
    bitind_paren = 0;
    bitind_perm = 2 * (a_ntaxa - 3) + 1;

    /* The first open paren is implied. */

    /* Create a taxon permutation from a_trn.  This requires maintaining a list
     * of the taxa that remain to be chosen from.  An array of taxon numbers is
     * maintained in compact form (already chosen taxa are removed, and trailing
     * space in the array is ignored). */
    unchosen = cw_opaque_alloc(mema_alloc_get(a_mema), mema_arg_get(a_mema),
			     (a_ntaxa - 2) * sizeof(cw_uint32_t));
    for (i = 0; i < (a_ntaxa - 2); i++)
    {
	unchosen[i] = i + 2;
    }
    nunchosen = a_ntaxa - 2;

    /* Recurse if the tree has an internal node.  Taxon 0 must be avoided
     * during the traversal, since it is implicit in the canonical tree
     * encoding.  Additionally, the first internal node must be handled here
     * rather than simply recursing to it, since the first '(' is implied. */
    node = root->neighbors[0];
    if (node != NULL && node->taxon_num == CW_TRN_TAXON_NONE)
    {
	for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    if (node->neighbors[i] != NULL && node->neighbors[i] != root)
	    {
		/* Recurse. */
		tr_p_new_recurse(retval, &bitind_paren, &bitind_perm,
				 node->neighbors[i], node,
				 unchosen, &nunchosen);
	    }
	}
    }

    /* Insert close paren. */
    TR_BIT_SET(retval, bitind_paren, 1);
//    fprintf(stderr, ")\n");
    cw_assert(bitind_paren == 2 * (a_ntaxa - 3));

    /* Clean up. */
    cw_opaque_dealloc(mema_dealloc_get(a_mema), mema_arg_get(a_mema),
		      unchosen, (a_ntaxa - 2) * sizeof(cw_uint32_t));

    return retval;
}

void
tr_delete(cw_tr_t *a_tr, cw_mema_t *a_mema, cw_uint32_t a_ntaxa)
{
    cw_check_ptr(a_mema);
    cw_check_ptr(mema_dealloc_get(a_mema));
    cw_dassert(tr_p_validate(a_tr));
    cw_assert(tr_p_sizeof(a_tr) == tr_ntaxa2sizeof(a_ntaxa));

    cw_opaque_dealloc(mema_dealloc_get(a_mema),
		      mema_arg_get(a_mema),
		      a_tr, tr_ntaxa2sizeof(a_ntaxa));
}

/* Determine the number of taxa in a_tr. */
cw_uint32_t
tr_ntaxa(const cw_tr_t *a_tr)
{
    cw_uint32_t retval;
    cw_uint32_t i, npairs, curdepth;

    /* Determine how many pairs of parentheses there are in the bit vector. */
    for (i = 0, npairs = curdepth = 1; curdepth > 0; i++)
    {
	if (TR_BIT_GET(a_tr, i))
	{
	    curdepth--;
	}
	else
	{
	    npairs++;
	    curdepth++;
	}
    }

    /* Using the first term of the formula at the top of the file, we know that:
     *
     *   nparens == 2(n-3) + 1
     *
     *   npairs == (nparens + 1) / 2
     *
     *   npairs == ((2(n-3) + 1) + 1) / 2
     *
     *   npairs == (2(n-3) + 2) / 2
     *
     *   npairs == (2(n-2)) / 2
     *
     *   npairs == n - 2
     *
     *   n == npairs + 2
     *
     * Solve for n (retval). */
    retval = npairs + 2;

    return retval;
}

/* Calculate bit vector byte count, given the number of taxa in the tree. */
cw_uint32_t
tr_ntaxa2sizeof(cw_uint32_t a_ntaxa)
{
    cw_uint32_t retval;
    cw_uint32_t ceil;

    /* Use the formula at the beginning of this file to calculate the number of
     * bytes required to store a tree with a_ntaxa. */

    retval = tr_p_ntaxa2nbits(a_ntaxa);

    /* Ceiling of denominator division. */
    ceil = !!(retval & 0x7);

    /* Divide. */
    retval >>= 3;
    retval += ceil;

    return retval;
}

void
tr_memcopy(cw_uint8_t *a_dest, cw_tr_t *a_tr)
{
    memcpy(a_dest, a_tr, tr_ntaxa2sizeof(tr_ntaxa(a_tr)));
}

cw_trn_t *
tr_trn(cw_tr_t *a_tr, cw_mema_t *a_mema, cw_uint32_t a_ntaxa)
{
    cw_dassert(tr_p_validate(a_tr));
    cw_assert(tr_ntaxa(a_tr) == a_ntaxa);

    /* Unimplemented for now, since it isn't needed for anything yet. */
    cw_error("XXX Not implemented");
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
    cw_dassert(tr_p_validate(tr));

    len = tr_p_sizeof(tr);
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
    cw_dassert(tr_p_validate(tr1));
    cw_dassert(tr_p_validate(tr2));

    len1 = tr_p_sizeof(tr1);
    len2 = tr_p_sizeof(tr2);
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
#ifdef CW_DBG
/* Validate an individual trn. */
static cw_bool_t
trn_p_validate(cw_trn_t *a_trn)
{
    cw_check_ptr(a_trn);
    cw_assert(a_trn->magic == CW_TRN_MAGIC);

    /* A trn must either be a leaf node or an internal node.  This check allows
     * an internal node to have less than 3 neighbors, since multiple calls are
     * necessary to set up all the neighbor pointers.  Likewise, a leaf node is
     * not required to have a neighbor. */
    if (a_trn->taxon_num != CW_TRN_TAXON_NONE)
    {
	cw_uint32_t i, nneighbors;

	for (i = nneighbors = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    if (a_trn->neighbors[i] != NULL)
	    {
		nneighbors++;
	    }
	}
	cw_assert(nneighbors <= 1);
    }

    return TRUE;
}
#endif

cw_trn_t *
trn_new(cw_trn_t *a_trn, cw_mema_t *a_mema)
{
    cw_trn_t *retval;

    if (a_trn != NULL)
    {
	retval = a_trn;
	memset(a_trn, 0, sizeof(cw_trn_t));
	a_trn->mema = NULL;
    }
    else
    {
	cw_check_ptr(a_mema);
	cw_check_ptr(mema_calloc_get(a_mema));
	cw_check_ptr(mema_dealloc_get(a_mema));

	retval = cw_opaque_calloc(mema_calloc_get(a_mema),
				  mema_arg_get(a_mema), 1, sizeof(cw_trn_t));
	a_trn->mema = a_mema;
    }

    retval->taxon_num = CW_TRN_TAXON_NONE;

#ifdef CW_DBG
    retval->magic = CW_TRN_MAGIC;
#endif

    return retval;
}

void
trn_delete(cw_trn_t *a_trn)
{
    cw_dassert(trn_p_validate(a_trn));

    if (a_trn->mema != NULL)
    {
	cw_opaque_dealloc(mema_dealloc_get(a_trn->mema),
			  mema_arg_get(a_trn->mema), a_trn, sizeof(cw_trn_t));
    }
#ifdef CW_DBG
    else
    {
	memset(a_trn, 0x5a, sizeof(cw_trn_t));
    }
#endif
}

cw_uint32_t
trn_taxon_num_get(cw_trn_t *a_trn)
{
    cw_dassert(trn_p_validate(a_trn));

    return a_trn->taxon_num;
}

void
trn_taxon_num_set(cw_trn_t *a_trn, cw_uint32_t a_taxon_num)
{
    cw_dassert(trn_p_validate(a_trn));

    a_trn->taxon_num = a_taxon_num;
}

cw_trn_t *
trn_neighbor_get(cw_trn_t *a_trn, cw_uint32_t a_i)
{
    cw_dassert(trn_p_validate(a_trn));
    cw_assert(a_i < CW_TRN_MAX_NEIGHBORS);

    return a_trn->neighbors[a_i];
}

void
trn_neighbors_swap(cw_trn_t *a_trn, cw_uint32_t a_i, cw_uint32_t a_j)
{
    cw_trn_t *t_trn;

    cw_dassert(trn_p_validate(a_trn));
    cw_assert(a_i < CW_TRN_MAX_NEIGHBORS);
    cw_assert(a_j < CW_TRN_MAX_NEIGHBORS);
    cw_assert(a_i != a_j);

    t_trn = a_trn->neighbors[a_i];
    a_trn->neighbors[a_i] = a_trn->neighbors[a_j];
    a_trn->neighbors[a_j] = t_trn;
}

void
trn_join(cw_trn_t *a_a, cw_trn_t *a_b)
{
    cw_uint32_t i, j;

    cw_dassert(trn_p_validate(a_a));
    cw_dassert(trn_p_validate(a_b));
    cw_assert(a_a != a_b);
#ifdef CW_DBG
    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	cw_assert(a_a->neighbors[i] != a_b);
	cw_assert(a_b->neighbors[i] != a_a);
    }
#endif

    /* Find an empty slot in a_a. */
    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_a->neighbors[i] == NULL)
	{
	    break;
	}
    }
    cw_assert(i < CW_TRN_MAX_NEIGHBORS);

    /* Find an empty slot in a_b. */
    for (j = 0; j < CW_TRN_MAX_NEIGHBORS; j++)
    {
	if (a_b->neighbors[j] == NULL)
	{
	    break;
	}
    }
    cw_assert(j < CW_TRN_MAX_NEIGHBORS);

    /* Join the two nodes. */
    a_a->neighbors[i] = a_b;
    a_b->neighbors[j] = a_a;

    cw_dassert(trn_p_validate(a_a));
    cw_dassert(trn_p_validate(a_b));
}

void
trn_detach(cw_trn_t *a_a, cw_trn_t *a_b)
{
    cw_uint32_t i, j;

    cw_dassert(trn_p_validate(a_a));
    cw_dassert(trn_p_validate(a_b));

    /* Find the slot in a_a that points to a_b. */
    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_a->neighbors[i] == a_b)
	{
	    break;
	}
    }
    cw_assert(i < CW_TRN_MAX_NEIGHBORS);

    /* Find the slot in a_b that points to a_a. */
    for (j = 0; j < CW_TRN_MAX_NEIGHBORS; j++)
    {
	if (a_b->neighbors[j] == a_a)
	{
	    break;
	}
    }
    cw_assert(j < CW_TRN_MAX_NEIGHBORS);

    /* Detach the two nodes. */
    a_a->neighbors[i] = NULL;
    a_b->neighbors[j] = NULL;

    cw_dassert(trn_p_validate(a_a));
    cw_dassert(trn_p_validate(a_b));
}

void *
trn_aux_get(cw_trn_t *a_trn)
{
    cw_dassert(trn_p_validate(a_trn));

    return a_trn->aux;
}

void
trn_aux_set(cw_trn_t *a_trn, void *a_aux)
{
    cw_dassert(trn_p_validate(a_trn));

    a_trn->aux = a_aux;
}

/* trn_tree. */
#ifdef CW_DBG
/* Return the number of taxa with number a_taxon_num in the subtree rooted at
 * a_trn. */
static cw_uint32_t
trn_p_tree_validate_recurse(cw_trn_t *a_trn, cw_trn_t *a_prev,
			    cw_uint32_t a_taxon_num)
{
    cw_uint32_t retval;
    cw_uint32_t i;

    if (a_trn->taxon_num != CW_TRN_TAXON_NONE)
    {
	/* Leaf node. */
	cw_assert(a_trn->neighbors[0] != NULL);
	for (i = 1; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    cw_assert(a_trn->neighbors[i] == NULL);
	}

	if (a_trn->taxon_num == a_taxon_num)
	{
	    retval = 1;
	}
	else
	{
	    retval = 0;
	}
    }
    else
    {
	/* Internal node. */
	for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
	{
	    cw_assert(a_trn->neighbors[i] != NULL);
	}

	retval = 0;
    }

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    retval += trn_p_tree_validate_recurse(a_trn->neighbors[i],
						  a_trn, a_taxon_num);
	}
    }

    return retval;
}

/* Validate a tree constructed with trn's. */
static cw_bool_t
trn_p_tree_validate(cw_trn_t *a_trn)
{
    cw_trn_t *root;
    cw_uint32_t i, j, n, ntaxa;

    trn_p_validate(a_trn);

    /* Find the root. */
    root = trn_p_tree_root_get(a_trn, NULL, NULL);
    cw_check_ptr(root);

    /* Traverse the tree, and make sure that the following invariants hold:
     *
     * + Leaf nodes have a taxon number and precisely 1 neighbor.
     *
     * + Internal nodes have no taxon number and precisely 3 neighbors.
     *
     * + Each taxon number appears no more than once in the tree.
     *
     * These invariants allow gaps in the taxon numbering, which has the
     * potential to cause problems (not in the tr implementation but elswhere),
     * but requiring contiguous taxon numbering would make validating tree
     * bisections impossible.
     */
    for (i = j = 0, ntaxa = trn_p_tree_ntaxa_get(root, NULL);
	 j < ntaxa;
	 i++, j += n)
    {
	n = trn_p_tree_validate_recurse(root, NULL, i);
	cw_assert(n <= 1);
    }

    return TRUE;
}
#endif

static void
trn_p_tree_delete(cw_trn_t *a_trn)
{
    cw_uint32_t i;
    cw_trn_t *trn;

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL)
	{
	    trn = a_trn->neighbors[i];
	    trn_detach(a_trn, trn);
	    trn_p_tree_delete(trn);
	}
    }

    trn_delete(a_trn);
}

/* Generate a random tree. */
/* XXX Take a function pointer to a random number function. */
cw_trn_t *
trn_tree_random(cw_mema_t *a_mema, cw_uint32_t a_ntaxa)
{
    cw_trn_t *trn, **subtrees;
    cw_uint32_t i, a, b, t;

    cw_assert(a_ntaxa >= 2);

    /* Allocate an array that is large enough to hold a pointer to each taxon
     * trn, and populated it. */
    subtrees = (cw_trn_t **) cw_opaque_alloc(mema_alloc_get(a_mema),
					     mema_arg_get(a_mema),
					     a_ntaxa * sizeof(cw_trn_t *));
    for (i = 0; i < a_ntaxa; i++)
    {
	subtrees[i] = trn_new(NULL, a_mema);
	trn_taxon_num_set(subtrees[i], i);
    }

    /* Iteratively randomly select two items from the array, join them, and
     * insert the result back into the array.  Stop when there are two subtrees
     * left, and join them directly, in order to create an unrooted tree as the
     * final result. */
    for (i = a_ntaxa; i > 2; i--)
    {
	/* Choose two elements randomly, such that a < b. */
	a = random() % i;
	b = random() % (i - 1);
	if (b >= a)
	{
	    b++;
	}
	cw_assert(a != b);
	if (a > b)
	{
	    t = a;
	    a = b;
	    b = t;
	}

	/* Allocate a new internal node and join. */
	trn = trn_new(NULL, a_mema);
	trn_join(trn, subtrees[a]);
	trn_join(trn, subtrees[b]);

	/* Insert the new subtree and shorten the array. */
	subtrees[a] = trn;
	if (b < (i - 1))
	{
	    subtrees[b] = subtrees[i - 1];
	}
    }

    /* Join the final two subtrees. */
    trn_join(subtrees[0], subtrees[1]);
    trn = trn_tree_root_get(subtrees[0]);

    /* Clean up. */
    cw_opaque_dealloc(mema_dealloc_get(a_mema), mema_arg_get(a_mema),
		      subtrees, a_ntaxa * sizeof(cw_trn_t *));

    return trn;
}

void
trn_tree_delete(cw_trn_t *a_trn)
{
    cw_dassert(trn_p_validate(a_trn));

    trn_p_tree_delete(a_trn);
}

/* Recursively traverse the tree and count the number of taxa. */
static cw_uint32_t
trn_p_tree_ntaxa_get(cw_trn_t *a_trn, cw_trn_t *a_prev)
{
    cw_uint32_t retval;
    cw_uint32_t i;

    cw_dassert(trn_p_validate(a_trn));

    if (a_trn->taxon_num != CW_TRN_TAXON_NONE)
    {
	/* Leaf node. */
	retval = 1;
    }
    else
    {
	/* Internal node. */
	retval = 0;
    }

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    retval += trn_p_tree_ntaxa_get(a_trn->neighbors[i], a_trn);
	}
    }

    return retval;
}

cw_uint32_t
trn_tree_ntaxa_get(cw_trn_t *a_trn)
{
    cw_dassert(trn_p_tree_validate(a_trn));

    return trn_p_tree_ntaxa_get(a_trn, NULL);
}

cw_uint32_t
trn_tree_nedges_get(cw_trn_t *a_trn)
{
    cw_dassert(trn_p_tree_validate(a_trn));

    return ((trn_p_tree_ntaxa_get(a_trn, NULL) * 2) - 3);
}

/* Recursively traverse the tree and find the lowest numbered taxon. */
static cw_trn_t *
trn_p_tree_root_get(cw_trn_t *a_trn, cw_trn_t *a_prev, cw_trn_t *a_root)
{
    cw_trn_t *retval, *root, *troot;
    cw_uint32_t i;

    if (a_trn->taxon_num != CW_TRN_TAXON_NONE
	&& (a_root == NULL || a_trn->taxon_num < a_root->taxon_num))
    {
	retval = a_trn;
	root = a_trn;
    }
    else
    {
	retval = NULL;
	root = a_root;
    }

    /* Iterate over neighbors. */
    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    troot = trn_p_tree_root_get(a_trn->neighbors[i], a_trn, root);
	    if (troot != NULL)
	    {
		retval = troot;
		root = troot;
	    }
	}
    }

    return retval;
}

cw_trn_t *
trn_tree_root_get(cw_trn_t *a_trn)
{
    cw_dassert(trn_p_tree_validate(a_trn));

    return trn_p_tree_root_get(a_trn, NULL, NULL);
}

/* Convert a tree to canonical form by re-ordering the neighbors array such that
 * subtrees are in increasing order of minimum taxon number contained. */
static cw_uint32_t
trn_p_tree_canonize(cw_trn_t *a_trn, cw_trn_t *a_prev)
{
    cw_uint32_t retval;
    cw_uint32_t i, j, t;
    cw_uint32_t subtree_mins[CW_TRN_MAX_NEIGHBORS - 1];
    cw_uint32_t subtree_inds[CW_TRN_MAX_NEIGHBORS - 1];
    cw_bool_t swapped;

    if (a_trn->taxon_num != CW_TRN_TAXON_NONE)
    {
	/* Leaf node. */
	retval = a_trn->taxon_num;
    }
    else
    {
	/* Internal node. */
	retval = CW_TRN_TAXON_NONE;
    }

    /* Iteratively canonize subtrees, keeping track of the minimum taxon
     * number seen overall, as well as for each subtree. */
    for (i = j = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    cw_assert(j < (CW_TRN_MAX_NEIGHBORS - 1));
	    subtree_mins[j] = trn_p_tree_canonize(a_trn->neighbors[i], a_trn);
	    if (subtree_mins[j] < retval)
	    {
		retval = subtree_mins[j];
	    }
	    subtree_inds[j] = i;
	    j++;
	}
    }

    /* Bubble sort the subtrees.  This algorithm works in the general case, and
     * in the case this code is actually designed for (bifurcating trees), it
     * only requires a couple of extra branches. */
    do
    {
	swapped = FALSE;

	for (i = 0; i + 1 < j; i++)
	{
	    if (subtree_mins[i] > subtree_mins[i + 1])
	    {
		swapped = TRUE;

		/* Swap subtrees. */
		trn_neighbors_swap(a_trn, subtree_inds[i], subtree_inds[i + 1]);

		/* Swap subtree_* arrays. */
		t = subtree_mins[i];
		subtree_mins[i] = subtree_mins[i + 1];
		subtree_mins[i + 1] = t;
	    
		t = subtree_inds[i];
		subtree_inds[i] = subtree_inds[i + 1];
		subtree_inds[i + 1] = t;
	    }
	}
    } while (swapped);

    return retval;
}

static cw_bool_t
trn_p_tree_edge_get_recurse(cw_trn_t *a_trn, cw_uint32_t a_edge,
			    cw_trn_t *a_prev, cw_uint32_t *r_edge_count,
			    cw_trn_t **r_trn, cw_uint32_t *r_neighbor)
{
    cw_bool_t retval;
    cw_uint32_t i;

    for (i = 0; i < CW_TRN_MAX_NEIGHBORS; i++)
    {
	if (a_trn->neighbors[i] != NULL && a_trn->neighbors[i] != a_prev)
	{
	    /* Increment edge count before recursing.  If the edge count has
	     * reached the desired value, return this trn and neighbor index,
	     * and terminate recursion. */
	    *r_edge_count++;
	    if (*r_edge_count > a_edge)
	    {
		cw_assert(*r_edge_count == a_edge + 1);
		*r_trn = a_trn;
		*r_neighbor = i;

		retval = TRUE;
		goto RETURN;
	    }

	    /* Iteratively recurse into neighbor subtrees. */
	    if (trn_p_tree_edge_get_recurse(a_trn->neighbors[i], a_edge,
					    a_trn, r_edge_count,
					    r_trn, r_neighbor))
	    {
		retval = TRUE;
		goto RETURN;
	    }
	}
    }

    retval = FALSE;
    RETURN:
    return retval;
}

/* Do an in-order traversal of the tree and return the node that neighbors the
 * edge being sought, along with which neighbor of that node the edge is
 * between. */
static void
trn_p_tree_edge_get(cw_trn_t *a_trn, cw_uint32_t a_edge, cw_trn_t **r_trn,
		    cw_uint32_t *r_neighbor)
{
    cw_uint32_t edge_count = 0;

    trn_p_tree_edge_get_recurse(a_trn, a_edge, NULL, &edge_count,
				r_trn, r_neighbor);
}


/* a_trn must be the root of the tree. */
void
trn_tree_bisect(cw_trn_t *a_trn, cw_uint32_t a_edge, cw_trn_t **r_trn_a,
		cw_trn_t **r_trn_b, cw_trn_t **r_spare_a, cw_trn_t **r_spare_b)
{
    cw_trn_t *trn_a, *trn_b;
    cw_uint32_t edge;

    cw_dassert(trn_p_tree_validate(a_trn));
    cw_assert(trn_p_tree_root_get(a_trn, NULL, NULL) == a_trn);

    /* Get the nodes to either side of the edge where the bisection will be
     * done. */
    trn_p_tree_edge_get(a_trn, a_edge, &trn_a, &edge);
    trn_b = trn_a->neighbors[edge];

    /* Detach the two nodes. */
    trn_detach(trn_a, trn_b);

    /* There are two cases possible for each of the nodes.  Each is either a
     * leaf node or an internal node.  For a leaf node, do nothing.  For an
     * internal node, join its neighbors together, and return the node as a
     * spare. */

    /* trn_a. */
    if (trn_a->taxon_num == CW_TRN_TAXON_NONE)
    {
	cw_trn_t *a, *b;
	cw_uint32_t i;

	/* Get trn_a's neighbors. */
	for (i = 0, a = b = NULL; b == NULL; i++)
	{
	    cw_assert(i < CW_TRN_MAX_NEIGHBORS);

	    if (trn_a->neighbors[i] != NULL)
	    {
		if (a == NULL)
		{
		    a = trn_a->neighbors[i];
		}
		else
		{
		    b = trn_a->neighbors[i];
		}
	    }
	}

	/* Detach. */
	trn_detach(trn_a, a);
	trn_detach(trn_a, b);

	/* Join. */
	trn_join(a, b);

	*r_spare_a = trn_a;
    }
    else
    {
	*r_spare_a = NULL;
    }

    /* trn_b. */
    if (trn_b->taxon_num == CW_TRN_TAXON_NONE)
    {
	cw_trn_t *a, *b;
	cw_uint32_t i;

	/* Get trn_b's neighbors. */
	for (i = 0, a = b = NULL; b == NULL; i++)
	{
	    cw_assert(i < CW_TRN_MAX_NEIGHBORS);

	    if (trn_b->neighbors[i] != NULL)
	    {
		if (a == NULL)
		{
		    a = trn_b->neighbors[i];
		}
		else
		{
		    b = trn_b->neighbors[i];
		}
	    }
	}

	/* Detach. */
	trn_detach(trn_b, a);
	trn_detach(trn_b, b);

	/* Join. */
	trn_join(a, b);

	*r_spare_b = trn_b;
    }
    else
    {
	*r_spare_b = NULL;
    }

    /* Move *r_spare_b to *r_spare_a if *r_spare_a is NULL. */
    if (*r_spare_a == NULL)
    {
	*r_spare_a = *r_spare_b;
	*r_spare_b = NULL;
    }
}

void
trn_tree_connect(cw_trn_t *a_trn_a, cw_uint32_t a_edge_a,
		 cw_trn_t *a_trn_b, cw_uint32_t a_edge_b,
		 cw_trn_t *a_spare_a, cw_trn_t *a_spare_b)
{
    cw_dassert(trn_p_tree_validate(a_trn_a));
    cw_assert(a_edge_a < trn_p_tree_ntaxa_get(a_trn_a, NULL));
    cw_dassert(trn_p_tree_validate(a_trn_b));
    cw_assert(a_edge_b < trn_p_tree_ntaxa_get(a_trn_b, NULL));

    cw_error("XXX not implemented");
}
