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
typedef struct cw_tri_s cw_tri_t;
typedef struct cw_trl_s cw_trl_t;

/* Unrooted bifurcating phylogenetic tree.  The feasible duration of each
 * iteration of tabu search is bounded by the number of trees (and assoicated
 * tabu lists) that can be stored.  As such, the internal tree representation is
 * as compact as reasonably possible without resorting to
 * compression/decompression.
 *
 * Since tabu search depends on being able to determine if a tree has been
 * previously visited, tree comparison must be fast.  Therefore, a canonical
 * form is defined for trees.  Since the internal representation is effectively
 * string, it can be used as the input to a string hashing function, which makes
 * searching for a specific tree among previously visited trees a constant time
 * operation.
 *
 * A tree is stored as a list of edges.  The leaf edges come first (in
 * taxon-sorted order), followed by the internal edges.  Each edge specifies
 * which other edges it is attached to, which means that efficient traversal in
 * any direction is possible (i.e. the data structure is doubly linked).
 *
 * The following example is in canonical form:
 *
 *            A           C
 *             \         /
 *              \0     2/
 *               \     /
 *                \   /
 *                 \ /
 *                  V
 *                  |
 *                  |
 *                  |6
 *                  |
 *                  |
 *                 / \
 *                /   \
 *              8/     \7
 *        5     /       \     1
 *   F---------/         \---------B
 *             |         |
 *             |         |
 *             |3       4|
 *             |         |
 *             |         |
 *             D         E
 *
 *   Taxon | Edge | Neigbors
 *   ------+------+---------
 *       A |    0 | 6 2
 *       B |    1 | 7 4
 *       C |    2 | 0 6
 *       D |    3 | 8 5
 *       E |    4 | 1 7
 *       F |    5 | 8 3
 *         |    6 | 0 2 7 8
 *         |    7 | 6 8 1 4
 *         |    8 | 6 7 3 5
 *
 * In the above table, only the neighbors column is actually stored.  There are
 * (2n - 3) edges, of which the first n are listed first in taxon-sorted order,
 * and the rest are listed in in-order traversal order (starting from the first
 * taxon).  The ordering of the neighbor lists are in an order defined by the
 * following rules:
 *
 *   1) The neighbors attached to one end of an edge are listed consecutively.
 *
 *   2) The order of the consecutive listing defined in step 1 is determined by
 *      which neighbor is visited first in an in-order traversal of the tree
 *      (starting from the first taxon).
 *
 * The neighors list is stored as an array of numbers, where each number uses
 * log_2(2n - 3) bits of memory.  The total number of bytes necessary for a tree
 * with n taxa is given by the formula:
 *
 *   ceil(ceil(log_2(2n - 3)) x (6n - 12) / 8)
 *
 * So for example, the above example would require 96 bits (12 bytes) of memory
 * to store the full tree.  Storing a tree with 500 taxa would require 3735
 * bytes of memory.
 */
struct cw_tr_s
{
#ifdef CW_DBG
    cw_uint32_t magic;
#define CW_TR_MAGIC 0x37d478a3
#endif

};
