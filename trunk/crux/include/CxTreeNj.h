/******************************************************************************
 *
 * <Copyright = jasone>
 * <License>
 *
 ******************************************************************************
 *
 * Version: Crux <Version = crux>
 *
 ******************************************************************************/

/* Create a tree from a symmetric pair-wise distance matrix, using the
 * neighbor joining (NJ) algorithm. */
PyObject *
CxTreeNj(CxtTreeObject *self, PyObject *args);

/* Create a tree from a symmetric pair-wise distance matrix, using the
 * relaxed neighbor joining (RNJ) algorithm. */
PyObject *
CxTreeRnj(CxtTreeObject *self, PyObject *args);
