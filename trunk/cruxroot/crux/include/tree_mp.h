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

PyObject *
tree_mp_prepare(TreeObject *self, PyObject *args);

PyObject *
tree_mp_finish(TreeObject *self);

PyObject *
tree_mp(TreeObject *self);

PyObject *
tree_tbr_best_neighbors_mp(TreeObject *self, PyObject *args);

PyObject *
tree_tbr_better_neighbors_mp(TreeObject *self, PyObject *args);

PyObject *
tree_tbr_all_neighbors_mp(TreeObject *self);

PyObject *
tree_nheld_get(TreeObject *self);

PyObject *
tree_held_get(TreeObject *self, PyObject *args);
