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
CxTreeMpPrepare(CxtCxtTreeObject *self, PyObject *args);

PyObject *
CxTreeMpFinish(CxtCxtTreeObject *self);

PyObject *
CxTreeMp(CxtCxtTreeObject *self);

PyObject *
CxTreeTbrBestNeighborsMp(CxtCxtTreeObject *self, PyObject *args);

PyObject *
CxTreeTbrBetterNeighborsMp(CxtCxtTreeObject *self, PyObject *args);

PyObject *
CxTreeTbrAllNeighborsMp(CxtCxtTreeObject *self);

PyObject *
CxTreeNheldGet(CxtCxtTreeObject *self);

PyObject *
CxTreeheldGet(CxtCxtTreeObject *self, PyObject *args);
