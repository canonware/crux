//==============================================================================
//
// <Copyright = jasone>
// <License>
//
//==============================================================================
//
// Version: Crux <Version = crux>
//
//==============================================================================

PyObject *
CxTreeMpPrepare(CxtTreeObject *self, PyObject *args);

PyObject *
CxTreeMpFinish(CxtTreeObject *self);

PyObject *
CxTreeMp(CxtTreeObject *self);

PyObject *
CxTreeTbrBestNeighborsMp(CxtTreeObject *self, PyObject *args);

PyObject *
CxTreeTbrBetterNeighborsMp(CxtTreeObject *self, PyObject *args);

PyObject *
CxTreeTbrAllNeighborsMp(CxtTreeObject *self);

PyObject *
CxTreeNHeldGet(CxtTreeObject *self);

PyObject *
CxTreeHeldGet(CxtTreeObject *self, PyObject *args);
