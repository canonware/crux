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

bool
CxTreeTbr(CxtTreeObject *self, CxtEdgeObject *aBisect,
	  CxtEdgeObject *aReconnectA, CxtEdgeObject *aReconnectB);
PyObject *
CxTreeTbrPargs(CxtTreeObject *self, PyObject *args);

bool
CxTreeTbrNneighborsGet(CxtTreeObject *self, unsigned *rNneighbors);
PyObject *
CxTreeTbrNneighborsGetPargs(CxtTreeObject *self);

bool
CxTreeTbrNeighborGet(CxtTreeObject *self, unsigned aNeighbor,
		     CxtEdgeObject **rBisect, CxtEdgeObject **rReconnectA,
		     CxtEdgeObject **rReconnectB);
PyObject *
CxTreeTbrNeighborGetPargs(CxtTreeObject *self, PyObject *args);
