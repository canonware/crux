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
CxTreeTbrBEdgeSetsGet(CxtTreeObject *self, CxtEdgeObject *aEdge,
		      CxtEdgeObject ***rSetA, unsigned *rNSetA,
		      CxtRingObject **rRingA,
		      CxtEdgeObject ***rSetB, unsigned *rNSetB,
		      CxtRingObject **rRingB);

bool
CxTreeTbr(CxtTreeObject *self, CxtEdgeObject *aBisect,
	  CxtEdgeObject *aReconnectA, CxtEdgeObject *aReconnectB);
PyObject *
CxTreeTbrPargs(CxtTreeObject *self, PyObject *args);

bool
CxTreeTbrNEdgesGet(CxtTreeObject *self, unsigned *rNEdges);

bool
CxTreeTbrEdgeGet(CxtTreeObject *self, unsigned aEdge, CxtEdgeObject **rEdge);

bool
CxTreeTbrEdgeOffset(CxtTreeObject *self, unsigned aEdge, unsigned *rOffset);

bool
CxTreeTbrNNeighborsGet(CxtTreeObject *self, unsigned *rNNeighbors);
PyObject *
CxTreeTbrNNeighborsGetPargs(CxtTreeObject *self);

bool
CxTreeTbrNeighborGet(CxtTreeObject *self, unsigned aNeighbor,
		     CxtEdgeObject **rBisect, CxtEdgeObject **rReconnectA,
		     CxtEdgeObject **rReconnectB);
PyObject *
CxTreeTbrNeighborGetPargs(CxtTreeObject *self, PyObject *args);
