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

// Create a tree from a symmetric pair-wise distance matrix, using the
// neighbor joining (NJ) algorithm.
PyObject *
CxDistMatrixNj(CxtDistMatrixObject *self, PyObject *args);

// Create a tree from a symmetric pair-wise distance matrix, using the
// relaxed neighbor joining (RNJ) algorithm.
PyObject *
CxDistMatrixRnj(CxtDistMatrixObject *self, PyObject *args);
