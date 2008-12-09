cdef extern from "sys/types.h":
    ctypedef unsigned long size_t

cdef extern from "CxDistMatrix.h":
    ctypedef double CxtDMDist
    ctypedef unsigned long CxtDMSize
    ctypedef enum CxtDistMatrixLexerInputMode:
        CxeDistMatrixLexerInputModeFd
        CxeDistMatrixLexerInputModeStr
    ctypedef enum CxtDistMatrixLexerTokType:
        CxeDistMatrixLexerTokTypeInt
        CxeDistMatrixLexerTokTypeDist
        CxeDistMatrixLexerTokTypeLabel
    cdef struct struct_input_s:
        char *s
        size_t len
    cdef union union_in:
        int fd
        struct_input_s s
    cdef union union_tok:
        unsigned long int_
        CxtDMDist dist
        char *label
    ctypedef struct CxtDistMatrixLexerExtra:
        unsigned line
        int ioerror
        CxtDistMatrixLexerInputMode inputMode
        union_in input
        CxtDistMatrixLexerTokType tokType
        union_tok tok

    cdef inline CxtDMSize CxDistMatrixNxy2i(CxtDMSize n, CxtDMSize x,
      CxtDMSize y)
    cdef CxtDMDist *CxDistMatrixNew(CxtDMSize ntaxa)
    cdef void CxDistMatrixDelete(CxtDMDist *matrix)

cdef extern from "CxDistMatrixLexer.h":
    ctypedef void *yyscan_t
    cdef int CxDistMatrixLexer_lex_init_extra(CxtDistMatrixLexerExtra *extra, \
      yyscan_t *scanner)
    cdef int CxDistMatrixLexer_lex_destroy(yyscan_t scanner)
    cdef int CxDistMatrixLexer_lex(yyscan_t scanner)
