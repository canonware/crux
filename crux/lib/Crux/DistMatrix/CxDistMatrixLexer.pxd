cdef extern from "CxDistMatrixLexer.h":
    ctypedef void *yyscan_t
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
        float dist
        char *label
    ctypedef struct CxtDistMatrixLexerExtra:
        unsigned line
        int ioerror
        CxtDistMatrixLexerInputMode inputMode
        union_in input
        CxtDistMatrixLexerTokType tokType
        union_tok tok

    cdef int CxDistMatrixLexer_lex_init_extra(CxtDistMatrixLexerExtra *extra, \
      yyscan_t *scanner)
    cdef int CxDistMatrixLexer_lex_destroy(yyscan_t scanner)
    cdef int CxDistMatrixLexer_lex(yyscan_t scanner)
    cdef char *CxDistMatrixLexer_get_text(yyscan_t scanner)
