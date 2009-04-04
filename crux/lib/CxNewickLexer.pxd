cdef extern from "CxNewickLexer.h":
    ctypedef void *yyscan_t
    ctypedef enum CxtNewickLexerInputMode:
        CxeNewickLexerInputModeFd
        CxeNewickLexerInputModeStr
    ctypedef enum CxtNewickLexerTokType:
        CxeNewickLexerTokTypeComment
        CxeNewickLexerTokTypeLparen
        CxeNewickLexerTokTypeRparen
        CxeNewickLexerTokTypeComma
        CxeNewickLexerTokTypeColon
        CxeNewickLexerTokTypeSemicolon
        CxeNewickLexerTokTypeBranchLength
        CxeNewickLexerTokTypeUnquotedLabel
        CxeNewickLexerTokTypeQuotedLabel
        CxeNewickLexerTokTypeWhitespace
    cdef struct struct_input_s:
        char *s
        size_t len
    cdef union union_in:
        int fd
        struct_input_s s
    cdef struct struct_tok:
        char *s
        size_t len
    ctypedef struct CxtNewickLexerExtra:
        unsigned line
        int ioerror
        CxtNewickLexerInputMode inputMode
        union_in input
        CxtNewickLexerTokType tokType
        struct_tok tok

    cdef int CxNewickLexer_lex_init_extra(CxtNewickLexerExtra *extra, \
      yyscan_t *scanner)
    cdef int CxNewickLexer_lex_destroy(yyscan_t scanner)
    cdef int CxNewickLexer_lex(yyscan_t scanner)
    cdef char *CxNewickLexer_get_text(yyscan_t scanner)
