#ifndef CORELIB___NCBIDBG__HPP
#define CORELIB___NCBIDBG__HPP

#  define NCBI_ASSERT(expr, mess)   ((void)0)

#define _ASSERT(expr)   NCBI_ASSERT(expr, NULL)

#endif  /* CORELIB___NCBIDBG__HPP */
