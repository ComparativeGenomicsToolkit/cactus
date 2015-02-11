#ifndef CORELIB___NCBISTL__HPP
#define CORELIB___NCBISTL__HPP

#include <common/ncbi_export.h>

#define BEGIN_SCOPE(ns) namespace ns {

#define END_SCOPE(ns) }

#define USING_SCOPE(ns) using namespace ns

#define NCBI_NS_STD  std

#define NCBI_USING_NAMESPACE_STD using namespace NCBI_NS_STD

#define NCBI_NS_NCBI ncbi

#define BEGIN_NCBI_SCOPE BEGIN_SCOPE(NCBI_NS_NCBI)

#define END_NCBI_SCOPE   END_SCOPE(NCBI_NS_NCBI)

#define USING_NCBI_SCOPE USING_SCOPE(NCBI_NS_NCBI)

namespace NCBI_NS_NCBI { /* the fake one, +"std" */ NCBI_USING_NAMESPACE_STD; }

#endif /* NCBISTL__HPP */
