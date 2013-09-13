#include "config.h"
#ifdef HAVE_UNORDERED_MAP
# include <unordered_map>
# define OUR_UNORDERED_MAP std::unordered_map
#elif defined ( HAVE_TR1_UNORDERED_MAP )
# include <tr1/unordered_map>
# define OUR_UNORDERED_MAP std::tr1::unordered_map
#else
# define OUR_UNORDERED_MAP std::map
#endif
