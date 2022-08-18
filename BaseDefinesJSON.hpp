#ifndef DEF_BASEDEF_HEADER_DEF_JSONCPP
#define DEF_BASEDEF_HEADER_DEF_JSONCPP



/*
nlohmann/json.hpp　を c_jsonとして活用
*/
#include "000_thirdparty/nlohmann/json.hpp"

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{

using c_json = nlohmann::json;

/* end of namespace */
};


#endif
