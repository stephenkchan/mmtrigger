// File gen_vd.C
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<int> >+;
#else
template class std::vector<std::vector<float> >;
template class std::vector<std::vector<int> >;
#endif
