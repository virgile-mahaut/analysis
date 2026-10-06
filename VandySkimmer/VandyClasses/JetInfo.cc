#include "JetInfo.h"

JetInfo::~JetInfo() = default;

void JetInfo::CopyTo(JetInfo* jet)
{
    jet->set_px(Px);
    jet->set_py(Py);
    jet->set_pz(Pz);
    jet->set_e(E);
    jet->set_pt(Pt);
    jet->set_pt_uncalib(Pt_uncalib);
    jet->set_hCaloFrac(HCaloFrac);
    std::vector<int> constit;
    for(auto cons : constituents)
    {
        constit.push_back(cons);
    }
    jet->set_constituents(constit);
    jet->set_parent(Parent);
    jet->set_isQuark(IsQuark);
    jet->set_parentPID(Parent_pID);

    return;
}
