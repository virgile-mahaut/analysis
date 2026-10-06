#include "SiliconCaloTrack_v1.h"

#include <phool/PHObject.h>  // for PHObject

#include <climits>
#include <map>
#include <vector>  // for vector

SiliconCaloTrack_v1::SiliconCaloTrack_v1(const SiliconCaloTrack& source)
{
  SiliconCaloTrack_v1::CopyFrom(source);
}

// have to suppress missingMemberCopy from cppcheck, it does not
// go down to the CopyFrom method where things are done correctly
// cppcheck-suppress missingMemberCopy
SiliconCaloTrack_v1::SiliconCaloTrack_v1(const SiliconCaloTrack_v1& source)
  : SiliconCaloTrack(source)
{
  SiliconCaloTrack_v1::CopyFrom(source);
}

SiliconCaloTrack_v1& SiliconCaloTrack_v1::operator=(const SiliconCaloTrack_v1& source)
{
  if (this != &source)
  {
    CopyFrom(source);
  }
  return *this;
}

void SiliconCaloTrack_v1::CopyFrom(const SiliconCaloTrack& source)
{
  // do nothing if copying onto oneself
  if (this == &source)
  {
    return;
  }

  // parent class method
  SiliconCaloTrack::CopyFrom(source);

  _track_id = source.get_id();
  _calo_id  = source.get_calo_id();

  set_pt(source.get_pt());
  set_phi(source.get_phi());
  set_eta(source.get_eta());

  set_cal_dphi(source.get_cal_dphi());
  set_cal_dz(  source.get_cal_dz());

}

void SiliconCaloTrack_v1::identify(std::ostream& os) const
{
  os << "SiliconCaloTrack_v1 Object ";
  os << "id: " << get_id() << " ";
  os << "calo id: " << get_calo_id() << " ";
  os << std::endl;

  os << "(px,py,pz) = ("
     << get_px() << ","
     << get_py() << ","
     << get_pz() << ")" << std::endl;

  //os << "(x,y,z) = (" << get_x() << "," << get_y() << "," << get_z() << ")" << std::endl;

  return;
}

int SiliconCaloTrack_v1::isValid() const
{
  return 1;
}


