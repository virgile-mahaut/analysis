#ifndef SILICONCALOTRACK_H
#define SILICONCALOTRACK_H

#include <phool/PHObject.h>

#include <limits.h>
#include <cmath>
#include <iostream>
#include <map>
#include <set>

class SiliconCaloTrack : public PHObject
{
 public:
  ~SiliconCaloTrack() override = default;

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override
  {
    os << "SiliconCaloTrack base class" << std::endl;
  }

  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  //! import PHObject CopyFrom, in order to avoid clang warning
  using PHObject::CopyFrom;

  //! copy content from base class
  virtual void CopyFrom(const SiliconCaloTrack&)
  {
  }

  //! copy content from base class
  virtual void CopyFrom(SiliconCaloTrack*)
  {
  }

  //
  // basic track information ---------------------------------------------------
  //

  virtual unsigned int get_id() const { return UINT_MAX; }
  virtual void set_id(unsigned int) {}

  virtual unsigned int get_calo_id() const { return UINT_MAX; }
  virtual void set_calo_id(unsigned int) {}

  virtual float get_x() const { return NAN; }
  virtual void set_x(float) {}

  virtual float get_y() const { return NAN; }
  virtual void set_y(float) {}

  virtual float get_z() const { return NAN; }
  virtual void set_z(float) {}

  virtual float get_pos(unsigned int) const { return NAN; }

  virtual float get_pt() const { return NAN; }
  virtual void  set_pt(float ) {}
  virtual float get_eta() const { return NAN; }
  virtual void  set_eta(float ) { }
  virtual float get_phi() const { return NAN; }
  virtual void  set_phi(float ) { }

  virtual float get_px() const { return NAN; }
  virtual float get_py() const { return NAN; }
  virtual float get_pz() const { return NAN; }
  virtual float get_p() const { return NAN; }


  virtual float get_error(int /*i*/, int /*j*/) const { return NAN; }
  virtual void set_error(int /*i*/, int /*j*/, float /*value*/) {}


  // cal info
  virtual float get_cal_dphi() const { return 0.; }
  virtual void  set_cal_dphi(float /*dphi*/) {}
  virtual float get_cal_dz() const { return 0.; }
  virtual void  set_cal_dz(float /*z*/) {}

 protected:
  SiliconCaloTrack() {}

  ClassDefOverride(SiliconCaloTrack, 1);
};

#endif
