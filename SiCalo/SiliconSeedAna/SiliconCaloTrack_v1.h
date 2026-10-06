#ifndef SILICONCALOTRACKV1_H
#define SILICONCALOTRACKV1_H

#include "SiliconCaloTrack.h"

#include <cmath>
#include <cstddef>  // for size_t
#include <iostream>
#include <map>
#include <utility>  // for pair

class PHObject;

class SiliconCaloTrack_v1 : public SiliconCaloTrack
{
 public:
  SiliconCaloTrack_v1(){};

  //* base class copy constructor
  SiliconCaloTrack_v1(const SiliconCaloTrack&);

  //* copy constructor
  SiliconCaloTrack_v1(const SiliconCaloTrack_v1& track);

  //* assignment operator
  SiliconCaloTrack_v1& operator=(const SiliconCaloTrack_v1& track);

  //~SiliconCaloTrack_v1() override;

  // The "standard PHObject response" functions...
  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { *this = SiliconCaloTrack_v1(); }
  int isValid() const override;
  PHObject* CloneMe() const override { return new SiliconCaloTrack_v1(*this); }

  // copy content from base class
  using PHObject::CopyFrom;  // avoid warning for not implemented CopyFrom methods
  void CopyFrom(const SiliconCaloTrack&) override;
  void CopyFrom(SiliconCaloTrack* source) override
  {
    CopyFrom(*source);
  }

  //
  // basic track information ---------------------------------------------------
  //

  unsigned int get_id() const override { return _track_id; }
  void set_id(unsigned int id) override { _track_id = id; }

  unsigned int get_calo_id() const override { return _calo_id; }
  void set_calo_id(unsigned int id) override { _calo_id = id; }

  //--float get_x() const override { return _states.find(0.0)->second->get_x(); }
  //--void set_x(float x) override { _states[0.0]->set_x(x); }

  //--float get_y() const override { return _states.find(0.0)->second->get_y(); }
  //--void set_y(float y) override { _states[0.0]->set_y(y); }

  //--float get_z() const override { return _states.find(0.0)->second->get_z(); }
  //--void set_z(float z) override { _states[0.0]->set_z(z); }

  //--float get_pos(unsigned int i) const override { return _states.find(0.0)->second->get_pos(i); }

  float get_pt()   const override { return _pt; }
  void  set_pt(float pt) override { _pt = pt; }
  float get_phi()    const override { return _phi; }
  void  set_phi(float phi) override { _phi = phi; }
  float get_eta()    const override { return _eta; }
  void  set_eta(float eta) override { _eta = eta; }

  float get_px() const override { return _pt * cos(_phi); }
  float get_py() const override { return _pt * sin(_phi); }
  float get_pz() const override { return _pt * sinh(_eta); }

  float get_p() const override { return sqrt(pow(get_px(), 2) + pow(get_py(), 2) + pow(get_pz(), 2)); }

  //
  // calo projection methods ---------------------------------------------------
  //
  float get_cal_dphi() const override { return _cal_dphi; }
  void set_cal_dphi(float dphi) override { _cal_dphi = dphi; }

  float get_cal_dz() const override { return _cal_dz; }
  void set_cal_dz(float dz) override { _cal_dz = dz; }


 private:
  // track information
  unsigned int _track_id = UINT_MAX;
  unsigned int _calo_id  = UINT_MAX;

  float _pt;
  float _phi;
  float _eta;

  // calorimeter matches
  float _cal_dphi;
  float _cal_dz;

  ClassDefOverride(SiliconCaloTrack_v1, 1)
};

#endif
