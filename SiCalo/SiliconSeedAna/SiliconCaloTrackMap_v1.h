#ifndef SILICONCALOTRACKMAPV1_H
#define SILICONCALOTRACKMAPV1_H

#include "SiliconCaloTrack.h"
#include "SiliconCaloTrackMap.h"

#include <cstddef>   // for size_t
#include <iostream>  // for cout, ostream

class PHObject;

class SiliconCaloTrackMap_v1 : public SiliconCaloTrackMap
{
 public:
  SiliconCaloTrackMap_v1();
  SiliconCaloTrackMap_v1(const SiliconCaloTrackMap_v1& trackmap);
  SiliconCaloTrackMap_v1& operator=(const SiliconCaloTrackMap_v1& trackmap);
  ~SiliconCaloTrackMap_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  // cppcheck-suppress virtualCallInConstructor
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new SiliconCaloTrackMap_v1(*this); }

  bool empty() const override { return _map.empty(); }
  size_t size() const override { return _map.size(); }
  size_t count(unsigned int idkey) const override { return _map.count(idkey); }
  void clear() override { Reset(); }

  const SiliconCaloTrack* get(unsigned int idkey) const override;
  SiliconCaloTrack* get(unsigned int idkey) override;
  SiliconCaloTrack* insert(const SiliconCaloTrack* track) override;
  SiliconCaloTrack* insertWithKey(const SiliconCaloTrack* track, unsigned int index) override;
  size_t erase(unsigned int idkey) override
  {
    delete _map[idkey];
    return _map.erase(idkey);
  }

  ConstIter begin() const override { return _map.begin(); }
  ConstIter find(unsigned int idkey) const override { return _map.find(idkey); }
  ConstIter end() const override { return _map.end(); }

  Iter begin() override { return _map.begin(); }
  Iter find(unsigned int idkey) override { return _map.find(idkey); }
  Iter end() override { return _map.end(); }

 private:
  TrackMap _map;

  ClassDefOverride(SiliconCaloTrackMap_v1, 1);
};

#endif
