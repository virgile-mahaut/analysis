#ifndef SILICONCALOTRACKMAP_H
#define SILICONCALOTRACKMAP_H

#include "SiliconCaloTrack.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>

class SiliconCaloTrackMap : public PHObject
{
 public:
  typedef std::map<unsigned int, SiliconCaloTrack*> TrackMap;
  typedef std::map<unsigned int, SiliconCaloTrack*>::const_iterator ConstIter;
  typedef std::map<unsigned int, SiliconCaloTrack*>::iterator Iter;

  ~SiliconCaloTrackMap() override {}

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SiliconCaloTrackMap base class" << std::endl;
  }
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int /*idkey*/) const { return 0; }
  virtual void clear() {}

  virtual const SiliconCaloTrack* get(unsigned int /*idkey*/) const { return nullptr; }
  virtual SiliconCaloTrack* get(unsigned int /*idkey*/) { return nullptr; }
  virtual SiliconCaloTrack* insert(const SiliconCaloTrack* /*cluster*/) { return nullptr; }
  virtual SiliconCaloTrack* insertWithKey(const SiliconCaloTrack*, unsigned int) { return nullptr; }
  virtual size_t erase(unsigned int /*idkey*/) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  SiliconCaloTrackMap() {}

 private:
  ClassDefOverride(SiliconCaloTrackMap, 1);
};

#endif
