#include "SiliconCaloTrackMap_v1.h"

#include "SiliconCaloTrack.h"

#include <phool/PHObject.h>  // for PHObject

#include <iterator>  // for reverse_iterator
#include <map>       // for _Rb_tree_const_iterator, _Rb_tree_iterator
#include <ostream>   // for operator<<, endl, ostream, basic_ostream, bas...
#include <utility>   // for pair, make_pair

SiliconCaloTrackMap_v1::SiliconCaloTrackMap_v1()
  : _map()
{
}

SiliconCaloTrackMap_v1::SiliconCaloTrackMap_v1(const SiliconCaloTrackMap_v1& trackmap)
  : _map()
{
  for (auto iter : trackmap)
  {
    auto track = static_cast<SiliconCaloTrack*>(iter.second->CloneMe());
    _map.insert(std::make_pair(track->get_id(), track));
  }
}

SiliconCaloTrackMap_v1& SiliconCaloTrackMap_v1::operator=(const SiliconCaloTrackMap_v1& trackmap)
{
  // do nothing if same  copying map onto itself
  if (&trackmap == this)
  {
    return *this;
  }

  Reset();
  for (auto iter : trackmap)
  {
    auto track = static_cast<SiliconCaloTrack*>(iter.second->CloneMe());
    _map.insert(std::make_pair(track->get_id(), track));
  }
  return *this;
}

SiliconCaloTrackMap_v1::~SiliconCaloTrackMap_v1()
{
  Reset();
}

void SiliconCaloTrackMap_v1::Reset()
{
  for (auto& iter : _map)
  {
    SiliconCaloTrack* track = iter.second;
    delete track;
  }
  _map.clear();
}

void SiliconCaloTrackMap_v1::identify(std::ostream& os) const
{
  os << "SiliconCaloTrackMap_v1: size = " << _map.size() << std::endl;
  return;
}

const SiliconCaloTrack* SiliconCaloTrackMap_v1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

SiliconCaloTrack* SiliconCaloTrackMap_v1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

SiliconCaloTrack* SiliconCaloTrackMap_v1::insert(const SiliconCaloTrack* track)
{
  unsigned int index = 0;
  if (!_map.empty())
  {
    index = _map.rbegin()->first + 1;
  }
  auto copy = static_cast<SiliconCaloTrack*>(track->CloneMe());
  copy->set_id(index);

  const auto result = _map.insert(std::make_pair(index, copy));
  if (!result.second)
  {
    std::cout << "SiliconCaloTrackMap_v1::insert - duplicated key. track not inserted" << std::endl;
    delete copy;
    return nullptr;
  }
  else
  {
    return copy;
  }
}

SiliconCaloTrack* SiliconCaloTrackMap_v1::insertWithKey(const SiliconCaloTrack* track, unsigned int index)
{
  auto copy = static_cast<SiliconCaloTrack*>(track->CloneMe());
  copy->set_id(index);
  const auto result = _map.insert(std::make_pair(index, copy));
  if (!result.second)
  {
    std::cout << "SiliconCaloTrackMap_v1::insertWithKey - duplicated key. track not inserted" << std::endl;
    delete copy;
    return nullptr;
  }
  else
  {
    return copy;
  }
}
