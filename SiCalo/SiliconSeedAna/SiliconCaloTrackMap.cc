#include "SiliconCaloTrackMap.h"

namespace
{
  SiliconCaloTrackMap::TrackMap DummyTrackMap;
}

SiliconCaloTrackMap::ConstIter SiliconCaloTrackMap::begin() const
{
  return DummyTrackMap.end();
}

SiliconCaloTrackMap::ConstIter SiliconCaloTrackMap::find(unsigned int /*idkey*/) const
{
  return DummyTrackMap.end();
}

SiliconCaloTrackMap::ConstIter SiliconCaloTrackMap::end() const
{
  return DummyTrackMap.end();
}

SiliconCaloTrackMap::Iter SiliconCaloTrackMap::begin()
{
  return DummyTrackMap.end();
}

SiliconCaloTrackMap::Iter SiliconCaloTrackMap::find(unsigned int /*idkey*/)
{
  return DummyTrackMap.end();
}

SiliconCaloTrackMap::Iter SiliconCaloTrackMap::end()
{
  return DummyTrackMap.end();
}
