#include "GeometryCache.h"
#include "printutils.h"
#include "geometry.h"
#ifdef DEBUG
  #include "CGAL_Nef_polyhedron.h"
#endif

GeometryCache *GeometryCache::inst = NULL;

bool GeometryCache::insert(const std::string &id, const shared_ptr<const Geometry> &geom)
{
	bool inserted = this->cache.insert(id, new cache_entry(geom), geom ? geom->memsize() : 0);
#ifdef DEBUG
	assert(!dynamic_cast<const CGAL_Nef_polyhedron*>(geom.get()));
	if (inserted) PRINTB("Geometry Cache insert: %s (%d bytes)", 
                         id.substr(0, 40) % (geom ? geom->memsize() : 0));
	else PRINTB("Geometry Cache insert failed: %s (%d bytes)",
                id.substr(0, 40) % (geom ? geom->memsize() : 0));
#endif
	return inserted;
}

size_t GeometryCache::maxSize() const
{
	return this->cache.maxCost();
}

void GeometryCache::setMaxSize(size_t limit)
{
	this->cache.setMaxCost(limit);
}

void GeometryCache::print()
{
	PRINTB("Geometries in cache: %d", this->cache.size());
	PRINTB("Geometry cache size in bytes: %d", this->cache.totalCost());
}

GeometryCache::cache_entry::cache_entry(const shared_ptr<const Geometry> &geom)
	: geom(geom)
{
	if (print_messages_stack.size() > 0) this->msg = print_messages_stack.back();
}
