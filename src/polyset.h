#ifndef POLYSET_H_
#define POLYSET_H_

#include "Geometry.h"
#include "system-gl.h"
#include "grid.h"
#include "linalg.h"
#include  "renderer.h"
#include <vector>
#include <string>

class PolySet : public Geometry
{
public:
	typedef std::vector<Vector3d> Polygon;
	std::vector<Polygon> polygons;
	std::vector<Polygon> borders;
	Grid3d<void*> grid;

	bool is2d;

	PolySet();
	virtual ~PolySet();

	virtual size_t memsize() const;
	virtual BoundingBox getBoundingBox() const;
	virtual std::string dump() const;
	virtual unsigned int getDimension() const { return this->is2d ? 2 : 3; }

	bool empty() const { return polygons.size() == 0; }
	void append_poly();
	void append_vertex(double x, double y, double z = 0.0);
	void append_vertex(Vector3d v);
	void insert_vertex(double x, double y, double z = 0.0);
	void insert_vertex(Vector3d v);
	void append(const PolySet &ps);

	void render_surface(Renderer::csgmode_e csgmode, const Transform3d &m, GLint *shaderinfo = NULL) const;
	void render_edges(Renderer::csgmode_e csgmode) const;
};

#endif
