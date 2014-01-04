#ifndef POLYGON2D_H_
#define POLYGON2D_H_

#include "Geometry.h"
#include "linalg.h"
#include <vector>

struct Outline2d {
	Outline2d() : positive(true) {}
	std::vector<Vector2d> vertices;
	bool positive;
};

class Polygon2d : public Geometry
{
public:
	Polygon2d() : sanitized(false) {}
	virtual size_t memsize() const;
	virtual BoundingBox getBoundingBox() const;
	virtual std::string dump() const;
	virtual unsigned int getDimension() const { return 2; }
	virtual bool isEmpty() const;

	void addOutline(const Outline2d &outline) { this->theoutlines.push_back(outline); }
	class PolySet *tessellate() const;

	typedef std::vector<Outline2d> Outlines2d;
	const Outlines2d &outlines() const { return theoutlines; }

	void transform(const Transform2d &mat);
	void resize(Vector2d newsize, const Eigen::Matrix<bool,2,1> &autosize);

	bool isSanitized() const { return this->sanitized; }
	void setSanitized(bool s) { this->sanitized = s; }
private:
	Outlines2d theoutlines;
	bool sanitized;
};

#endif
