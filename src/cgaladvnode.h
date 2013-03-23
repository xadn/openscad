#ifndef CGALADVNODE_H_
#define CGALADVNODE_H_

#include "node.h"
#include "visitor.h"
#include "value.h"
#include "linalg.h"
#include <set>
#include <string>

enum cgaladv_type_e {
	MINKOWSKI,
	GLIDE,
	SUBDIV,
	HULL,
	RESIZE
};

enum cgaladv_subdiv_type_e {
	SUBDIV_CATMULL_CLARK,
	SUBDIV_LOOP,
	SUBDIV_DOO_SABIN,
	SUBDIV_SQRT3
};

class CgaladvNode : public AbstractNode
{
public:
	CgaladvNode(const ModuleInstantiation *mi, cgaladv_type_e type) : AbstractNode(mi), type(type) {
		convexity = 1;
	}
	virtual ~CgaladvNode() { }
  virtual Response accept(class State &state, Visitor &visitor) const {
		return visitor.visit(state, *this);
	}
	virtual std::string toString() const;
	virtual std::string name() const;
	PolySet *evaluate_polyset(class PolySetEvaluator *ps) const;

	Value path;
	cgaladv_subdiv_type_e subdiv_type;
	int convexity, subdiv_level;
	Vector3d newsize;
	Eigen::Matrix<bool,3,1> autosize;
	cgaladv_type_e type;
};

#endif
