/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *  Copyright (C) 2013 Felipe Sanches <fsanches@metamaquina.com.br>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <math.h>
#include <svgdata.h>
#include <dxfdata.h>
#include <polyset.h>
#include "dxftess.h"
#include "handle_dep.h"
#include "printutils.h"
#include <fstream>
#include <boost/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <Eigen/Core>

SVGData::SVGData(double fn, double fs, double fa, std::string filename, std::string layername) : fn(fn), fs(fs), fa(fa), filename(filename), layername(layername) {
  handle_dep(filename); // Register ourselves as a dependency
  this->fa=0;
  dxfdata = new DxfData();
  p = new PolySet();
  rapid_polyset = new PolySet();
  grid = new Grid2d<int>(GRID_COARSE);

  try{
    this->rapid_svgfile = new rapidxml::file<>( filename.c_str() );
    std::cout << "RapidXML read n chars:" << rapid_svgfile->size() << "\n";
    rapid_rootdoc.parse<0>(rapid_svgfile->data());
  }
  catch(const std::exception& ex)
  {
    std::cout << "RapidXML SVG parse error: " << ex.what() << std::endl;
  }

}

SVGData::~SVGData(){
  free(dxfdata);
  free(grid);
  free(rapid_svgfile);
}

void SVGData::add_point(float x, float y){
  if (!(layername.empty() || layername == layer))
   return;

  ctm.applyTransform(&x, &y);
  double px = (double) x;
  double py = (double) (document_height - y);

  dxfdata->addPathPoint( this->first_point, this->last_point, px, py, *this->grid );
/*
  if (grid->has(x, y)) {
	  this_point = grid->align(px, py);
  } else {
	  this_point = grid->align(px, py) = dxfdata->points.size();
	  dxfdata->points.push_back(Vector2d(px, py));
  }
  if (first_point < 0) {
	  dxfdata->paths.push_back(DxfData::Path());
	  first_point = this_point;
  }
  if (this_point != last_point) {
	  dxfdata->paths.back().indices.push_back(this_point);
	  last_point = this_point;
  }
*/
}

void SVGData::start_path(){
  first_point = -1;
  last_point = -1;
}

void SVGData::close_path(){
  if (first_point >= 0) {
	  dxfdata->paths.back().is_closed = 1;
	  dxfdata->paths.back().indices.push_back(first_point);
  }
}

std::vector<float> SVGData::get_path_params(std::string str){
  std::cout << "get params: " << str << "\n";
  std::vector<float> values;
  std::list<std::string> strs;
  boost::split( strs, str, boost::is_any_of(", ") );
  BOOST_FOREACH( std::string s, strs )
    if (s.size())
      try { values.push_back( boost::lexical_cast<float>(s) ); }
      catch (...) { PRINT("WARNING: SVG Parse error"); };
  return values;
}

/*!
	Returns the number of subdivision of a bezier curve, given its total length and the special variables $fn and $fs
*/
static int get_fragments_from_length(double length, double fn, double fs){
	if (fn > 0.0)
		return (int)fn;
	return (int)ceil(fmax(length / fs, 5));
}

void SVGData::add_arc_points(float xc, float yc, float rx, float ry, float start, float end){
  for (int i=0; i<=fn; i++){
    float t = ((float) i) / fn;
    float angle = start + (end-start)*t;
    float x = xc + rx*cos(angle);
    float y = yc + ry*sin(angle);
    add_point(x, y);
  }
}

void SVGData::render_rect(float x, float y, float width, float height, float rx, float ry){
  std::cout << "rect x=" << x << " y=" << y << " rx=" << rx << " ry=" << ry << " width=" << width << " height=" << height << std::endl;
  start_path();
/*  add_point(x+rx,y);
  add_point(x+width-rx,y);
  //add_arc_points(x+width-rx, y+ry, rx, ry, 3*M_PI/2, 2*M_PI);
  add_point(x+width, y+ry);
  add_point(x+width, y+height-ry);
  //add_arc_points(x+width-rx, y+height-ry, rx, ry, 0, M_PI/2);
  add_point(x+width-rx, y+height);
  add_point(x+rx, y+height);
  //add_arc_points(x+rx, y+height-ry, rx, ry, M_PI/2, M_PI);
  add_point(x, y+height-ry);
  add_point(x, y+ry);
  //add_arc_points(x+rx, y+ry, rx, ry, M_PI, 3*M_PI/2);
*/
  add_point(x,y);
  add_point(x+width,y);
  add_point(x+width,y+height);
  add_point(x, y+height);
  close_path();
}

void SVGData::render_line_to(float x0, float y0, float x1, float y1){
  std::cout << "render line: x0:" << x0 << " y0:" << y0 << " x1:" << x1 << " y1:" << y1 << std::endl;
  add_point(x0,y0);
  add_point(x1,y1);
}

float SVGData::quadratic_curve_length(float x0, float y0, float x1, float y1, float x2, float y2){
  // based on math described at
  // http://segfaultlabs.com/docs/quadratic-bezier-curve-length

  float ax, ay, bx, by;
  ax = x0 - 2*x1 + x2;
  ay = y0 - 2*y1 + y2;
  bx = 2*x1 - 2*x0;
  by = 2*y1 - 2*y0;

  float A,B,C;
  A = 4*(ax*ax+ay*ay);
  B = 4*(ax*bx + ay*by);
  C = bx*bx + by*by;
  float S = sqrt(A+B+C);
  float A32 = A*sqrt(A);

  return (1/(8.0*A32)) * (4*A32*S + 2*sqrt(A)*B*(S-sqrt(C)) + (4*C*A-B*B)*log(abs((2*sqrt(A)+B/sqrt(A)+2*S)/(B/sqrt(A) +2*sqrt(C)))));
}

void SVGData::render_quadratic_curve_to(float x0, float y0, float x1, float y1, float x2, float y2){
  int segments = get_fragments_from_length(quadratic_curve_length(x0, y0, x1, y1, x2, y2), fn, fs);
  for (double i=0; i<=segments; i++){
    double t = i/segments;
    add_point(pow(1-t,2)*x0 + 2*(1-t)*x1*t + x2*pow(t,2),
              pow(1-t,2)*y0 + 2*(1-t)*y1*t + y2*pow(t,2));
  }
}

void SVGData::render_cubic_curve_to(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3){

  int segments = 5;
  if (fn > 0.0)
    segments = fn;

  //TODO: This only deals with $fn. Add some logic to support $fs
  //int segments = get_fragments_from_length(cubic_curve_length(x0, y0, x1, y1, x2, y2, x3, y3), fn, fs);

  for (double i=0; i<=segments; i++){
    double t = i/segments;
    add_point(pow(1-t,3)*x0 + 3*pow(1-t,2)*x1*t + 3*(1-t)*x2*pow(t,2) + x3*pow(t,3),
              pow(1-t,3)*y0 + 3*pow(1-t,2)*y1*t + 3*(1-t)*y2*pow(t,2) + y3*pow(t,3));
  }
}

void SVGData::render_elliptical_arc(float x0, float y0, float rx, float 
ry, float x_axis_rotation, int large_arc_flag, int sweep_flag, float x, 
float y){

  std::cout << "elliptical arc x0 " << x0 << ","
    << "y0 " << y0 << ","
    << "radiusx " << rx << ","
    << "radiusy " << ry << ","
    << "x_axis_rot " << x_axis_rotation << ","
    << "large arc flag " << large_arc_flag << ","
    << "swap flag " << sweep_flag << ","
    << "endx " << x << ","
    << "endy " << y << "]\n";
  render_rect(x0, y0, x, y,0,0);
  return;

  // Conversion from endpoint to center parameterization
  // http://www.w3.org/TR/SVG11/implnote.html#ArcImplementationNotes
  float cx, cy, cx_, cy_, x0_, y0_, start_angle, end_angle;

  x0_ = cos(x_axis_rotation)*(x0-x)/2 - sin(x_axis_rotation)*(y0-y)/2;
  y0_ = -sin(x_axis_rotation)*(x0-x)/2 + cos(x_axis_rotation)*(y0-y)/2;

  float k = sqrt((rx*rx*ry*ry - rx*rx*y0_*y0_ - ry*ry*x0_*x0_) / (rx*rx*y0_*y0_ + ry*ry*x0_*x0_));
  cx_ = k * rx*y0_/ry;
  cy_ = -k * ry*x0_/rx;

  if (large_arc_flag==sweep_flag){
    cx_ = -cx_;
    cy_ = -cy_;
  }

  cx = cos(x_axis_rotation) * cx_ - sin(x_axis_rotation) * cy_ + (x0+x)/2;
  cy = sin(x_axis_rotation) * cx_ + cos(x_axis_rotation) * cy_ + (y0+y)/2;

  //TODO: calculate start_angle and end_angle

  //TODO: This only deals with $fn. Add some logic to support $fa and $fs
  for (int i=0; i<=fn; i++){
    double t = i/fn;
    //TODO: generate points for the elliptical curve by calling add_point(x, y)
  }
}

void SVGData::parse_path_description(std::string d){
  std::string pattern = "(([mMzZlLhHvVcCsSqQtTaA])([^mMzZlLhHvVcCsSqQtTaA]*))";
  boost::regex regexPattern(pattern);
  boost::match_results<std::string::const_iterator> result;

  std::string::const_iterator start, end;
  start = d.begin();
  end = d.end();

  float x=0, y=0;
  while(boost::regex_search(start, end, result, regexPattern)){
    std::cout << "result: " << std::endl << result[1] << std::endl << std::endl;

    char instruction_code = ((std::string) result[2]).c_str()[0];
    std::vector<float> params;
    params = get_path_params(((std::string) result[3]).c_str());
    std::cout << "instruction_code=" << instruction_code << std::endl;

    int idx=0;
    switch (instruction_code){
      case 'm':
        std::cout << "moveto command - relative\n";

        //Since we're doing a move, we must close the path we may have drawn so far
        close_path();

        //And prepare to start drawing something else from our new starting point after the move command
        start_path();

        while (params.size() - idx >= 2){
          if (idx>0)
            render_line_to(x, y, x+params[idx], y+params[idx+1]);
          x += params[idx];
          y += params[idx+1];
          std::cout << "m: x=" << x << " y=" << y << std::endl;
          idx+=2;
        }
        break;
      case 'M':
        std::cout << "moveto command - absolute\n";

        //Since we're doing a move, we must close the path we may have drawn so far
        close_path();

        //And prepare to start drawing something else from our new starting point after the move command
        start_path();

        while (params.size() - idx >= 2){
          if (idx>0)
            render_line_to(x, y, params[idx], params[idx+1]);
          x = params[idx];
          y = params[idx+1];
          std::cout << "M: x=" << x << " y=" << y << std::endl;
          idx+=2;
        }
        break;
      case 'z':
      case 'Z':
        std::cout << "closepath command\n";
        //do nothing?
        break;
      case 'l':
        std::cout << "lineto command - relative\n";
        while (params.size() - idx >= 2){
          render_line_to(x, y, x+params[idx], y+params[idx+1]);
          x += params[idx];
          y += params[idx+1];
          std::cout << "l: x=" << x << " y=" << y << std::endl;
          idx+=2;
        }
        break;
      case 'L':
        std::cout << "lineto command - absolute\n";
        while (params.size() - idx >= 2){
          render_line_to(x, y, params[idx], params[idx+1]);
          x = params[idx];
          y = params[idx+1];
          std::cout << "L: x=" << x << " y=" << y << std::endl;
          idx+=2;
        }
        break;
      case 'h':
        std::cout << "horizontal lineto command - relative\n";
        while (params.size() - idx >= 1){
          render_line_to(x, y, x+params[idx], y);
          x += params[idx];
          std::cout << "h: x=" << x << " y=" << y << std::endl;
          idx++;
        }
        break;
      case 'H':
        std::cout << "horizontal lineto command - absolute\n";
        while (params.size() - idx >= 1){
          render_line_to(x, y, params[idx], y);
          x = params[idx];
          std::cout << "H: x=" << x << " y=" << y << std::endl;
          idx++;
        }
        break;
      case 'v':
        std::cout << "vertical lineto command - relative\n";
        while (params.size() - idx >= 1){
          render_line_to(x, y, x, y+params[idx]);
          y += params[idx];
          std::cout << "v: x=" << x << " y=" << y << std::endl;
          idx++;
        }
        break;
      case 'V':
        std::cout << "vertical lineto command - absolute\n";
        while (params.size() - idx >= 1){
          render_line_to(x, y, x, params[idx]);
          y = params[idx];
          std::cout << "V: x=" << x << " y=" << y << std::endl;
          idx++;
        }
        break;
      case 'c':
        std::cout << "curveto cubic Bézier command - relative\n";
        while (params.size() - idx >= 6){
          render_cubic_curve_to(x, y, x+params[idx], y+params[idx+1], x+params[idx+2], y+params[idx+3], x+params[idx+4], y+params[idx+5]);
          x += params[idx+4];
          y += params[idx+5];
          std::cout << "c: x=" << x << " y=" << y << std::endl;
          idx+=6;
        }
        break;
      case 'C':
        std::cout << "curveto cubic Bézier command - absolute\n";
        while (params.size() - idx >= 6){
          render_cubic_curve_to(x, y, params[idx], params[idx+1], params[idx+2], params[idx+3], params[idx+4], params[idx+5]);
          x = params[idx+4];
          y = params[idx+5];
          std::cout << "C: x=" << x << " y=" << y << std::endl;
          idx+=6;
        }
        break;
      case 's':
        std::cout << "shorthand/smooth curveto cubic Bézier command - relative\n";
        while (params.size() - idx >= 4){
          std::cout << "this is wrong! TODO: implement reflection of previous control point as described in the SVG spec.\n";
          render_cubic_curve_to(x, y, x, y, x+params[idx], y+params[idx+1], x+params[idx+2], y+params[idx+3]);
          x += params[idx+2];
          y += params[idx+3];
          std::cout << "s: x=" << x << " y=" << y << std::endl;
          idx+=4;
        }
        break;
      case 'S':
        std::cout << "shorthand/smooth curveto cubic Bézier command - absolute\n";
        while (params.size() - idx >= 4){
          std::cout << "this is wrong! TODO: implement reflection of previous control point as described in the SVG spec.\n";
          render_cubic_curve_to(x, y, x, y, params[idx], params[idx+1], params[idx+2], params[idx+3]);
          x = params[idx+2];
          y = params[idx+3];
          std::cout << "S: x=" << x << " y=" << y << std::endl;
          idx+=4;
        }
        break;
      case 'q':
        std::cout << "curveto quadratic Bézier command - relative\n";
        while (params.size() - idx >= 4){
          render_quadratic_curve_to(x, y, x+params[idx], y+params[idx+1], x+params[idx+2], y+params[idx+3]);
          x += params[idx+2];
          y += params[idx+3];
          std::cout << "q: x=" << x << " y=" << y << std::endl;
          idx+=4;
        }
        break;
      case 'Q':
        std::cout << "curveto quadratic Bézier command - absolute\n";
        while (params.size() - idx >= 4){
          render_quadratic_curve_to(x, y, params[idx], params[idx+1], params[idx+2], params[idx+3]);
          x = params[idx+2];
          y = params[idx+3];
          std::cout << "Q: x=" << x << " y=" << y << std::endl;
          idx+=4;
        }
        break;
      case 't':
        std::cout << "shorthand/smooth curveto quadratic Bézier command - relative\n";
        while (params.size() - idx >= 2){
          //this is wrong! TODO: implement reflection of previous control point as described in the SVG spec.
          render_quadratic_curve_to(x, y, x, y, x+params[idx], y+params[idx+1]);
          x += params[idx];
          y += params[idx+1];
          //std::cout << "t: x=" << x << " y=" << y << std::endl;
          idx+=2;
        }
        break;
      case 'T':
        std::cout << "shorthand/smooth curveto quadratic Bézier command - absolute\n";
        while (params.size() - idx >= 2){
          //this is wrong! TODO: implement reflection of previous control point as described in the SVG spec.
          render_quadratic_curve_to(x, y, x, y, params[idx], params[idx+1]);
          x = params[idx];
          y = params[idx+1];
          //std::cout << "T: x=" << x << " y=" << y << std::endl;
          idx+=2;
        }
        break;
      case 'a':
        std::cout << "elliptical arc command - relative\n";
        while (params.size() - idx >= 7){
          render_elliptical_arc(x, y, params[idx], params[idx+1], params[idx+2], params[idx+3], params[idx+4], x+params[idx+5], y+params[idx+6]);
          x += params[idx+5];
          y += params[idx+6];
          //std::cout << "a: x=" << x << " y=" << y << std::endl;
          idx+=7;
        }
        break;
      case 'A':
        std::cout << "elliptical arc command - absolute\n";
        while (params.size() - idx >= 7){
          render_elliptical_arc(x, y, params[idx], params[idx+1], params[idx+2], params[idx+3], params[idx+4], params[idx+5], params[idx+6]);
          x = params[idx+5];
          y = params[idx+6];
          //std::cout << "A: x=" << x << " y=" << y << std::endl;
          idx+=7;
        }
        break;
    }

    start = result[1].second;
  }
}

#define NUMBER_REGEX "-?[0-9]+(\\.[0-9]+)?(e-?[0-9]+)?"
TransformMatrix SVGData::parse_transform(std::string transform){
  std::cout << "Parsing SVG 'transform' atrtibute = '" << transform << "'" << std::endl;

  TransformMatrix tm;
  tm.setIdentity();

  boost::regex matrix_regex("\\s*matrix\\(\\s*(" NUMBER_REGEX ")\\s*,\\s*(" NUMBER_REGEX ")\\s*,\\s*(" NUMBER_REGEX ")\\s*,\\s*(" NUMBER_REGEX ")\\s*,\\s*(" NUMBER_REGEX ")\\s*,\\s*(" NUMBER_REGEX ")\\s*\\)");

  boost::regex translate_regex("\\s*translate\\(\\s*(" NUMBER_REGEX ")\\s*,\\s*(" NUMBER_REGEX ")?\\s*\\)");

  boost::regex scale_regex("\\s*scale\\(\\s*(" NUMBER_REGEX ")\\s*,\\s*(" NUMBER_REGEX ")?\\s*\\)");

  boost::regex rotate_regex("\\s*rotate\\(\\s*(" NUMBER_REGEX ")\\s*,\\s*((" NUMBER_REGEX ")\\s*,\\s*(" NUMBER_REGEX "))?\\s*\\)");

  boost::regex skewX_regex("\\s*skewX\\(\\s*(" NUMBER_REGEX ")\\s*\\)");

  boost::regex skewY_regex("\\s*skewY\\(\\s*(" NUMBER_REGEX ")\\s*\\)");

  boost::match_results<std::string::const_iterator> result;

  std::string::const_iterator start, end;
  start = transform.begin();
  end = transform.end();

  bool continue_parsing = true;
  while(continue_parsing){
    continue_parsing = false;

    if (boost::regex_search(start, end, result, matrix_regex)){
      float a = atof(((std::string) result[1]).c_str());
      float b = atof(((std::string) result[4]).c_str());
      float c = atof(((std::string) result[7]).c_str());
      float d = atof(((std::string) result[10]).c_str());
      float e = atof(((std::string) result[13]).c_str());
      float f = atof(((std::string) result[16]).c_str());
      TransformMatrix transform_matrix(a, b, c, d, e, f);
      //std::cout << "found a matrix transform!" << std::endl;
      tm = tm * transform_matrix;
      start = result[0].second;
      continue_parsing = true;
    }

    if (boost::regex_search(start, end, result, translate_regex)){
      float tx = atof(((std::string) result[1]).c_str());
      float ty = atof(((std::string) result[4]).c_str());
      //std::cout << "found a translate transform!" << std::endl;
      tm.translate(tx, ty);
      start = result[0].second;
      continue_parsing = true;
    }

    if (boost::regex_search(start, end, result, scale_regex)){
      float sx = atof(((std::string) result[1]).c_str());
      float sy = atof(((std::string) result[4]).c_str());
      //std::cout << "found a scale transform!" << std::endl;
      tm.scale(sx, sy);
      start = result[0].second;
      continue_parsing = true;
    }

    if (boost::regex_search(start, end, result, rotate_regex)){
      float angle = atof(((std::string) result[1]).c_str());
      float cx = atof(((std::string) result[5]).c_str());
      float cy = atof(((std::string) result[8]).c_str());
      //std::cout << "found a rotate transform!" << std::endl;
      tm.rotate(angle, cx, cy);
      start = result[0].second;
      continue_parsing = true;
    }

    if (boost::regex_search(start, end, result, skewX_regex)){
      float angle = atof(((std::string) result[1]).c_str());
      //std::cout << "found a skewX transform!" << std::endl;
      tm.skewX(angle);
      start = result[0].second;
      continue_parsing = true;
    }

    if (boost::regex_search(start, end, result, skewY_regex)){
      float angle = atof(((std::string) result[1]).c_str());
      //std::cout << "found a skewY transform!" << std::endl;
      tm.skewY(angle);
      start = result[0].second;
      continue_parsing = true;
    }

  }

  return tm;
}

void SVGData::rapid_traverse_subtree(TransformMatrix parent_matrix, rapidxml::xml_node<> * node)
{
  TransformMatrix tm = parent_matrix;

  std::string nodename;
  if ( node ) nodename = node->name();

  if ( node->type() == rapidxml::node_element ) {
    rapidxml::xml_attribute<> *inkscape_label = NULL;
    rapidxml::xml_attribute<> *inkscape_groupmode = NULL;
    rapidxml::xml_attribute<> *transform = NULL;
    inkscape_label = node->first_attribute( "inkscape:label" );
    inkscape_groupmode = node->first_attribute( "inkscape:groupmode" );
    transform = node->first_attribute( "transform" );

    if (transform)
      tm = parent_matrix * parse_transform(transform->value());

    if (nodename == "g")
      if (inkscape_label && inkscape_groupmode)
        if (!strcmp(inkscape_groupmode->value(),"layer"))
          this->layer = std::string(inkscape_label->value());
  }

  setCurrentTransformMatrix(tm);

  if ( nodename == "svg" ){
    if( node ) {
      rapidxml::xml_attribute<> *height = node->first_attribute( "height" );
      if (height) document_height = atof(height->value());
    }
  }

  if ( nodename == "path" ){
    if( node ) {
      rapidxml::xml_attribute<> *d = node->first_attribute( "d" );
      if (d) {
        std::cout << "path description = " << d->value() << std::endl;
        start_path();
        parse_path_description(d->value());
        close_path();
      }
    }
  }

  if (nodename == "rect") {
    if (node) {
      rapidxml::xml_attribute<>* width_attr = node->first_attribute("width");
      rapidxml::xml_attribute<>* height_attr = node->first_attribute("height");
      rapidxml::xml_attribute<>* x_attr = node->first_attribute("x");
      rapidxml::xml_attribute<>* y_attr = node->first_attribute("y");
      rapidxml::xml_attribute<>* rx_attr = node->first_attribute("rx");
      rapidxml::xml_attribute<>* ry_attr = node->first_attribute("ry");

      float x, y, width, height, rx, ry;

      x = x_attr ? atof(x_attr->value()) : 0;
      y = y_attr ? atof(y_attr->value()) : 0;

      if(!rx_attr && !ry_attr){
        rx=0; ry=0;
      }

      if(rx_attr && !ry_attr){
        rx=ry=atof(rx_attr->value());
      }

      if(!rx_attr && ry_attr){
        rx=ry=atof(ry_attr->value());
      }

      if(rx_attr && ry_attr){
        rx=atof(rx_attr->value());
        ry=atof(ry_attr->value());
      }

      if(width_attr && height_attr){
        width = atof(width_attr->value());
        height = atof(height_attr->value());

        if (rx > width/2)
          rx = width/2;

        if (ry > height/2)
          ry = height/2;

        if(width > 0 && height > 0){
          render_rect(x, y, width, height, rx, ry);
          std::cout << "[svg:rect] x:" << x << " y:" << y << " width:" << width << " height:" << height << " rx:" << rx << " ry:" << ry << std::endl;
        }
      }
    }
  }

  if (node) {
    rapidxml::xml_node<> *i=NULL;
    for ( i = node->first_node(); i!=0 ; i=i->next_sibling() ) {
      rapid_traverse_subtree( tm, i );
    }
  }

  layer.empty();
}

PolySet* SVGData::convertToPolyset(){
  TransformMatrix tm;
  tm.setIdentity();

  rapid_traverse_subtree(tm, rapid_rootdoc.first_node() );
  if (rapid_polyset) {
    rapid_polyset->is2d = true;
    dxf_tesselate(rapid_polyset, *dxfdata, 0, Vector2d(1,1), true, false, 0);
    dxf_border_to_ps(rapid_polyset, *dxfdata);
  }

  std::cout << rapid_polyset->dump();

  return rapid_polyset;
}
