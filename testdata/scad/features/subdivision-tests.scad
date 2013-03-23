/*
to add:
holes
things touching at a point
rings
interscections / differences
*/

module basica() {
	cube(4);
	translate([2,2,2]) cube(4);
}

translate([0,0,0]) subdiv(0) basica();
translate([0,10,0]) subdiv(1) basica();
translate([0,20,0]) subdiv(2) basica();
translate([0,30,0]) subdiv(3) basica();
translate([0,40,0]) subdiv(4) basica();

translate([0,0,10]) subdiv("loop",0) basica();
translate([0,10,10]) subdiv("loop",1) basica();
translate([0,20,10]) subdiv("loop",2) basica();
translate([0,30,10]) subdiv("loop",3) basica();
translate([0,40,10]) subdiv("loop",4) basica();

translate([0,0,20]) subdiv("doosabin",0) basica();
translate([0,10,20]) subdiv("doosabin",1) basica();
translate([0,20,20]) subdiv("doosabin",2) basica();
translate([0,30,20]) subdiv("doosabin",3) basica();
translate([0,40,20]) subdiv("doosabin",4) basica();

translate([0,0,30]) subdiv("sqrt3",0) basica();
translate([0,10,30]) subdiv("sqrt3",1) basica();
translate([0,20,30]) subdiv("sqrt3",2) basica();
translate([0,30,30]) subdiv("sqrt3",3) basica();
translate([0,40,30]) subdiv("sqrt3",4) basica();

translate([0,0,-10]) subdiv() basica();
