/*
to add:
holes
things touching at a point
rings
interscections / differences
*/

module shape() {
	cube(4);
	translate([2,2,2]) cube(4);
}

translate([0,0,0]) subdiv(0) shape();
translate([0,10,0]) subdiv(1) shape();
translate([0,20,0]) subdiv(2) shape();
translate([0,30,0]) subdiv(3) shape();
translate([0,40,0]) subdiv(4) shape();

translate([0,0,10]) subdiv("loop",0) shape();
translate([0,10,10]) subdiv(1,"loop") shape();
translate([0,20,10]) subdiv("loop",2) shape();
translate([0,30,10]) subdiv(3,"loop") shape();
translate([0,40,10]) subdiv("loop",4) shape();

translate([0,0,20]) subdiv("doosabin",0) shape();
translate([0,10,20]) subdiv("doosabin",1) shape();
translate([0,20,20]) subdiv("doosabin",2) shape();
translate([0,30,20]) subdiv(3,"doosabin") shape();
translate([0,40,20]) subdiv("doosabin",4) shape();

translate([0,0,30]) subdiv("sqrt3",0) shape();
translate([0,10,30]) subdiv("sqrt3",1) shape();
translate([0,20,30]) subdiv("sqrt3",2) shape();
translate([0,30,30]) subdiv("sqrt3",3) shape();
translate([0,40,30]) subdiv(4,"sqrt3") shape();

translate([0,0,-10]) subdiv() shape();
translate([0,10,-10]) subdiv(2) difference() { shape(); cube(4); }
translate([0,20,-10]) subdiv(2) difference() { shape(); sphere(4); }
translate([0,30,-10]) subdiv() square();
translate([0,40,-10]) subdiv() subdiv() shape();

translate([0,0,-20]) difference() { subdiv() shape(); translate([1,1,1]) cube(4); }
translate([0,10,-20]) difference() { subdiv() shape(); sphere(5); }
translate([0,20,-20]) subdiv("asjdifo","asdifo",13) shape();
translate([0,30,-20]) subdiv(-123,32,3) shape();
translate([0,40,-20]) subdiv("-123","1",1,"loop",32,3) shape();


