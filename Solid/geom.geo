Point(1) = {0,0, 0,1};
Point(2) = {0,0, 1,1};
Point(3) = {0,64,1,1};
Point(4) = {0,64,0,1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};
Transfinite Line {1,3} = 2 Using Progression 1;
Transfinite Line {4,2} = 128 Using Progression 1;
Transfinite Surface {6};
