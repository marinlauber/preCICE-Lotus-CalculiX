Point(1) = { 2,0,-0.5,1};
Point(2) = { 2,0, 0.5,1};
Point(3) = { 3, 0, 0.5,1};
Point(4) = { 3, 0,-0.5,1};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};
Transfinite Line {1,3} = 2 Using Progression 1;
Transfinite Line {4,2} = 64 Using Progression 1;
Transfinite Surface {6};
