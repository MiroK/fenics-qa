Point(1) = {0, 0, 0, 1};
Point(2) = {10, 0, 0, 1};
Point(3) = {20, 0, 0, 1};
Point(4) = {0, 10, 0, 1};
Point(5) = {0, 20, 0, 1};
Line(1) = {2, 3};
Line(2) = {4, 5};
Circle(3) = {4, 1, 2};
Circle(4) = {5, 1, 3};
Line Loop(6) = {2, 4, -1, -3};
Plane Surface(7) = {6};
Physical Surface(1) = {7};
Physical Line(1) = {2};
Physical Line(2) = {1};
Physical Line(3) = {3};
Physical Line(4) = {4};
