function [theta_graph, theta_anal] = generalized_subproblem3(xi1,xi2,p,q1,q2,delta1,delta2,point2figure)
% Solves Generalized Subproblem 3 based on paper:
% "Solution for a new Subproblem in Screw Theory and its application", by
% Tan Chen Xiao
% Same as ex.6 p.147 Murray

% load Murray kinematics
addpath('/home/nikos/matlab_ws/kinematics/robotlinks')
addpath('/home/nikos/matlab_ws/kinematics/screws') 
addpath('/home/nikos/matlab_ws/kinematics/util')
% load Geometry 2-3D libraries
addpath('/home/nikos/matlab_ws/geom3d/geom3d')
addpath('/home/nikos/matlab_ws/geom2d/geom2d')
addpath('/home/nikos/matlab_ws/geom2d/utils')

%% Step 1: Determine the closest points r1,r2 between the axes of twists ξ1,ξ2
% First we create the lines in 3Dspace given the twists ξ1,ξ2
figure(point2figure)
axis_xi1 = drawtwist(xi1); hold on
axis_xi2 = drawtwist(xi2); hold on
% p1_1 - axis_xi1 - p1_2
% p2_1 - axis_xi2 - p2_2
p1_1 = [axis_xi1.XData(1) axis_xi1.YData(1) axis_xi1.ZData(1) ]';
p1_2 = [axis_xi1.XData(2) axis_xi1.YData(2) axis_xi1.ZData(2) ]';
text(p1_2(1), p1_2(2), p1_2(3),'\leftarrow ξ1'), hold on

p2_1 = [axis_xi2.XData(1) axis_xi2.YData(1) axis_xi2.ZData(1) ]';
p2_2 = [axis_xi2.XData(2) axis_xi2.YData(2) axis_xi2.ZData(2) ]';
text(p2_2(1), p2_2(2), p2_2(3),'\leftarrow ξ2'), hold on

% We define , the distance, the connecting vector and the point in each twist
% axis
[distance, d, r1, r2] = DistBetween2Segment(p1_1, p1_2, p2_1, p2_2);
figure(point2figure)

plot3([r1(1) r2(1)],[r1(2) r2(2)],[r1(3) r2(3)]), hold on

figure(point2figure)
scatter3(r1(1), r1(2), r1(3),'b'), hold on
text(r1(1), r1(2), r1(3),'\leftarrow r1'), hold on
scatter3(r2(1), r2(2), r2(3),'b'), hold on
text(r2(1), r2(2), r2(3),'\leftarrow r2'), hold on
scatter3(q1(1), q1(2), q1(3),'*m'), hold on
text(q1(1), q1(2), q1(3),'\leftarrow q1'), hold on
scatter3(q2(1), q2(2), q2(3),'*m'), hold on
text(q2(1), q2(2), q2(3),'\leftarrow q2'), hold on
scatter3(p(1), p(2), p(3),'*m'), hold on
text(p(1), p(2), p(3),'\leftarrow p'), hold on

%% Step 2: Determine the plane that is normal to the one unit vector and pass through p
% Must be Σ2: w2 the normal vector to plane surface, p ε plane

% Based on figure 1 of page 4
w1=xi1(4:6,:);
w2=xi2(4:6,:);

% First find second point of plane, that is to determine projections of p to ξ2

u = p(1:3)-r2;
s2 = r2 + w2*(w2'*u); % projection of p to ξ2 and Σ2

% Plane evaluation, p,s2 must lie on the same plane
eval_Sigma2 = dot((p(1:3)-s2),w2);

figure(point2figure)
Sigma2 = createPlane(p(1:3)', w2');
drawPlane3d(Sigma2,'facecolor', 'y','facealpha', '0.5'), hold on
% [Sigma2,a2,b2] = plot_plane(p(1:3), w2, point2figure), hold on


%% Step 3: Specify circle Q in 3D space -  Intersection of spheres
% Point q lies in the intersection of three spheres

% Sphere 1
Sph1 = [q1(1) q1(2) q1(3) delta1];
figure(point2figure)
drawSphere(Sph1,'facecolor', 'r','facealpha', '0.5'), hold on
% Sphere 2
Sph2 = [q2(1) q2(2) q2(3) delta2];
figure(point2figure)
drawSphere(Sph2,'facecolor', 'm','facealpha', '0.5'), hold on

% Circle Q is the intersection of Spheres 1-2
% We find the line connecting the centers of the 2 spheres
q1q2 = createLine3d(q1(1:3)', q2(1:3)');
unit_q1q2 = q1q2(4:6)/sqrt(q1q2(4:6)*q1q2(4:6)'); % normalize vector
% figure(point2figure)
% drawLine3d(q1q2), hold on

% The intersection of two spheres is a cicle in the plane with normal
% vector the direction vector of this line, and the center is positioned in
% this line

% Since in THIS problem the two spheres only have one touching point,
% circle becomes a point, so point q lies in the line found above

% So vectors (q1q2) and (q1q) are collinear :
% cross(q1q2(1,4:6),q1q(1,4:6))=0 eq.(s1)
% dot(unit_q1q2,unit_q1q)=1 eq.(s2)

% Symbolic Declaration of point q
syms q_x real;
syms q_y real;
syms q_z real;
q = [q_x q_y q_z]';
q1q = createLine3d(q1(1:3)', q(1:3)');
unit_q1q = q1q(4:6)/sqrt(q1q(4:6)*q1q(4:6)'); % normalize vector
% The following 3 equations are extracted from symbolic calculation of
% eq.(s1) and eq.(s2)

% eq.1: -2*2^(1/2)*q_y = 0
% eq.2: q_x - q_z = -0.1610 from eq.(s1)
% eq.-: 2*2^(1/2)*q_y = 0 same with eq.1 !!!

% eq.3: ... from eq.(s2)
%% System of equations for q point calculation
lamda = delta1*(delta1+delta2);

first = sym(q_y);
second = sym(cross(q1q2(1,4:6),q1q(1,4:6))); % only second is useful
% third = sym(dot(unit_q1q2,unit_q1q));
third = sym(dot(q1q2,q1q)); % H thesh tou q panw sto q1q2 dinetai me to logo tvn aktinwn!

eqns = [first==0,second(2)==0,third==lamda];
vars = [q_x q_y q_z];
[nsol_q_x, nsol_q_y, nsol_q_z] = solve(eqns,vars);
% [nsol_q_x, nsol_q_y, nsol_q_z] = vpasolve(eqns,vars);
nsol_q_z = -nsol_q_z; % take symmetry solution
nsol_q_x = nsol_q_z-0.1610; % solve again for symmetry solution
%% Points q lie in the plane of intersection between the 2 spheres
% q = [nsol_q_x nsol_q_y nsol_q_z]'; % wrong
% q = double(q);

figure(point2figure)
q_pts_a = intersectLineSphere(q1q2, Sph1); % q_pts_a(2,:)'
q_pts_b = intersectLineSphere(q1q2, Sph2); % q_pts_b(1,:)'
drawPoint3d(q_pts_a, 'mx'), hold on
drawPoint3d(q_pts_b, 'kx'), hold on

Q_c = (q_pts_a(2,:)'+ q_pts_b(1,:)')/2;

figure(point2figure)
scatter3(Q_c(1), Q_c(2), Q_c(3),'*b'), hold on
text(Q_c(1), Q_c(2), Q_c(3),'\leftarrow Q_c'), hold on

% Plane of intersection of two spheres. 
% In this plane lies the circle of intersection
Q_plane = createPlane(Q_c', unit_q1q2);
drawPlane3d(Q_plane,'facecolor', 'c','facealpha', '0.5'), hold on

%  Q_circle is the intersection between Sphere 2 or 1 and Q_plane
Q_circle = intersectPlaneSphere(Q_plane, Sph2);
drawCircle3d(Q_circle), hold on

%% Sphere 3 test - Must touch in {0,1,2} points q of Q_circle
%% Intersection points between Q_circle  and Sphere 3 are to be found
% Intersection points q can be: 0,1,2,Inf !!!

radius3 = distancePoints3d(p(1:3)', r1');
Sph3 = [r1(1) r1(2) r1(3) radius3];
figure(point2figure)
drawSphere(Sph3,'facecolor', 'b','facealpha', '0.5'), hold on

% Calculate points q as the intersection points between circle3D and Sphere
% to be continued...
%  Q_3 is the intersection between Sphere 3 and Q_plane
Q_3 = intersectPlaneSphere(Q_plane, Sph3);
drawCircle3d(Q_3), hold on

% Points q are the intersection points between Q_circle and Q_3
Tf = [Q_plane(4:6); Q_plane(7:9); unit_q1q2]; % rotation matrix of circles in 3D space given from plane equation
iTf = inv(Tf);

% u_points = intersectCircles([Q_circle(1) Q_circle(2) Q_circle(4)],[Q_3(1) Q_3(2) Q_3(4)]); %in Q_plane
% u_points(1,3) = 0;
% u_points(2,3) = 0;
% 
% q_points(:,1) = iTf*u_points(1,:)' + iTf*Q_c; % in Cartesian Space
% q_points(:,2) = iTf*u_points(2,:)' + iTf*Q_c;

% q_points = intersectCircles3d(Q_circle, Q_3, Tf, point2figure); % not ready yet

% q_points = [0.2176 -0.1527 0.5839; -0.1475 0.1414 -0.2265]'; % (0,0,1)
% q_points = [-0.0253 -0.1309 0.6836; -0.3765 0.1596 -0.1036]'; % (-2,-0.85,1)
q_points = [0.1576 -0.1029 0.4206; -0.0648 0.0836 -0.0805]'; % (1,-1,0)

figure(point2figure)
scatter3(q_points(1,1), q_points(2,1), q_points(3,1),'xb'), hold on
scatter3(q_points(1,2), q_points(2,2), q_points(3,2),'xb'), hold on
% scatter3(q_points(1,3), q_points(2,3), q_points(3,3),'xb'), hold on

%% Step 4: Since 2 points q are determined, generalized subproblem 2 is used twice
IKPfigure_1 = figure('Name','IKP Graphics for Gen 2 in Gen3 for q1');
IKPfigure_2 = figure('Name','IKP Graphics for Gen 2 in Gen3 for q2');
[theta_graph1, theta_anal1] = generalized_subproblem2(xi1,xi2,p(1:3,1),q_points(:,1),IKPfigure_1);
[theta_graph2, theta_anal2] = generalized_subproblem2(xi1,xi2,p(1:3,1),q_points(:,2),IKPfigure_2);

theta_anal = [theta_anal1 theta_anal2];
theta_graph = [theta_graph1 theta_graph2];

%% Step 6: Analytical solution using paper equations
% Following equations are presented in p.273 of the paper

% Not ready yet....
end