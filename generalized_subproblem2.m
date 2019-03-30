function [theta_graph, theta_anal] = generalized_subproblem2(xi1,xi2,p,q,point2figure)
% Based on paper of Dimovski, "Algoritmic approach to geometric solution of
% generalized PK subproblem and its extension", methodology presented in page 4
% Gives solution (θ1,θ2) to the equation 
% exp(^ξ1*θ1)*exp(^ξ2*θ2)*p = q; (eq.1)
% for ξ1, ξ2 disjoint lines!
% for ξ1, ξ2 not parallel - > to be fixed!!!

% load Murray kinematics
addpath('/home/nikos/matlab_ws/kinematics/robotlinks')
addpath('/home/nikos/matlab_ws/kinematics/screws') 
addpath('/home/nikos/matlab_ws/kinematics/util')
% load Geometry 2-3D libraries
addpath('/home/nikos/matlab_ws/geom3d/geom3d')
addpath('/home/nikos/matlab_ws/geom2d/utils')

%% Step1: Determine the closest points r1,r2 between the axes of twists ξ1,ξ2
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
scatter3(q(1), q(2), q(3),'*m'), hold on
text(q(1), q(2), q(3),'\leftarrow q'), hold on
scatter3(p(1), p(2), p(3),'*m'), hold on
text(p(1), p(2), p(3),'\leftarrow p'), hold on

%% Step2: Determine the two planes that are normal to the unit vectors and pass through p,q and their intersecting line l(t)
% Must be Σ1: w1 the normal vector to plane surface, q ε plane
% Must be Σ2: w2 the normal vector to plane surface, p ε plane

% Based on figure 1 of page 4
w1=xi1(4:6,:);
w2=xi2(4:6,:);

% First find second point of plane, that is to determine projections of p,q to ξ2,ξ1
v = q-r1;
s1 = r1 + w1*(w1'*v); % projection of q to ξ1 and Σ1
u = p-r2;
s2 = r2 + w2*(w2'*u); % projection of p to ξ2 and Σ2

% Plane evaluation, q,s1 must lie on the same plane
eval_Sigma1 = dot((q-s1),w1);
% Plane evaluation, p,s2 must lie on the same plane
eval_Sigma2 = dot((p-s2),w2);

% We now have 2 known points and the normal vector

% fig = gcf;
% ax = fig.CurrentAxes;
% Sigma11 = createPlane(q', w1');
% h_Sigma1 = drawPlane3d(ax,Sigma11), hold on

figure(point2figure)
% [Sigma1,a1,b1] = plot_plane(q, w1, point2figure), hold on
Sigma1 = createPlane(q(1:3)', w1');
drawPlane3d(Sigma1,'facecolor', 'g','facealpha', '0.5'), hold on

% [Sigma2,a2,b2] = plot_plane(p, w2, point2figure), hold on
Sigma2 = createPlane(p(1:3)', w2');
drawPlane3d(Sigma2,'facecolor', 'y','facealpha', '0.5'), hold on

% Computes intersecting line of the form: l(t)= l0+t*D;
[l0,lD,check]=plane_intersect(w1,q,w2,p);
syms t;
l_t = l0' + t*lD;
figure(point2figure)
fplot3(l_t(1), l_t(2), l_t(3), [-1000,1000],'LineWidth',2), hold on

%% Step 3:Determine projections of p,q to ξ2,ξ1 --> dp,dq respectively
% [qx,qy,qz] = project_point2plane(a1(1),a1(2),a1(3),b1,q(1),q(2),q(3)) % projection of q to ξ1,Σ2
% dq = [qx qy qz]';
point2line = q';
line_xi1(1,1:3) = p1_2; %p1_1 before
line_xi1(1,4:6) = w1;
dq = projPointOnLine3d(point2line, line_xi1); % inputs: point2line 1x3, line 1x6
% [px,py,pz] = project_point2plane(a2(1),a2(2),a2(3),b2,p(1),p(2),p(3)) % projection of p to ξ2,Σ2
% dp = [px py pz]';
point2line = p';
line_xi2(1,1:3) = p2_2; %p2_1 before
line_xi2(1,4:6) = w2;
dp = projPointOnLine3d(point2line, line_xi2); % inputs: point2line 1x3, line 1x6

figure(point2figure)
axis([-10, 10, -10, 10, -10, 10]);
scatter3(dq(1), dq(2), dq(3),'k'), hold on
text(dq(1), dq(2), dq(3),'\leftarrow dq'), hold on
scatter3(dp(1), dp(2), dp(3),'k'), hold on
text(dp(1), dp(2), dp(3),'\leftarrow dp'), hold on

% Den bgainoun ta swsta, pairnw s1->dq s2->dp
% figure(point2figure)
% scatter3(s1(1), s1(2), s1(3),'k'), hold on
% text(s1(1), s1(2), s1(3),'\leftarrow s1'), hold on
% scatter3(s2(1), s2(2), s2(3),'k'), hold on
% text(s2(1), s2(2), s2(3),'\leftarrow s2'), hold on

% R1 = norm(q'-dq);
R1 = distancePoints3d(q', dq);
% Determine the intersection between l(t) and red circle (s1,R1), red
% circle lies in the plane normal to twist ξ1
figure(point2figure)
% fig = gcf;
% ax = fig.CurrentAxes;
Sph1 = [dq(1) dq(2) dq(3) R1];
Circ1 = intersectPlaneSphere(Sigma1, Sph1);
drawCircle3d(Circ1, 'LineWidth', 2, 'Color', 'r'), hold on
% H1 = drawCircle3d([dq(1) dq(2) dq(3) R1], [angZ1 angX1 angY1], 'LineWidth', 2, 'Color', 'r'), hold on
% H1 = circlePlane3D( [dq(1) dq(2) dq(3)], w1', R1, 0.05, 0, 'r', point2figure), hold on

line1(1,1:3) = l0;
line1(1,4:6) = lD;
% circle1 = [dq(1) dq(2) dq(3) R1];
c12 = intersectLineSphere(line1, Circ1); %  Returns the two points 2x3 array: 1x3 each point
figure(point2figure) % shows points of intersection between l(t) and red circle
scatter3(c12(1,1), c12(1,2), c12(1,3),'go'), hold on
text(c12(1,1), c12(1,2), c12(1,3),'\leftarrow c1'), hold on
scatter3(c12(2,1), c12(2,2), c12(2,3),'go'), hold on
text(c12(2,1), c12(2,2), c12(2,3),'\leftarrow c2'), hold on


% R2 = norm(p'-dp);
R2 = distancePoints3d(p', dp);
% Determine the intersection between l(t) and red circle (s2,R2)
figure(point2figure)
Sph2 = [dp(1) dp(2) dp(3) R2];
Circ2 = intersectPlaneSphere(Sigma2, Sph2);
drawCircle3d(Circ2, 'LineWidth', 1, 'Color', 'r'), hold on
% H2 = drawCircle3d([dp(1) dp(2) dp(3) R2], [angZ2 angX2 angY2], 'LineWidth', 1, 'Color', 'r'), hold on
% H2 = circlePlane3D( [dp(1) dp(2) dp(3)], w2', R2, 0.05, 0,'r--', point2figure), hold on

line2(1,1:3) = l0;
line2(1,4:6) = lD;
% circle2 = [dp(1) dp(2) dp(3) R2];
c34 = intersectLineSphere(line2, Circ2); %  Returns the two points 2x3 array: 1x3 each point
figure(point2figure) % shows points of intersection between l(t) and red circle
scatter3(c34(1,1), c34(1,2), c34(1,3),'m'), hold on
text(c34(1,1), c34(1,2), c34(1,3),'\leftarrow c3'), hold on
scatter3(c34(2,1), c34(2,2), c34(2,3),'m'), hold on
text(c34(2,1), c34(2,2), c34(2,3),'\leftarrow c4'), hold on

%% Step 4 Conditions for solution existence - Graphical Solution
% Find the closest points from s1,s2 to l(t)
% First we compute the line twist
line_twist(1,1:3) = l0;
line_twist(1,4:6) = lD;

point2line = s1';
ds1 = projPointOnLine3d(point2line, line_twist); % inputs: point2line 1x3, line 1x6
figure(point2figure)
scatter3(ds1(1), ds1(2), ds1(3),'y*'), hold on
text(ds1(1), ds1(2), ds1(3),'\leftarrow ds1'), hold on

point2line = s2';
ds2 = projPointOnLine3d(point2line, line_twist); % inputs: point2line 1x3, line 1x6
figure(point2figure)
scatter3(ds2(1), ds2(2), ds2(3),'y*'), hold on
text(ds2(1), ds2(2), ds2(3),'\leftarrow ds2'), hold on

% Condition for θ1 solution
% eval_sol_t1 = |dq,ds1| <= |q-dq| = R1
% We draw the circle exist_sol_t1(dq,eval_sol_t2), lies on Σ1

% eval_sol_t1 = norm(s1-ds1);
eval_sol_t1 = distancePoints3d(dq, ds1);

figure(point2figure)
% H3 = circlePlane3D( [dq(1) dq(2) dq(3)], w1', eval_sol_t1, 0.05, 0,'b', point2figure), hold on
H3 = drawCircle3d([dq(1) dq(2) dq(3) eval_sol_t1], [Circ1(5) Circ1(6) Circ1(7)], 'LineWidth', 2, 'Color', 'b'), hold on

% For θ1 existence blue circle must be INSIDE RED circle!!!

% Condition for θ2 solution
% eval_sol_t2 = |dp,ds2| <= |p-dp| = R2
% We draw the circle exist_sol_t2(dp,eval_sol_t2), lies on Σ1
% eval_sol_t2 = norm(s2-ds2);
eval_sol_t2 = distancePoints3d(dp, ds2); 

figure(point2figure)
% H4 = circlePlane3D( [dp(1) dp(2) dp(3)], w2', eval_sol_t2, 0.05, 0,'b--', point2figure), hold on
H4 = drawCircle3d([dp(1) dp(2) dp(3) eval_sol_t2], [Circ2(5) Circ2(6) Circ2(7)], 'LineWidth', 1, 'Color', 'b'), hold on

% For θ2 existence blue circle must be INSIDE RED circle!!!

% For existance of solution set (θ1,θ2) MUST at least one point to be equal
% c12 = c34 | at least 1 !!!

%% Step 5: Analytical solution using c1,c2,c3,c4 from graphs

%% Analytical using Murray-TanXiao
% d: common normal is calculated in l.
% u = p-r2 is in l.
% v = q-r1 is in l.
alpha = ( ((w1'*w2)*w2'*u) - (w1'*v)) / ( (w1'*w2)^2-1 );
beta = ( ((w1'*w2)*w1'*v) - (w2'*u)) / ( (w1'*w2)^2-1 );
gamma2 = (norm(u)^2 - alpha^2 - beta^2 - 2*alpha*beta*w1'*w2 ) / ( norm(cross(w1,w2))^2 );
    if gamma2<0 % No point c exists
        % solution desn't exist
        disp('Solution doesnt exist');
        theta_anal(:,1) = [NaN NaN]';
        theta_anal(:,2) = [NaN NaN]';
%         theta_anal(:,3) = [NaN NaN]';
%         theta_anal(:,4) = [NaN NaN]';
    elseif  gamma2==0 % NOT READY, Only 1 point c exists! => ????? MUST INTERSECT ?????
        % this ELSE IF not ready!
        z1 = alpha*w1 + beta*w2 - norm(d)*cross(w1,w2);
        z2 = alpha*w1 + beta*w2;
        c1 = z1 + r1 ; % c1 must be equal to c2!!!!
        c2 = z2 + r2 ;
        %  We solve two times (eq.2)
        %  exp(^ξ2*θ2)*p = c1  
        theta2(1) = subproblem1(xi2,p,c1,r2);
        theta2(2) = subproblem1(xi2,p,c2,r2);
    %     theta2(2) = theta2(1);
        %  exp(^-ξ1*θ1)*q = c1; 
        theta1(1) = subproblem1(-xi1,q,c1,r1);
        theta1(2) = subproblem1(-xi1,q,c2,r1);
    %     theta1(2) = theta2(1);
        theta_anal(:,1) = [theta1(1) theta2(1)]';
        
    else % 2 point c exist: c1,2
        z1(:,1)= alpha*w1 + beta*w2 + (sqrt(gamma2)-norm(d))*cross(w1,w2); % z11
        z1(:,2)= alpha*w1 + beta*w2 + sqrt(gamma2)*cross(w1,w2); % z12
        c(:,1) = z1(:,1) + r1 ; 
%         c(:,2) = z1(:,2) + r2 ; % c(:,1)=c(:,2) if -norm(d)
        
%         figure(point2figure) % shows points of intersection between l(t) and red circle
%         scatter3(c(1,1), c(2,1), c(3,1),'go'), hold on
%         text(c(1,1), c(2,1), c(3,1),'\leftarrow ca1'), hold on
        
        z2(:,1)= alpha*w1 + beta*w2 + (-sqrt(gamma2)-norm(d))*cross(w1,w2); % z21
        z2(:,2)= alpha*w1 + beta*w2 - sqrt(gamma2)*cross(w1,w2); % z22
        c(:,2) = z2(:,1) + r1 ; 
%         c(:,4) = z2(:,2) + r2 ; % c(:,3)=c(:,4) if -norm(d)

        % Which to trust test graph vs analytic => both error~10^(-2,-3)
%         c = [c12' c34'];
        
%         figure(point2figure) % shows points of intersection between l(t) and red circle
%         scatter3(c(1,2), c(2,2), c(3,2),'go'), hold on
%         text(c(1,2), c(2,2), c(3,2),'\leftarrow ca2'), hold on

        for i=1:2
            %  Solves exp(^ξ2*θ2)*p = c1,2 && exp(^-ξ1*θ1)*q = c1,2
            theta_anal(2,i) = subproblem1(xi2,p,c(:,i),r2);
            theta_anal(1,i) = subproblem1(-xi1,q,c(:,i),r1);
        end

    end
%% Exploiting graphical extraction of ci points, not ready yet!
% Solution exists only if for, A={c12} , B={c34} : A tomh B <> kenou sunolou
A = c12
B = c34
[C,ia,ib] = intersect(A,B,'rows')

    if isempty(C)
        % solution desn't exist
        disp('Solution doesnt exist');
        theta_graph(:,1) = [NaN NaN]';
        theta_graph(:,2) = [NaN NaN]';
    else
        % solutions exist
        
        %defife c1,c2
        c1 = C(1,:)';
        % if 1 intersection point
        theta11=subproblem1(-xi1,q,c1,r1);
        theta21=subproblem1(xi2,p,c1,r2);
        theta_graph(:,1) = [theta11 theta21]';
        if length(C(:,1))==2
            disp('Two solutions exist');
            c2 = C(2,:)';
            % if 2 intersection points
            theta12=subproblem1(-xi1,q,c2,r1);
            theta22=subproblem1(xi2,p,c2,r2);
            theta_graph(:,2) = [theta12 theta22]';
        else
            disp('One solution exists');
        end
        
    end

end % function