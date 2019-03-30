function angles = inverse_MMD_2(gd,pi,xi_ai,Js,xi_pi,tpi,exp_pi,ref2figure,point2figure,gst0)
% computes the IKP for 3 DoF MMD with 6 pseudojoints
% gd: is user-desired position and orientation of TCP
% pi: specified transformed points on axes, given by fkp, 3x9 matrix
% tpi: pseudojoint angles, 6x1 vector
% expi: pseudojoint exponenetial matrices, (4x4)x6 matrix

% load Murray kinematics
addpath('/home/nikos/matlab_ws/kinematics')
% load Generalized Subproblem 2 solution 
addpath('/home/nikos/matlab_ws/project2')

%% Points on reference anatomy and configuration
pi_ref = [           0   -0.0001   -0.0103   -0.0105   -0.0108   -0.0212   -0.0218   -0.0223   -0.0325   -0.0084    0.0966    0.1711;...
                     0         0   -0.0690   -0.0000   -0.0000   -0.0570   -0.0000   -0.0000   -0.0000    0.0001    0.0002    0.0002;...
                0.1610    0.2250    0.2900    0.3790    0.4430    0.5079    0.5969    0.6609    0.6928    0.7398    0.7406    0.7412];

%% Analytical Solution for 3 DoF Metamorphic Manipulator
%  based on solution written in Thursday 21.2.2019

% Specify FK map for MMD: eq.(I)
% gd(θi,θpi) =
% exp(ξa1*θ1)*Π(j,n1){exp(ξpj*θpj)}*exp(ξa2*θ2)*Π(j,n2){exp(ξpj*θpj)}*exp(ξa1*θ1)*Π(j,n3){exp(ξpj*θpj)}*gst0
% 
Pi1 = exp_pi(:,:,1)*exp_pi(:,:,2);
Pi2 = exp_pi(:,:,3)*exp_pi(:,:,4);
Pi3 = exp_pi(:,:,5)*exp_pi(:,:,6);

% We isolate exponenetials  eq.(II)
% gd*gst0^{-1} = exp(ξa1*θ1)*Π(j,n1){exp(ξpj*θpj)}*exp(ξa2*θ2)*Π(j,n2){exp(ξpj*θpj)}*exp(ξa3*θ3)*Π(j,n3){exp(ξpj*θpj)}

% Since Π(j,n3){exp(ξpj*θpj)} = Π3 is a known 4x4 matrix of exponentials
% multiplication eq.(II) becomes:
% gd*gst0^{-1}*Π3^{-1} =
% exp(ξa1*θ1)*Π(j,n1){exp(ξpj*θpj)}*exp(ξa2*θ2)*Π(j,n2){exp(ξpj*θpj)}*exp(ξa3*θ3) eq.(IΙI)
g = gd*inv(gst0)*inv(Pi3);

% So eq.(III) is g = exp(ξa1*θ1)*Pi1*exp(ξa2*θ2)*Pi2*exp(ξa3*θ3) eq.(IΙI)

% We transform ξa3 to ξa3' using the Adjoint transform
% g*Pi2^{-1} = exp(ξa1*θ1)*Pi1*exp(ξa2*θ2)*( Pi2*exp(ξa3*θ3)*Pi2^{-1} )
%                                           =========Adjoint==========
g = g*inv(Pi2);

% Re-create the ξa3 using the correct rotation vector and position
figure(point2figure) % Draw ξa3 in test anatomy
% xi_ai(:,3) = createtwist(xi_ai(4:6,3),pi(:,6));
% xi_ai(:,3) = createtwist(Js(4:6,3),pi(:,6)); %test11
axis_xi_a3 = drawtwist(xi_ai(:,3)); hold on

xi_a3_new = ad(Pi2)*xi_ai(:,3);
% xi_a3_new = xi_ai(:,3); % test sklhro

% Draw ξa3'
figure(point2figure)
axis_xi_a3_new = drawtwist(xi_a3_new); hold on
% pa3new_1 - axis_xi_a3_new - pa3new_2
pa3new_1 = [axis_xi_a3_new.XData(1) axis_xi_a3_new.YData(1) axis_xi_a3_new.ZData(1) ]';
pa3new_2 = [axis_xi_a3_new.XData(2) axis_xi_a3_new.YData(2) axis_xi_a3_new.ZData(2) ]';
text(pa3new_2(1), pa3new_2(2), pa3new_2(3),'\leftarrow ξa3_n'), hold on

% exp(ξa3'*θ3) = Pi2*exp(ξa3*θ3)*Pi2^{-1} eq.(IV)

% So, eq.(III) now becomes: 
% g*Pi2^{-1} = exp(ξa1*θ1)*Pi1*exp(ξa2*θ2)*exp(ξa3'*θ3) eq.(V)

%% ~Now starts the new challenge~
% We use the Reverse Adjoint transform to eliminate Pi1 in the right side
% and transform ξa1 to ξa1'
% Pi1^{-1}*g*Pi2^{-1} = ( Pi1^{-1}*exp(ξa1*θ1)*Pi1 )*exp(ξa2*θ2)*exp(ξa3'*θ3)
%                        ======Reverse Adjoint======
% eq.(V) now becomes:
% Pi1^{-1}*g*Pi2^{-1} = exp(ξa1'*θ1)*exp(ξa2*θ2)*exp(ξa3'*θ3) eq.(VI)
g1 = inv(Pi1)*g;
% From Murray page 58 we have:
% g^{-1)* (ξ) * g = [ Ad(g^{-1})*ξ ]^, so Exponential Reverse Adjoint = the Wedge
% transform of the Inverse Adjoint of the twist Exponential
% Applying in exponentials the same method of Murray as presented in page 94

% Re-create the ξa1 using the correct rotation vector and position
% xi_ai(:,1) = createtwist(xi_ai(4:6,1),pi(:,1));
% xi_ai(:,1) = createtwist(Js(4:6,1),pi(:,1)); %test11
axis_xi_a1 = drawtwist(xi_ai(:,1)); hold on % draw ξa1 before new
xi_a1_new = ad(inv(Pi1))*xi_ai(:,1);
% xi_a1_new = xi_ai(:,1); % test sklhro

% Draw ξa1'
figure(point2figure)
axis_xi_a1_new = drawtwist(xi_a1_new); hold on
% pa1new_1 - axis_xi_a1_new - pa1new_2
pa1new_1 = [axis_xi_a1_new.XData(1) axis_xi_a1_new.YData(1) axis_xi_a1_new.ZData(1) ]';
pa1new_2 = [axis_xi_a1_new.XData(2) axis_xi_a1_new.YData(2) axis_xi_a1_new.ZData(2) ]';
text(pa1new_2(1), pa1new_2(2), pa1new_2(3),'\leftarrow ξa1_n'), hold on

% We search the solution set for eq.(VI)
% g1 = exp(ξa1'*θ1)*exp(ξa2*θ2)*exp(ξa3'*θ3) eq.(VI)

%  ===================!!!!!!!!!!!==================
%% Up to this point it is the same as inverse_MMD.m
%  ===================!!!!!!!!!!!==================

% We multiply eq.(VI) with point p that does NOT belong to either ξai
% g1*p = exp(ξa1'*θ1)*exp(ξa2*θ2)*exp(ξa3'*θ3)*p eq.(VII)

% p = [pi_ref(:,12); 1]; % gst0 as the arbitrary point 
p = [pi(:,12); 1]; %  gst point in new anatomy

% We substract point p1' ε ξa1':
% g1*p - p1' = exp(ξa1'*θ1)*exp(ξa2*θ2)*exp(ξa3'*θ3)*p-p1'
% g1*p - p1' = exp(ξa1'*θ1)*( exp(ξa2*θ2)*exp(ξa3'*θ3)*p-p1')
% and by applying norms to both sides:

% δ1 = || ( exp(ξa2*θ2)*exp(ξa3'*θ3)*p-p1') || eq.(VIIIa)
dp1 = [pa1new_1; 1]; % define p1'
delta1 = norm(g1*p - dp1);

% We substract point p1'' ε ξa1':
% g1*p - p1'' = exp(ξa1'*θ1)*exp(ξa2*θ2)*exp(ξa3'*θ3)*p-p1''
% g1*p - p1'' = exp(ξa1'*θ1)*( exp(ξa2*θ2)*exp(ξa3'*θ3)*p-p1'')
% and by applying norms to both sides:

% δ2 = || ( exp(ξa2*θ2)*exp(ξa3'*θ3)*p-p1'') || eq.(VIIIb)
ddp1 = [pa1new_2; 1]; % define p1''
delta2 = norm(g1*p - ddp1);

%% Solve for (θ2,θ3')
% eq.(VIIIa) and .(VIIIb) are the Generalized Subproblem 3'
% xi1 = ξa2, xi2 = ξa3', p = p, q1=p1', q2=p1'' , delta1 = δ1, delta2 = δ2 

% Re-create the ξa2 using the correct rotation vector and position
% xi_ai(:,2) = createtwist(xi_ai(4:6,2),pi(:,3));
% xi_ai(:,2) = createtwist(Js(4:6,2),pi(:,3)); %test11

% Draw ξa2
figure(point2figure)
axis_xi_a2 = drawtwist(xi_ai(:,2)); hold on
% pa2_1 - axis_xi_a2_new - pa2_2
pa2_1 = [axis_xi_a2.XData(1) axis_xi_a2.YData(1) axis_xi_a2.ZData(1) ]';
pa2_2 = [axis_xi_a2.XData(2) axis_xi_a2.YData(2) axis_xi_a2.ZData(2) ]';
text(pa2_2(1), pa2_2(2), pa2_2(3),'\leftarrow ξa2'), hold on

IKPfigure = figure('Name','IKP Graphics');

[theta_graph, theta_anal] = generalized_subproblem3(xi_ai(:,2),xi_a3_new,p,dp1,ddp1,delta1,delta2,IKPfigure);

%% Solve for (θ1')
% 4 solution sets are extracted for (θ2,θ3), from the Generalized Paden Kahan Subproblem 3
% We return to equation (VI):
% g1 = exp(ξa1'*θ1)*exp(ξa2*θ2)*exp(ξa3'*θ3) eq.(VI)

% Solve for the first set
t1 = theta_anal(:,1);
% Since the exponentials of
% exp(ξa2*θ2)*exp(ξa3'*θ3) can be calculated, we get:
% g1*exp(-ξa3'*θ3)*exp(-ξa2*θ2) = exp(ξa1'*θ1) eq.(VII)
g2 = g1*twistexp(-xi_a3_new,t1(2))*twistexp(-xi_ai(:,2),t1(1));
% We multiply from right side with point p ( NOT ε )in ξa1'
% So, eq.(VII) becomes:
% g2*p =exp(ξa1'*θ1)*p eq.(VIII) which is Subproblem1
q = g2*p;
t1(3) = subproblem1(xi_a1_new,p(1:3,1),q(1:3,1),pa1new_1);
t1 = [t1(3) t1(1) t1(2)]';

% Solve for the second set (same procedure)
t2 = theta_anal(:,2);
g2 = g1*twistexp(-xi_a3_new,t2(2))*twistexp(-xi_ai(:,2),t2(1));
q = g2*p;
t2(3) = subproblem1(xi_a1_new,p(1:3,1),q(1:3,1),pa1new_1);
t2 = [t2(3) t2(1) t2(2)]';

% Solve for the third set (same procedure)
t3 = theta_anal(:,3);
g2 = g1*twistexp(-xi_a3_new,t3(2))*twistexp(-xi_ai(:,2),t3(1));
q = g2*p;
t3(3) = subproblem1(xi_a1_new,p(1:3,1),q(1:3,1),pa1new_1);
t3 = [t3(3) t3(1) t3(2)]';
% Solve for the fourth set (same procedure)
t4 = theta_anal(:,4);
g2 = g1*twistexp(-xi_a3_new,t4(2))*twistexp(-xi_ai(:,2),t4(1));
q = g2*p;
t4(3) = subproblem1(xi_a1_new,p(1:3,1),q(1:3,1),pa1new_1);
t4 = [t4(3) t4(1) t4(2)]';

%% Final results-must be 4 sets of [θ1 θ2 θ3]'
angles = [t1 t2 t3 t4];

end