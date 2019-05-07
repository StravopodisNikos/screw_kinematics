function [metJs, metJb, Error] = metamorphic_jacobians_6dof(xi_ai,exp_ai, Pi, gst0, gst) 
% Calculates Js = Jspatial and Jb = Jbody
% of 6 d.o.f. serial modular manipulator with 2x3 pseudojoints

% gst0(1:3,4) = gst0(1:3,4)*10^(-3); % [mm-->m]
% gst(1:3,4) = gst(1:3,4)*10^(-3); % [mm-->m]

% Jspatial calculation (Murray p.116)
metJs = zeros(6,6);
metJs(:,1)=xi_ai(:,1);
metJs(:,2)=ad(exp_ai(:,:,1)*Pi(:,:,1))*xi_ai(:,2);
metJs(:,3)=ad(exp_ai(:,:,1)*Pi(:,:,1)*exp_ai(:,:,2)*Pi(:,:,2))*xi_ai(:,3);
metJs(:,4)=ad(exp_ai(:,:,1)*Pi(:,:,1)*exp_ai(:,:,2)*Pi(:,:,2)*exp_ai(:,:,3)*Pi(:,:,3))*xi_ai(:,4);
metJs(:,5)=ad(exp_ai(:,:,1)*Pi(:,:,1)*exp_ai(:,:,2)*Pi(:,:,2)*exp_ai(:,:,3)*Pi(:,:,3)*exp_ai(:,:,4))*xi_ai(:,5);
metJs(:,6)=ad(exp_ai(:,:,1)*Pi(:,:,1)*exp_ai(:,:,2)*Pi(:,:,2)*exp_ai(:,:,3)*Pi(:,:,3)*exp_ai(:,:,4)*exp_ai(:,:,5))*xi_ai(:,6);

% Jbody calculation (Murray p.117)
% Xrhsimopoioume thn idiothta gia ton metasxhmatismo
% tou SE(3) apo sel 55/Murray
% (Adg)^(-1) = Ad(g)^(-1)
metJb = zeros(6,6);
metJb(:,1) = ad(inv(exp_ai(:,:,1)*Pi(:,:,1)*exp_ai(:,:,2)*Pi(:,:,2)*exp_ai(:,:,3)*Pi(:,:,3)*exp_ai(:,:,4)*exp_ai(:,:,5)*exp_ai(:,:,6)*gst0))*xi_ai(:,1);
metJb(:,2) = ad(inv(exp_ai(:,:,2)*Pi(:,:,2)*exp_ai(:,:,3)*Pi(:,:,3)*exp_ai(:,:,4)*exp_ai(:,:,5)*exp_ai(:,:,6)*gst0))*xi_ai(:,2);
metJb(:,3) = ad(inv(exp_ai(:,:,3)*Pi(:,:,3)*exp_ai(:,:,4)*exp_ai(:,:,5)*exp_ai(:,:,6)*gst0))*xi_ai(:,3);
metJb(:,4) = ad(inv(exp_ai(:,:,4)*exp_ai(:,:,5)*exp_ai(:,:,6)*gst0))*xi_ai(:,4);
metJb(:,5) = ad(inv(exp_ai(:,:,5)*exp_ai(:,:,6)*gst0))*xi_ai(:,5);
metJb(:,6) = ad(inv(exp_ai(:,:,6)*gst0))*xi_ai(:,6);

% Jb = Jb.*10^(-3);
% Js = Js.*10^(-3);

% Gia epalh8eush xrhsimopoioume th sxesh (3.56)/Murray
% Js = (Adg)*Jb
Adg = ad(gst);
Jepal = Adg*metJb;
Error = Jepal - metJs; % Eimaste swstoi
end