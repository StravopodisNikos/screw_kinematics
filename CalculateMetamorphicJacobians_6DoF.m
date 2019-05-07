function [metJs,metJb,metError] = CalculateMetamorphicJacobians_6DoF(config,xi_ai,tpi,xi_pi,Pi,gst0)
% Calculates Spatial and Body Jacobians
% Works for 6DoF with 2x3 PseudoJoints with the following structure:
% A1 - P1-P2 - A2 - P3-P4 - A3 - P5-P6 - A4 - A5 - A6

[gst,exp_ai] = MMD_POE_FKP_6DoF(config,xi_ai,tpi,xi_pi,gst0);
[metJs, metJb, metError] = metamorphic_jacobians_6dof(xi_ai,exp_ai, Pi, gst0, gst);

end