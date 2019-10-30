function Residuals = SelfContactRegionNormalPressureBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact at a single point.
% Specifically, we split the problem into three intervals and introduce
% matching conditions. This results in 21 boundary conditions.

% Obtain the relevant model parameters
w0 = modelParams.w0;
y0 = modelParams.y0;
A0 = modelParams.A0;
L = modelParams.L;

Residuals = [Ml(1, 1), ... % S(0) = 0
             Ml(2, 1), ... % x(0) = 0
             Ml(5, 1), ... % G(0) = 0
             Ml(6, 1), ... % theta(0) = 0
             Ml(9, 1), ... % A(0) = 0
             Mr(2, 1) - (w0), ... % contact condition x = w0
             Ml(2, 2) - (w0), ... % contact condition x = w0
             Mr(6, 1) + 0.5*pi, ... % contact condition theta = -pi/2
             Ml(6, 2) + 0.5*pi, ... % contact condition theta = -pi/2
             Mr(7, 1), ... % contact condition m = 0
             Ml(7, 2), ... % contact condition m = 0
             Ml(1, 2) - Mr(1, 1), ... 
             Ml(8, 2) - Mr(8, 1), ... % continuity in sc
             Ml(9, 2) - Mr(9, 1), ... % continuity in A
             Mr(9, 2) - A0, ... A(sc) = A0
             Ml(10, 2) - Mr(10, 1), ... % continuity in p
             Mr(1, 2) - L, ... S(L) = L
             Mr(2, 2) - 0.5, ... % x(L) = 0.5
             Mr(3, 2) - y0, ... % y(L) = y0
             Mr(6, 2)];  % theta(L) = 0

         
end