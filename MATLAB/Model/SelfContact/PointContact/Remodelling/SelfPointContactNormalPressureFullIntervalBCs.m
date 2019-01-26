function Residuals = SelfPointContactNormalPressureFullIntervalBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact at a single point.
% Obtain the relevant model parameters
w0 = modelParams.w0;
y0 = modelParams.y0;
L = modelParams.L;
A0 = modelParams.A0;

Residuals = [Ml(1, 1), ... % S1(0) = 0
             Ml(2, 1), ... % x1(0) = 0
             Ml(3, 1) - y0, ... % y1(0) = y0
             Ml(6, 1), ... % theta1(0) = 0
             Ml(9, 1), ...
             Ml(1, 2) - Mr(1, 1), ... % continuity in S
             Mr(2, 1) - (0.5*L - w0), ... % x1(S1) = L/2 - w; contact point
             Ml(2, 2) - (0.5*L - w0), ... % x2(S1) = L/2 - w; contact point
             Ml(3, 2) - Mr(3, 1), ... % y3(S2) = y2(S2)
             Ml(5, 2) - Mr(5, 1), ... % G2(S1) = G1(S2); continuity                
             Mr(6, 1) - 0.5*pi, ... % theta1(S1) = pi/2; contact point)
             Ml(6, 2) - 0.5*pi, ... % theta2(s1) = pi/2
             Ml(7, 2) - Mr(7, 1), ... % m(S2) = m(S2); continuity
             Ml(8, 2) - Mr(8, 1), ...
             Ml(9, 2) - Mr(9, 1), ... % continuity in A
             Ml(10, 2) - Mr(10, 1), ...
             Mr(9, 2) - A0, ... % continuity in p
             Ml(9, 3) - Mr(9, 2), ...
             Ml(1, 3) - Mr(1, 2), ...
             Mr(2, 2) - 0.5*L, ... % x3(L/2) = L/2
             Ml(2, 3) - Mr(2, 2), ...             
             Mr(5, 2), ... G(L/2) = 0
             Ml(5, 3), ...
             Mr(6, 2), ... theta(L/2) = 0
             Ml(6, 3), ...
             Ml(1, 4) - Mr(1, 3), ... % continuity in S
             Mr(2, 3) - (0.5*L + w0), ... % x1(S1) = L/2 - w; contact point
             Ml(2, 4) - (0.5*L + w0), ... % x2(S1) = L/2 - w; contact point
             Ml(3, 4) - Mr(3, 3), ... % y3(S2) = y2(S2)
             Ml(5, 4) - Mr(5, 3), ... % G2(S1) = G1(S2); continuity
             Mr(6, 3) + 0.5*pi, ... % theta1(S1) = pi/2; contact point)
             Ml(6, 4) + 0.5*pi, ... % theta2(s1) = pi/2
             Ml(7, 4) - Mr(7, 3), ... % m(S2) = m(S2); continuity
             Ml(9, 4) - Mr(9, 3), ... % continuity in A
             Ml(10, 4) - Mr(10, 3), ...
             Mr(1, 4) - L, ...
             Mr(2, 4) - L, ...
             Mr(3, 4) - y0, ...
             Mr(6, 4), ...
             Mr(9, 4)];  
         
end