function Residuals = KellerModelHalfIntervalContactRegionBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact along a region

L = modelParams.L;
P = modelParams.P;

Residuals = [Ml(1, 1), ... s = 0
             Ml(2, 1), ... x(0) = 0
             Ml(3, 1), ... y(0) = 0
             Ml(5, 1), ... N(0) = 0
             Ml(6, 1), ... theta(0) = 0
             Mr(2, 1), ... x(s1) = 0
             Ml(2, 2) - Mr(2, 1), ... continuity of x at s1
             Ml(4, 2) - Mr(4, 1), ... continuity of Q s1
             Mr(6, 1) - 0.5*pi, ... theta(s1) = pi/2
             Ml(6, 2) - Mr(6, 1), ... theta(s1) = pi/2
             Mr(7, 1), ... k(s1) = 0
             Ml(7, 2) - Mr(7, 1), ... continuity in k
             Ml(3, 2) - (Mr(3, 1) + 2*Mr(8, 1)), ... jump in y
             Mr(1, 2) - L, ... S(L) = L
             Mr(2, 2), ... x(L) = 0
             Mr(6, 2) - L]; %theta(L) = L
         
end