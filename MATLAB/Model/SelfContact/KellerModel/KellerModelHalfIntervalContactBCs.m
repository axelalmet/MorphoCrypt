function Residuals = KellerModelHalfIntervalContactBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact along a region

L = modelParams.L;

Residuals = [Ml(1, 1), ... S(0) = 0
             Ml(2, 1), ... x(0) = 0
             Ml(3, 1), ... y(0) = 0
             Ml(5, 1), ... N(0) = 0
             Ml(6, 1), ... theta(0) = 0
             Ml(1, 2) - Mr(1, 1), ... continuity of S
             Mr(2, 1), ... x(s1) = 0
             Ml(2, 2) - Mr(2, 1), ... continuity of x
             Ml(3, 2) - Mr(3, 1), ... continuity of y
             Ml(4, 2) - Mr(4, 1), ... continuity of Q
             Ml(5, 2) - (Mr(5, 1) - Mr(9, 1)), ... jump in N
             Mr(6, 1) - 0.5*pi, ... theta(s1) = pi/2
             Ml(6, 2) - Mr(6, 1), ... continuity of theta
             Ml(7, 2) - Mr(7, 1), ... continuity of k
             Ml(9, 2) - Mr(9, 1), ... continuity of fC
             Mr(1, 2) - L, ... S(L) = L
             Mr(2, 2), ... x(L) = 0
             Mr(6, 2) - L]; % theta(L) = L
         
end