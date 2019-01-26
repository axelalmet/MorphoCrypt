function Residuals = KellerModelHalfIntervalPreContactBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact at a single point.
% Specifically, we split the problem into three intervals and introduce
% matching conditions. This results in 21 boundary conditions.

% Obtain the relevant model parameters
L = modelParams.L;

Residuals = [Ml(1), Ml(2), Ml(3), Ml(5), Mr(2), Ml(6), Mr(6) - L];
         
end