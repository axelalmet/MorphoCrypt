function Residuals = KellerModelContactBCs(Ml, Mr, modelParams)
% Boundary conditions specified for self-contact at a single point.
% Specifically, we split the problem into three intervals and introduce
% matching conditions. This results in 21 boundary conditions.

% Obtain the relevant model parameters
L = modelParams.L;

Residuals = [Ml(1), Ml(2), Ml(3), Mr(2), Ml(5), Ml(6), Mr(6) - L]; % theta3(L) = 0
         
end