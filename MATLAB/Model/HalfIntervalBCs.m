function Residuals = HalfIntervalBCs(Ml, Mr, modelParams)
% Boundary conditions specified for the rod on the half-interval

% Obtain the relevant model parameters
y0 = modelParams.y0;
L = modelParams.L;

Residuals = [Mr(1) - L, Ml(2), Mr(2) - L, Ml(3) - y0, Mr(5), Ml(6), Mr(6)];

end