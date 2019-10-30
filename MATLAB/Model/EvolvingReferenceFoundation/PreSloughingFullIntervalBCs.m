function Residuals = PreSloughingFullIntervalBCs(Ml, Mr, parameters)

y0 = parameters.y0;
L = parameters.L;

Residuals = [Ml(1), Ml(2), Ml(3) - y0, Ml(6), Mr(2) - L, Mr(3) - y0, Mr(6)];
