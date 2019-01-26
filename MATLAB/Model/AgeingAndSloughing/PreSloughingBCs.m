function Residuals = PreSloughingBCs(Ml, Mr, parameters)

y0 = parameters.y0;

Residuals = [Ml(1), Ml(2), Ml(5), Ml(6), Mr(2) - 0.5, Mr(3) - y0, Mr(6)];
