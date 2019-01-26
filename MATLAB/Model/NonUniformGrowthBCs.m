function Residuals = NonUniformGrowthBCs(Ml, Mr, parameters)

L = parameters.L;
y0 = parameters.y0;

Residuals = [Ml(1), Ml(2), Mr(2) - L, Ml(3) - y0, Mr(3) - y0, Ml(6), Mr(6)];
