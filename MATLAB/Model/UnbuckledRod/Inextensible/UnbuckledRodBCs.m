function Residuals = UnbuckledRodBCs(Ml, Mr, modelParams)

l = modelParams.l;

Residuals = [Ml(1), Ml(2), Mr(3), Mr(4)];
