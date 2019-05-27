function Residuals = UnbuckledExtensibleRodPreSloughingClampedBCs(Ml, Mr, modelParams)

L = modelParams.L;

Residuals = [Ml(1), Ml(2), Mr(2) - L];
