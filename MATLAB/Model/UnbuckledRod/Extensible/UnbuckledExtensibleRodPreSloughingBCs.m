function Residuals = UnbuckledExtensibleRodPreSloughingClampedBCs(Ml, Mr, modelParams)

Residuals = [Ml(1), Ml(2), Mr(2) - 1];
