function Residuals = UnbuckledExtensibleRodBCs(Ml, Mr, modelParams)

Residuals = [Ml(1), Mr(1) - 1];
