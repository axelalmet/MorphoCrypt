function Residuals = UnbuckledExtensibleRodClampedBCsWithSloughing(Ml, Mr, modelParams)

Residuals = [Ml(1), Mr(1) - 1];  % standard clamped BCs for s(0) = 0, s(L) = 1,
