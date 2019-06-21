function Residuals = UnbuckledExtensibleRodClampedBCsWithSloughing(Ml, Mr, modelParams)

Residuals = [Ml(1), Ml(2), Mr(2) - 1];  % standard clamped BCs for s(0) = 0, s(L) = 1,
