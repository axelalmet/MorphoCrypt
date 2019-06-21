function Residuals = UnbuckledExtensibleRodSpringBCsWithSloughing(Ml, Mr, modelParams)

Ks = modelParams.Ks; % Stiffness of spring at left end

Residuals = [Ml(1), Ml(3) - Ks*Ml(1), Mr(2) - 1];
