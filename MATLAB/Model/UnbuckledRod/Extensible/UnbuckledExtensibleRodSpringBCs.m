function Residuals = UnbuckledExtensibleRodSpringBCs(Ml, Mr, modelParams)

Ks = modelParams.Ks; % Stiffness of spring at left end

Residuals = [Ml(2) - Ks*Ml(1), Mr(1) - 1];
