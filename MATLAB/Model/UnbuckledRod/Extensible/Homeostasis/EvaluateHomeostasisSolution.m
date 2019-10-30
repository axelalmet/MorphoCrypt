function Residuals = EvaluateHomeostasisSolution(solGuess, modelParams, solOptions)

% Set the derivative function
DerivFun = @(x, M) UnbuckledExtensibleRodHomeostasisOdes(x, M, modelParams);

% Integrate the system.
Sol = ode15s(DerivFun, [0, 1], solGuess, solOptions);

Residuals = norm(deval(Sol, 1, 2) - 0, 2);


