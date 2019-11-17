function Sol = GetHomeostasisSolution(initValues, solMesh, modelParams, solOptions)

% Set the derivative function
DerivFun = @(x, M) UnbuckledExtensibleRodHomeostasisOdes(x, M, modelParams);

% Integrate the system.
solution = ode15s(DerivFun, [0, 1], initValues, solOptions);

Sol.x = solMesh;
Sol.y = deval(solution, solMesh);



