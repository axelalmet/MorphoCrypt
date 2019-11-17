function Sol = GetBuckledExtensibleHomeostasisSolution(initValues, solMesh, modelParams, solOld, solOptions)

% Set the derivative function
DerivFun = @(x, M) BuckledExtensibleRodHomeostasisOdes(x, M, solOld, modelParams);

% Integrate the system.
solution = ode15s(DerivFun, solMesh, initValues, solOptions);

Sol.x = solMesh;
Sol.y = deval(solution, solMesh);



