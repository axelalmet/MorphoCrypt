function Sol = GetBuckledInextensibleHomeostasisSolution(initValues, solMesh, modelParams, solOld, solOptions)

% Set the derivative function
DerivFun = @(x, M) BuckledInextensibleRodHomeostasisOdes(x, M, solOld, modelParams);

% Integrate the system.
solution = ode15s(DerivFun, solMesh, initValues, solOptions);

Sol.x = solMesh;
Sol.y = deval(solution, solMesh);



