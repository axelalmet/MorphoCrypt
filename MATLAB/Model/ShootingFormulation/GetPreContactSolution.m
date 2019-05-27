function Sol = GetPreContactSolution(initValues, solMesh, solOld, modelParams, solOptions)

% Set the derivative function
DerivFun = @(x, M) RemodellingFoundationOdes(x, M, solOld, modelParams);

yL = initValues(1);
nxL = initValues(2);
mL = initValues(3);

initPoints = [0, 0, yL, nxL, 0, 0, mL];

% Integrate the system.
solution = ode45(DerivFun, [0, 1], initPoints, solOptions);

Sol.x = solMesh;
Sol.y = deval(solution, solMesh);



