function Residuals = EvaluatePreContactSolution(solGuess, solOld, modelParams, solOptions)

% Set the derivative function
DerivFun = @(x, M) RemodellingFoundationOdes(x, M, solOld, modelParams);

yL = solGuess(1);
nxL = solGuess(2);
mL = solGuess(3);

initPoints = [0, 0, yL, nxL, 0, 0, mL];

% Integrate the system.
Sol = ode45(DerivFun, [0, 1], initPoints, solOptions);

Residuals = norm(Sol.y([2, 3, 6], end) - [0.5; 0; 0], 2);


