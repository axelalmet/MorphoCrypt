function Residuals = EvaluateSelfContactAtAPointSolution(solGuess, solOld, modelParams, solOptions)

% Form the first initial guess
yL = solGuess(1);
nxL = solGuess(2);
mL = solGuess(3);
sC = solGuess(4);
fC = solGuess(5);

% Set the derivative function
InnerDerivFun = @(x, M) RemodellingFoundationSelfContactAtPointInnerOdes(x, M, solOld, modelParams, sC, fC);
% InnerDerivFun = @(x, M) RemodellingFoundationSelfContactAtPointOdes(x, M, solOld, modelParams, sC, fC);

initInnerPoints = [0, yL, nxL, 0, 0, mL];

% Integrate the system.
innerSol = ode45(InnerDerivFun, [0, 1], initInnerPoints, solOptions);

% Evaluate the solution at the contact point
w0 = modelParams.w0; % Contact width


initOuterPoints = deval(innerSol, 1)' - [0, 0, 0, fC, 0, 0, 0];

OuterDerivFun = @(x, M) RemodellingFoundationSelfContactAtPointOuterOdes(x, M, solOld, modelParams, sC);

outerSol = ode45(OuterDerivFun, [0, 1], initOuterPoints, solOptions);

% Evaluate the end boundary conditions
innerResiduals = innerSol.y([2, 6], end) - [w0; -0.5*pi];
outerResiduals = outerSol.y([2, 3, 6], end) - [0.5; 0; 0];

Residuals = norm([innerResiduals; outerResiduals], 2);
% Residuals = norm([innerResiduals; outerResiduals; (outerSol.y(4, 1) - innerSol.y(4, end)) - fC], 2);
