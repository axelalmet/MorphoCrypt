function Sol = GetSelfContactAtAPointSolution(initValues, solMesh, solOld, modelParams, solOptions)

% Split the solution mesh
meshSplitIndex = find(solMesh == 1, 1);
innerMesh = solMesh(1:meshSplitIndex);
outerMesh = solMesh((meshSplitIndex + 1):end) - 1;

% Form the first initial guess
yL = initValues(1);
nxL = initValues(2);
mL = initValues(3);
sC = initValues(4);
fC = initValues(5);

% Set the derivative function
InnerDerivFun = @(x, M) RemodellingFoundationSelfContactAtPointInnerOdes(x, M, solOld, modelParams, sC, fC);

initInnerPoints = [0, 0, yL, nxL, 0, 0, mL];

% Integrate the system.
innerSol = ode45(InnerDerivFun, [0, 1], initInnerPoints, solOptions);

% % Now integrate past the contact point
% initOuterPoints = innerSol.y(:, end)' - [0, 0, 0, fC, 0, 0, 0];
% 
% OuterDerivFun = @(x, M) RemodellingFoundationSelfContactAtPointOuterOdes(x, M, solOld, modelParams, sC);
% 
% outerSol = ode45(OuterDerivFun, [0, 1], initOuterPoints, solOptions);
% 
Sol.x = solMesh;
Sol.y = [deval(innerSol, innerMesh)];