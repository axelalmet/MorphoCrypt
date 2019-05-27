function solUnknowns = ShootForPreContactSolution(solGuess, solOld, modelParams, solOptions, tol)
% Shooting method to solve for a buckled morphoelastic rod attached to a
% remodelling foundation. We know what x, theta, and ny are at S = 0, and
% so we need to guess what y, nx and m are. 

% Define the objective function
OdeFunct = @(x) EvaluatePreContactSolution(x, solOld, modelParams, solOptions);

currentGuess = solGuess;

currentResiduals = EvaluatePreContactSolution(currentGuess, solOld, modelParams, solOptions);

while (currentResiduals > tol)
    
    % Solve for the unknowns
    currentLeftBoundaryValues = fminsearch(OdeFunct, currentGuess);
    
    currentGuess = currentLeftBoundaryValues;
    
    % Evaluate the residuals
    currentResiduals = EvaluatePreContactSolution(currentGuess, solOld, modelParams, solOptions);

end

solUnknowns = currentGuess;
    



