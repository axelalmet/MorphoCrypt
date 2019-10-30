function solUnknowns = ShootForHomeostasisSolution(solGuess, modelParams, solOptions, tol)
% Shooting method to solve for a buckled morphoelastic rod attached to a
% remodelling foundation. We know that n' = 0 at s = 0 and 1 (reasonably
% expect this). 

% Define the objective function
OdeFunct = @(x) EvaluateHomeostasisSolution(x, modelParams, solOptions);

currentGuess = solGuess;

currentResiduals = EvaluateHomeostasisSolution(currentGuess, modelParams, solOptions);

while (currentResiduals > tol)
    
    % Solve for the unknowns
    currentLeftBoundaryValues = fminsearch(OdeFunct, currentGuess);
    
    currentGuess = currentLeftBoundaryValues;
    
    % Evaluate the residuals
    currentResiduals = EvaluateHomeostasisSolution(currentGuess, modelParams, solOptions);

end

solUnknowns = currentGuess;
    



