function [solUnknowns, solResiduals, dtNew, numTries] = ShootForSelfContactAtAPointSolution(solGuess, solOld, modelParams, solOptions, tol, maxTries)
% Shooting method to solve for a buckled morphoelastic rod attached to a
% remodelling foundation. We know what x, theta, and ny are at S = 0, and
% so we need to guess what y, nx and m are.
% Define the objective function

dtOld = modelParams.dt;
dtNew = dtOld;

XOld = solOld.y(2,:);
YOld = solOld.y(3,:);
initPx = modelParams.Px;
initPy = modelParams.Py;
nu = modelParams.nu;
g = modelParams.g;
initGamma = modelParams.gamma;

% Evaluate the middle dt
gammaNew = initGamma + dtNew*g;
PxNew = initPx + nu*dtNew*(XOld - initPx);
PyNew = initPy + nu*dtNew*(YOld - initPy);

modelParams.gamma = gammaNew;
modelParams.Px = PxNew;
modelParams.Py = PyNew;

% Update growth and the foundations
dtNew = dtOld;
numTries = 0;
bestDt = dtOld;

OdeFunct = @(x) EvaluateSelfContactAtAPointSolution(x, solOld, modelParams, solOptions);
currentResiduals = OdeFunct(solGuess);
bestResiduals = currentResiduals;
bestGuess = solGuess;
previousResiduals = 0;

while (numTries < maxTries )
    
    if (currentResiduals > tol)
        
        numTries = numTries + 1;
        
        % Solve for the new guess
        OdeFunct = @(x) EvaluateSelfContactAtAPointSolution(x, solOld, modelParams, solOptions);
        
        currentLeftBoundaryValues = fminsearch(OdeFunct, solGuess, optimset('Display', 'Iter', 'TolX', 1e-6, 'TolFun', 1e-6));
        
        previousResiduals = currentResiduals;
        currentResiduals = OdeFunct(currentLeftBoundaryValues)
        
        if (currentResiduals < bestGuess)
            bestGuess = currentLeftBoundaryValues;
            bestResiduals = currentResiduals;
            bestDt = dtNew;
            solGuess = bestGuess;
        end
        
        % If it didn't work, we need to change the timestep, most likely
        if (currentResiduals > tol)
            
            if (currentResiduals < previousResiduals)
                
                dtNew = 1.25*dtOld;
                
            else
                
                dtNew = 0.5*dtOld;
                
            end
            dtOld = dtNew;
            
            % Re-update growth and the foundation
            gammaNew = initGamma + dtNew*g;
            PxNew = initPx + nu*dtNew*(XOld - initPx);
            PyNew = initPy + nu*dtNew*(YOld - initPy);
            
            modelParams.gamma = gammaNew;
            modelParams.Px = PxNew;
            modelParams.Py = PyNew;
            
        end
    else
        
        if (numTries < 4)
            dtNew = 1.25*dtOld;
        end
        break
    end
    
    
end
solUnknowns = bestGuess;
solResiduals = bestResiduals;
dtNew = bestDt;



