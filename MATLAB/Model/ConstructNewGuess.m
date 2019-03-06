function SolGuess = ConstructNewGuess(solCurrent, solPrevious)
% Script to construct a new solution based on a midpoint step method.

% Use the current solution mesh
SolGuess.x = solCurrent.x;

% Interpolate the previous solution onto the current solution
solPreviousInterp.x = solCurrent.x;

numSols = size(solPrevious.y, 1);

for i = 1:numSols
solPreviousInterp.y(i,:) =  interp1(solPrevious.x, solPrevious.y(i,:), solCurrent.x);
end

tangentVector = (solCurrent.y - solPreviousInterp.y)./norm((solCurrent.y - solPreviousInterp.y), 2);
SolGuess.y = solCurrent.y + tangentVector;