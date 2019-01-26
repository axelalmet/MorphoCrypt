function EvolveShapeWithRemodellingFoundationWithoutGrowth
% Load the solutions previously computed
outputDirectory = '../../Solutions/RemodellingFoundation/';
outputValues = 'Eb_1_nu_1_kf_0p01_L0_0p125_homoggrowth';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

timeSample = 1.0;
index = find(times == timeSample, 1);

% Initialise the current solution, gamma and the stresses
initSol.x = Sols{index}(1,:);
initSol.y = Sols{index}(2:end,:);

% initGamma = gammaSols{end - 1};
initGamma = parameters.gamma;
initP = foundationSols{index};

% parameters.gamma = initGamma(2, 1);
parameters.gamma = 1 + timeSample*parameters.g;
parameters.Px = initP(1,:);
parameters.Py = initP(2,:);

solMesh = initSol.x;

% Set the solver options for bvp4c
maxPoints = 1e6;

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');


%% Update the stresses until self-intersection or they tend to a steady state

TMax = 4.5; 
dt = parameters.dt;
newTimes = 0:dt:TMax;
numSols = length(newTimes);

SolsWithoutGrowth = cell(numSols, 1);
foundationSolsWithoutGrowth = cell(numSols, 1);

solOld = initSol;

SolsWithoutGrowth{1} = SolsWithoutGrowth{index};
foundationSolsWithoutGrowth{1} = foundationSolsWithoutGrowth{index};

tic

% Update the solutions in time
for i = 2:numSols
    
    % Update the solution
    [solNew, PxNew, PyNew] = UpdateRemodellingFoundationSolutionWithoutGrowth(solOld, parameters, solOptions);
    
    % Update the solutions, gamma, and the spring stresses
    PxOld = interp1(solOld.x, parameters.Px, solNew.x);
    PyOld = interp1(solOld.x, parameters.Py, solNew.x);

    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if ( ( ~isempty(InterX([solNew.y(2,:); solNew.y(3,:)])))||(norm([PxNew; PyNew] - [PxOld; PyOld], 2) < 1e-4) )
        
        SolsWithoutGrowth = SolsWithoutGrowth(1:(i - 1));
        newTimes = newTimes(1:(i - 1));
        foundationSolsWithoutGrowth = foundationSolsWithoutGrowth(1:(i - 1));
        
        break
    end
    
    parameters.Px = interp1(solOld.x, PxNew, solNew.x);
    parameters.Py = interp1(solOld.x, PyNew, solNew.x);
    
    solOld = solNew;
    
    SolsWithoutGrowth{i} = [solOld.x; solOld.y];
    foundationSolsWithoutGrowth{i} = [parameters.Px; parameters.Py];
        
end

toc
        
% Save the solutions
outputDirectory = '../../Solutions/RemodellingFoundation/';
outputValues = 'Eb_1_nu_1_kf_0p01_L0_0p125_relaxation';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'SolsWithoutGrowth') % Solutions
save([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationSolsWithoutGrowth') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
