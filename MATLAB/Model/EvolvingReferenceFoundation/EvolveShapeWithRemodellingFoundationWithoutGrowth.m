function EvolveShapeWithRemodellingFoundationWithoutGrowth
nuValues = [0.1, 1, 5, 10, 25, 50];
nuValueNames = {'0p1', '1', '5', '10', '25', '50'};
outputDirectory = '../../Solutions/RemodellingFoundation/';

for j = 1:length(nuValues)
    
% Load the solutions previously computed
outputValues = ['Eb_1_nu_', nuValueNames{j}, '_kf_0p01_L0_0p125_homoggrowth_w0_0p02'];

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
% load([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

timeSample = 1.0;
index = find(times == timeSample, 1);

% Initialise the current solution, gamma and the stresses
initSol.x = Sols{index}(1,:);
initSol.y = Sols{index}(2:(end - 2),:);

% initGamma = gammaSols{end - 1};
initGamma = parameters.gamma;
initPx = Sols{index}(end - 1, :);
initPy = Sols{index}(end, :);

parameters.uHat = Sols{index}(end, :);

dt = parameters.dt;
g = parameters.g;

% parameters.gamma = initGamma(2, 1);
% parameters.gamma = (1 + g*dt)^(index - 1);
parameters.gamma = 1 + g*timeSample;
parameters.Px = initPx;
parameters.Py = initPy;

solMesh = initSol.x;

% Set the solver options for bvp4c
maxPoints = 1e6;

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');

% Update the stresses until self-intersection or they tend to a steady state

TMax = 2.0;
dt = parameters.dt;
newTimes = timeSample:dt:TMax;
numSols = length(newTimes);

SolsWithoutGrowth = cell(numSols, 1);
% foundationSolsWithoutGrowth = cell(numSols, 1);

solOld = initSol;

SolsWithoutGrowth{1} = Sols{index};
% foundationSolsWithoutGrowth{1} = foundationSols{index};

tic

% Update the solutions in time
for i = 2:numSols
    
    % Update the solution
    [solNew, PxNew, PyNew, uHatNew] = UpdateRemodellingFoundationSolutionWithoutGrowth(solOld, parameters, solOptions);
    
    
    % self-intersects
    %     if ( ( ~isempty(InterX([solNew.y(2,:); solNew.y(3,:)])))||(norm([PxNew; PyNew] - [PxOld; PyOld], 2) < 1e-4) )
    if (~isempty(InterX([solNew.y(2,:); solNew.y(3,:)])))
        SolsWithoutGrowth = SolsWithoutGrowth(1:(i - 1));
        newTimes = newTimes(1:(i - 1));
%         foundationSolsWithoutGrowth = foundationSolsWithoutGrowth(1:(i - 1));
        
        break
    end
    
    % Stop the solution if net growth drops below unity or the curve
    parameters.Px = interp1(solOld.x, PxNew, solNew.x);
    parameters.Py = interp1(solOld.x, PyNew, solNew.x);
    
    parameters.uHat = interp1(solOld.x, uHatNew, solOld.x);
    
    solOld = solNew;
    
    SolsWithoutGrowth{i} = [solOld.x; solOld.y; parameters.Px; parameters.Py];
%     foundationSolsWithoutGrowth{i} = [parameters.Px; parameters.Py];
    
end

toc

%
% 
% SolsWithoutGrowth = SolsWithoutGrowth(1:(i - 1));
% newTimes = newTimes(1:(i - 1));
% foundationSolsWithoutGrowth = foundationSolsWithoutGrowth(1:(i - 1));

% Save the solutions
outputDirectory = '../../Solutions/RemodellingFoundation/';
outputValues = ['Eb_1_nu_', nuValueNames{j}, '_kf_0p01_L0_0p125_relaxation_homoggrowth_w0_0p02'];
save([outputDirectory, 'sols_', outputValues, '.mat'], 'SolsWithoutGrowth') % Solutions
% save([outputDirectory, 'foundationshapes_', outputValues,'.mat'], 'foundationSolsWithoutGrowth') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

end