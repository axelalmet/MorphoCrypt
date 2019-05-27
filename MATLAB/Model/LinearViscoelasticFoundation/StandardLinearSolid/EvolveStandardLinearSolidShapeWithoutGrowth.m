function EvolveStandardLinearSolidShapeWithoutGrowth
% Load the solutions previously computed
outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/StandardLinearSolid/';
outputValues = 'Eb_1_beta_0p75_nu_10_kf_0p01_L0_0p125_homoggrowth';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
% load([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
load([outputDirectory, 'stresses_', outputValues,'.mat'], 'stressSols') % Foundation stresses
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

timeSample = 1.0;
index = find(times == timeSample, 1);

% Initialise the current solution, gamma and the stresses
initSol.x = Sols{index}(1,:);
initSol.y = Sols{index}(2:end,:);

% initGamma = gammaSols{end - 1};
initGamma = parameters.gamma;
initP = stressSols{index};

% parameters.gamma = initGamma(2, 1);
parameters.gamma = 1 + timeSample*parameters.g;
parameters.P = initP(2,:);

solMesh = initSol.x;

% Set the solver options for bvp4c
maxPoints = 1e6;

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');
%% Update the stresses until self-intersection or they tend to a steady state

TMax = 4.5; 
dt = parameters.dt;
L = parameters.L;
newTimes = 0:dt:TMax;
numSols = length(newTimes);

SolsWithoutGrowth = cell(numSols, 1);
stressSolsWithoutGrowth = cell(numSols, 1);

solOld = initSol;

SolsWithoutGrowth{1} = Sols{index};
stressSolsWithoutGrowth{1} = stressSols{index};

%%

tic

for i = 77:numSols
        
    % Update the solution
    [solMeshNew, solNew, PNew] = UpdateStandardLinearSolidSolutionWithoutGrowth(solMesh, solOld, parameters, solOptions);
    
    % Stop the solution the curve self-intersects or the stress tends to a
    % steady state
    if ( (HasRodHitSelfContact(solNew, parameters))||(norm(PNew - parameters.P, 2) < 1e-4) )
        
        SolsWithoutGrowth = SolsWithoutGrowth(1:(i - 1));
        newTimes = newTimes(1:(i - 1));
        stressSolsWithoutGrowth = stressSolsWithoutGrowth(1:(i - 1));
        
        break
    end
    
    % Update the solutions and the spring stresses
    parameters.P = interp1(solOld.x, PNew, solNew.x);
    solOld = solNew;

    SolsWithoutGrowth{i} = [solOld.x; solOld.y];
    stressSolsWithoutGrowth{i} = [L.*solNew.x; parameters.P];
    
    solMesh = solMeshNew;

                        
end

toc

%%
SolsWithoutGrowth = SolsWithoutGrowth(1:(i - 1));
newTimes = newTimes(1:(i - 1));
stressSolsWithoutGrowth = stressSolsWithoutGrowth(1:(i - 1));
        
outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/StandardLinearSolid/';
outputValues = 'Eb_1_beta_0p75_nu_10_kf_0p01_L0_0p125_relaxation';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'SolsWithoutGrowth') % Solutions
save([outputDirectory, 'stresses_', outputValues,'.mat'], 'stressSolsWithoutGrowth') % Foundation stresses
save([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
