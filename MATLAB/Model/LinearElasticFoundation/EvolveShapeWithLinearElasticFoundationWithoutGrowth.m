function EvolveShapeWithLinearElasticFoundationWithoutGrowth
% Load the solutions previously computed
outputDirectory = '../../Solutions/LinearElasticFoundation/';
outputValues = 'Eb_1_kf_0p01_L0_0p125_homoggrowth_chi_0p1';

load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
% load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

timeSample = 1.0;
index = find(times == timeSample, 1);

% Set the parameters
kf = 0.01; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod`
L = 1;
y0 = 0;
K = 12*kf*L0^4/(w*h^3);
Eb = 1; % Bending stiffness
dt = 1e-2;

eta = 1.0/24; % Growth timescale (24 hours)
etaK = eta^(-1);

g = 1;

parameters.g = g;
parameters.K = K;% Foundation stiffness
parameters.L = 0.5*L; % Rod length
parameters.y0 = y0; % Rod centreline
parameters.Eb = Eb; % Bending stiffness
parameters.ext = 0; % Exstensibility
parameters.chi = 0.1*etaK*eta; % Curvature relaxation timescale
parameters.dt = dt; % Time step

% Initialise the current solution, gamma and the stresses
initSol.x = Sols{index}(1,:);
initSol.y = Sols{index}(2:(end - 1),:);

% parameters.gamma = initGamma(2, 1);
parameters.gamma = 1 + timeSample*parameters.g;
parameters.uHat = Sols{end}(end, :);


% Set the solver options for bvp4c
maxPoints = 1e6;

% Set the tolerances and max. number of mesh points
solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');

%
%% Evolve the intrinsic curvature

TMax = 10.0;
dt = parameters.dt;
newTimes = 0:dt:TMax;
numSols = length(newTimes);

SolsWithoutGrowth = cell(numSols, 1);

solOld = initSol;

SolsWithoutGrowth{1} = Sols{index};

tic

% Update the solutions in time
for i = 2:numSols
    
    % Update the solution
    [solNew, uHatNew] = UpdateLinearElasticSolutionWithoutGrowth(solOld, parameters, solOptions);
    
    % Update the solutions and intrinsic curvature
    uHatOld = parameters.uHat;
    
    % Stop the solution if net growth drops below unity or the curve
    % self-intersects
    if ( ( ~isempty(InterX([solNew.y(2,:); solNew.y(3,:)])))||(norm(uHatNew - uHatOld, 2) < 1e-4) )
        
        SolsWithoutGrowth = SolsWithoutGrowth(1:(i - 1));
        newTimes = newTimes(1:(i - 1));
        
        break
    end
    
    parameters.uHat = interp1(solOld.x, uHatNew, solNew.x);
    
    solOld = solNew;
    
    SolsWithoutGrowth{i} = [solOld.x; solOld.y; parameters.uHat];
    
end

toc

%%
SolsWithoutGrowth = SolsWithoutGrowth(1:(i - 1));
newTimes = newTimes(1:(i - 1));

% Save the solutions
outputDirectory = '../../Solutions/LinearElasticFoundation/';
outputValues = 'Eb_1_kf_0p01_L0_0p125_chi_0p1_relaxation';
save([outputDirectory, 'sols_', outputValues, '.mat'], 'SolsWithoutGrowth') % Solutions
save([outputDirectory, 'times_', outputValues, '.mat'], 'newTimes') % Times
