function EvolveShapeWithUnbuckledExtensibleRodInHomeostasis
% Set the parameters
kf = 0.01; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
L = 1;
K = 12*kf*L0^4/(w*h^3);
dt = 0.01; % Time step
Es = 1; % stretching stiffness

% % Define the Wnt function
% W = @(S, width) exp(-((S - S(1))./width).^2);
W = @(S, width) exp(-((S)./width).^2);

parameters.W = W;

eta = 1.0/24; % Growth timescale (24 hours)
etaF = eta^(-1);

g = 1;

parameters.g = g;
parameters.K = 1000;
parameters.L = L; % Rod length
parameters.nu = 10*etaF*eta; % Foundation relaxation timescale
parameters.dt = dt; % Time step
parameters.Es = Es;
parameters.Ks = Es;
parameters.ns = -0.5;
parameters.mu = 10;

sigma = 2*w/L0;
parameters.sigma = sigma;

% Load the initial solutions from mathematica
solData = load('../../../../../Data/initsol_homeostaticgrowth_t_1.mat');

% Define the initial solution
solMesh = solData.Expression1(:, 1)';
sFromData = solData.Expression1(:, 2)';
nFromData = solData.Expression1(:, 3)';

solFromData.x = solMesh;
solFromData.y = [sFromData; nFromData];

% Get the initial solutions for gamma, P, g(s)
PFromData = solData.Expression1(:, 4)';
gammaFromData = solData.Expression1(:, 5)';
gFromData = solData.Expression1(:, 6)';
g = @(x) interp1(sFromData, gFromData, x); % Define incremental growth as incremental function
parameters.g = g;

% Load the sloughing boundaries at the specified times.
LMuData = load('../../../../../Data/lmu_homeostaticgrowth_t_1.mat');
times = LMuData.Expression1(:, 1)';
LMuData = LMuData.Expression1(:, 2)';

LMu = @(t) interp1(times, LMuData, t);
parameters.LMu = LMu;

TValues = [0, 1, 2, 3, 4];
TValueNames = {'0', '1', '2', '3', '4'};

for j = 1:length(TValues)
    parameters.T = TValues(j);
    
    t = times(2);
    
    parameters.t = t;
    parameters.L = LMu(t);
    %         Solve the initial bvp to obtain a structure for the first solution.
    %         Initialise growth
    firstGamma = gammaFromData.*(1 + dt*(g(sFromData)));
    parameters.gamma = firstGamma;
    
    %         Initialise foundation shape
    parameters.P = PFromData;
    
    % Set rod stiffness
    
    %         Define the ODEs and BCs
    DerivFun = @(x, M) UnbuckledExtensibleRodInHomeostasisODEs(x, M, solFromData, parameters);
    
    %         Set the boundary conditions
    BcFun = @(Ml, Mr) UnbuckledExtensibleRodHomeostasisClampedBCs(Ml, Mr, parameters);
    %         BcFun = @(Ml, Mr) UnbuckledExtensibleRodSpringBCs(Ml, Mr, parameters);
    %
    maxPoints = 1e6;
    
    tic
    %         Set the tolerances and max. number of mesh points
    solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');
    
    %         Solve the system.
    numSol = bvp4c(DerivFun, BcFun, solFromData, solOptions);
    
    toc
    
    initSol.x = solMesh;
    initSol.y = deval(numSol, solMesh);
    
    %         Update the initial data
    POld = interp1(solFromData.x, parameters.P, initSol.x);
    XOld = initSol.y(1,:);
    parameters.P = POld + parameters.dt*parameters.nu.*(XOld - POld);
    
    parameters.gamma = interp1(solFromData.x, parameters.gamma, initSol.x);
    
    parameters.Vg = cumtrapz(initSol.x, (parameters.gamma - gammaFromData)./dt);
    parameters.Vc = cumtrapz(initSol.x, (initSol.y(1,:) - solFromData.y(1,:))./dt);
    
    numSols = length(times);
    
    %         Initialise the solutions
    Sols = cell(numSols, 1);
    
    %         First non-trivial solution
    Sols{1} = [solFromData.x; solFromData.y; PFromData; gammaFromData;...
        gFromData; zeros(2, length(solFromData.x))];
    
    Sols{2} = [initSol.x; initSol.y; parameters.P; parameters.gamma;...
        (parameters.gamma - gammaFromData)./(dt*gammaFromData); parameters.Vg; parameters.Vc];
    
    %         Update the solutions in time
    
    solOld = initSol;
    
    tic
    for i = 3:numSols
        
        parameters.t = times(i);
        
        %             Update the solution
        [solNew, gammaNew, PNew, VgNew, VcNew] = UpdateUnbuckledExtensibleRodHomeostasisSolution(solOld, parameters, solOptions);
        
        %             Update the incremental growth
        gammaOld = parameters.gamma;
        gammaIncremental = (gammaNew - gammaOld)./(dt*gammaOld);
        
        %             Update gamma to the new mesh
        parameters.gamma = gammaNew;
        
        
        %             Update the foundation shape
        parameters.P = PNew;
        
        %             Update the growth velocity
        parameters.Vg = VgNew;
        
        %             Update the current velocity
        parameters.Vc = VcNew;
        
        solOld = solNew;
        
        %             Sols{i} = [solOld.x; solOld.y; gammaIncremental];
        Sols{i} = [solOld.x; solOld.y; parameters.P; parameters.gamma; gammaIncremental;...
            gammaIncremental; parameters.Vg; parameters.Vc];
        
    end
    
    toc
    
    %         Save the solutions
    outputDirectory = '../../../../Solutions/UnbuckledRod/';
    outputValues = ['Es_1_nu_10_k_1000_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_clampedbcs_homeostasis_blockedprolif_T_', TValueNames{j}];
    save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
    save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
    save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
end
%% Load the solutions
outputDirectory = '../../../../Solutions/UnbuckledRod/';
outputValues = 'Es_1_nu_10_k_1000_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_clampedbcs_homeostasis_blockedprolif_T_1';
load([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
load([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
load([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times

for i = 1:length(times)

    figure(1)
    hold on
    plot(Sols{i}(2,:), Sols{i}(end - 1, :))

    figure(2)
    hold on
    plot(Sols{i}(2,:), Sols{i}(end, :))

end
