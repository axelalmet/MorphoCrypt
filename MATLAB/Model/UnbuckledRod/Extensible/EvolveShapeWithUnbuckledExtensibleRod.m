function EvolveShapeWithUnbuckledExtensibleRod
% Set the parameters
kf = 0.01; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
L = 1;
K = 12*kf*L0^4/(w*h^3);
dt = 0.025; % Time step
Es = 1; % stretching stiffness

% Get the initial solution from AUTO
solData = load('../../../../Data/planarmorphorodsinextkf0p01L1sol_1');
% solData = load('../../../Data/seashellsK50L1sol_1');

solFromData.x = solData(:,1)';
solFromData.y = solData(:,[2, 5])';

% % Define the Wnt function
% W = @(S, width) exp(-((S - S(1))./width).^2);
W = @(S, width) exp(-((S)./width).^2);

parameters.W = W;

eta = 1.0/24; % Growth timescale (24 hours)
etaF = eta^(-1);

g = 1;

parameters.g = g;
parameters.currentArcLength = solFromData.y(1,:);
% parameters.K = K;% Foundation stiffness
% parameters.K = 1000;
parameters.L = L; % Rod length
parameters.nu = 10*etaF*eta; % Foundation relaxation timescale
parameters.dt = dt; % Time step
parameters.Es = Es;
parameters.Ks = Es;
parameters.ns = -0.4;
% parameters.mu = 10;

sigma = 2*w/L0;
parameters.sigma = sigma;

% Iterate over different parameter names
%
stiffnessValues = [1, 10, 100, 1000, 10000];
stiffnessValueNames = {'1', '10', '100', '1000', '10000'};

muValues = [0.1, 1, 10, 100];
muValueNames = {'0p1', '1', '10', '100'};
%

for j = 1:length(stiffnessValueNames)
    
    K = stiffnessValues(j);
    parameters.K = K;
    
    for k = 1:length(muValueNames)
        
        mu = muValues(k);
        parameters.mu = mu;
        
        %         for l = 1:length(stressThresholdNames)
        %
        %             ns = stressThresholds(l);
        %             parameters.ns = ns;
        
        %         Solve the initial bvp to obtain a structure for the first solution.
        SOld = solFromData.x;
        
        %         Initialise growth
        firstGamma = 1 + dt*(W(SOld, sigma));
        %     firstGamma = 1 + dt*g;
        parameters.gamma = firstGamma;
        
        solFromData.y(2,:) = (1 - firstGamma)./firstGamma;
        
        %         Initialise foundation shape
        parameters.P = SOld;
        
        % Set rod stiffness
        %     parameters.Es = 1 - b.*W(SOld, sigma);
        
        %         Define the ODEs and BCs
        DerivFun = @(x, M) UnbuckledExtensibleRodWithRemodellingFoundationOdes(x, M, solFromData, parameters);
        
        %         Set the boundary conditions
        BcFun = @(Ml, Mr) UnbuckledExtensibleRodClampedBCs(Ml, Mr, parameters);
        %         BcFun = @(Ml, Mr) UnbuckledExtensibleRodSpringBCs(Ml, Mr, parameters);
        %
        maxPoints = 1e6;
        
        tic
        %         Set the tolerances and max. number of mesh points
        solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');
        
        %         Solve the system.
        numSol = bvp4c(DerivFun, BcFun, solFromData, solOptions);
        
        toc
        
        initSol = numSol;
        
        
        %         Update the initial data
        POld = interp1(solFromData.x, parameters.P, initSol.x);
        XOld = initSol.y(1,:);
        parameters.P = POld + parameters.dt*parameters.nu.*(XOld - POld);
        
        parameters.gamma = interp1(solFromData.x, parameters.gamma, initSol.x);
        %         parameters.gamma = 1 + dt*W(XOld, sigma);
        
        % Update age
        parameters.A = dt.*ones(1, length(initSol.x));
        
        parameters.Vg = cumtrapz(initSol.x, (parameters.gamma - 1)./dt);
        %     parameters.Vg = g.*ones(1, length(initSol.x));
        parameters.Vc = cumtrapz(initSol.x, (XOld - initSol.x)./dt);
        
        TMax = 5.0;
        times = 0:dt:TMax;
        numSols = length(times);
        
        %         Initialise the solutions
        Sols = cell(numSols, 1);
        
        %         First non-trivial solution
        Sols{1} = [initSol.x; initSol.y; parameters.P; parameters.gamma.*ones(1, length(initSol.x));...
            (parameters.gamma - 1)./dt.*ones(1, length(initSol.x)); parameters.Vg; parameters.Vc; parameters.Es.*ones(1, length(initSol.x)); ...
            parameters.A.*ones(1, length(initSol.x))];
        
        %         Update the solutions in time
        
        solOld = initSol;
        
        tic
        for i = 2:numSols
            
            parameters.t = times(i);
            
            %             Update the solution
            [solNew, gammaNew, EsNew, PNew, VgNew, ANew] = UpdateUnbuckledExtensibleRodPreSloughingSolution(solOld, parameters, solOptions);
            
            %             Update the incremental growth
            gammaOld = parameters.gamma;
            gammaInc = (gammaNew - gammaOld)./(dt*gammaOld);
            gammaIncremental = interp1(solOld.x, gammaInc, solNew.x);
            %         gammaIncremental = (gammaNew - gammaOld)./(dt*gammaOld);
            
            %             Update gamma to the new mesh
            parameters.gamma = interp1(solOld.x, gammaNew, solNew.x);
            %         parameters.gamma = gammaNew;
            
            %         parameters.Es = interp1(solOld.x, EsNew, solNew.x);
            parameters.A = interp1(solOld.x, ANew, solNew.x);
            
            %             Update the foundation shape
            parameters.P = interp1(solOld.x, PNew, solNew.x);
            
            %             Update the growth velocity
            parameters.Vg = interp1(solOld.x, VgNew, solNew.x);
            %         parameters.Vg = VgNew.*ones(1, length(solNew.x));
            
            %             Update the current velocity
            parameters.Vc = (solNew.y(1,:) - interp1(solOld.x, solOld.y(1,:), solNew.x))./dt;
            
            solOld = solNew;
            
            %             Sols{i} = [solOld.x; solOld.y; gammaIncremental];
            Sols{i} = [solOld.x; solOld.y; parameters.P; parameters.gamma.*ones(1, length(solOld.x)); ...
                gammaIncremental.*ones(1, length(solOld.x)); parameters.Vg; ...
                parameters.Vc; parameters.Es.*ones(1, length(solOld.x)); parameters.A.*ones(1, length(solOld.x))];
            
        end
        
        toc
        
        %         Save the solutions
        outputDirectory = '../../../Solutions/UnbuckledRod/';
        outputValues = ['Es_1_nu_10_k_', stiffnessValueNames{j}, '_L0_0p125_current_sigma_2w_presloughing_T_5_clamped_bcs_localmechano_wntsensitivity_mu_', ...
            muValueNames{k}, '_ns_-0p4'];
        save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
        save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
        save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
        
    end
end
% end