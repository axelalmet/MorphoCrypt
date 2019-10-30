function EvolveShapeWithStandardLinearSolidFoundation
% Set the parameters
kf = 0.01; % Dimensional foundational stiffness
h = 0.015; % Thickness of the rod cross section
w = 0.01; % Width of the rod cross section
L0 = 0.125; % Dimensional length of the rod
% L = 2*sqrt(3)*L0/h; %Dimensionless length
L = 1;
y0  = 6*L;
% K = kf*h/(12*w); % Dimensionless foundation stiffness
K = 12*kf*L0^4/(w*h^3);
l = 12*L0^4/(w*h^3);
n3s = 0; % Target axial tension
Es = 1; % Stretching stiffness
b1 = 0.0; % Bending stiffness
dt = 1e-2; % Time step

% Get the initial solution from AUTO
solData = load('../../../../Data/planarmorphorodsinextkf0p01L1sol_1'); %
% Cartesian basis

solFromData.x = solData(:,1)';
solFromData.y = solData(:,2:end)';

% Translate solution up
solFromData.y(3,:) = y0 + solFromData.y(3,:);

%
sigma = 2*sqrt(3)*2*w/h; % "Width" of Wnt gradient
% % Define the Wnt function
W = @(S, width) exp(-((S - 0.5*L)./width).^2);

eta = 24; % Define eta such that \eta^{-1} = 24 hours
mu = 0;
beta = 1;

g = 1;

parameters.w0 = 0;

parameters.beta = beta;
parameters.g = g;
parameters.K = K;
parameters.L = 0.5*L; % Rod length
parameters.l = l;
parameters.y0 = y0; % Rod centreline
parameters.sigma = sigma; % Width of wnt gradient
parameters.mu = mu; % Rate of mechanical inhibition
parameters.eta = eta; % Rate of chemical change
parameters.n3s = n3s; % Target axial stress
parameters.Es = Es; % Stretch stiffness
parameters.b1 = b1;
parameters.Eb = 1; % Bending stiffness
parameters.ext = 0; % Exstensibility
% parameters.nu = etaV/(kf*eta); % Foundation relaxation timescale
parameters.etaK = 0; % Curvature relaxation timescale
parameters.dt = dt; % Time step

% Shift the solutions so that we only focus on the 2nd half.
solFromData.x = solFromData.x(300:end) - 0.5;
solFromData.x(1) = 0;
solFromData.x = solFromData.x/(solFromData.x(end));
solFromData.y = solFromData.y(:, 300:end);
solFromData.y(4,:) = solData(1:302,5)';
solFromData.y(1, :) = solFromData.y(1, :) - 0.5;
solFromData.y(2, :) = solFromData.y(2, :) - 0.5;

solFromData.y([1; 2], 1) = [0; 0];

% Solve the initial bvp to obtain a structure for the first solution.
initS = solFromData.y(1,:);
initX = solFromData.y(2,:);
initY = solFromData.y(3,:);

initDelta = sqrt((initX - initS).^2 + (initY).^2);

nuValues = [0.01, 0.05, 0.1, 0.5, 1];
nuValueNames = {'0p01', '0p05', '0p1', '0p5', '1'};

parameters.uHat = zeros(1, length(solFromData.x));

for j = 1:length(nuValues)
    
    parameters.nu = nuValues(j);
    % parameters.nu = nuValues(1);
    gammaOld = 1;
    firstGamma = gammaOld + g*dt;
    parameters.gamma = firstGamma;
    
    parameters.P = zeros(1, length(solFromData.x));
    
    % Define the ODEs and BCs - all first steps are from a linearly elastic
    % foundation
    DerivFun = @(x, M) StandardLinearSolidFoundationOdes(x, M, solFromData, parameters);
    
    % Set the boundary conditions
    % BcFun = @(Ml, Mr) NonUniformGrowthBCs(Ml, Mr, parameters);
    BcFun = @(Ml, Mr) PreSloughingBCs(Ml, Mr, parameters);
    
    maxPoints = 1e6;
    
    tic
    % Set the tolerances and max. number of mesh points
    solOptions = bvpset('RelTol', 1e-4,'AbsTol', 1e-4, 'NMax', maxPoints, 'Vectorized', 'On');
    
    % Solve the system.
    numSol = bvp4c(DerivFun, BcFun, solFromData, solOptions);
    
    toc
    
    % initSol.x = solFromData.x;
    % initSol.y = deval(numSol, solFromData.x);
    
    initSol = numSol;
    
    % gammaOld = interp1(solFromData.x, firstGamma, initSol.x);
    % parameters.gamma = gammaOld;
    
    %%
    solOld = initSol;
    
    initS = initSol.y(1,:);
    initX = initSol.y(2,:);
    initY = initSol.y(3,:);
    
    initDelta = sqrt((initX - initS).^2 + initY.^2);
    
    parameters.P = dt*beta*(1 - beta)*K*parameters.nu.*(initDelta - y0) ...
        + (K).*(initDelta - y0);
    
    solMesh = initSol.x;
    
    
    % Set the times we want to solve the problem for
    
    TMax = 1000.0;
    times = 0:dt:TMax;
    numSols = length(times);
    
    % Initialise the solutions
    Sols = cell(numSols, 1);
    % gammaSols = cell(numSols, 1);
    stressSols = cell(numSols, 1);
    
    %  The first solution is always flat
    flatSol.x = initSol.x;
    flatSol.y = initSol.y;
    
    flatSol.y(3,:) = 0.*flatSol.y(3,:) + y0;
    flatSol.y(5:end,:) = 0.*flatSol.y(5:end,:);
    
    Sols{1} = [flatSol.x; flatSol.y];
    % gammaSols{1} = [L.*flatSol.x; ones(1, length(flatSol.x))];
    stressSols{1} = [flatSol.x; 0.*flatSol.y(1, :); 0*flatSol.x;];
    
    
    % First non-trivial solution
    Sols{2} = [solOld.x; solOld.y; parameters.P];
    % gammaSols{2} = [L.*initSol.x; gammaOld];
    % stressSols{2} = [L.*initSol.x; parameters.P];
    
    tic
    
    % Update the solutions in time
    for i = 3:numSols
        
        % Update the solution
        [solMeshNew, solNew, gammaNew, EbNew, PNew] = UpdateStandardLinearSolidSolution(solMesh, solOld, W, parameters, solOptions);
        
        % Update the solutions, gamma, and the spring stresses
        %     gammaOld = interp1(solOld.x, gammaNew, solNew.x);
        gammaOld = gammaNew;
        parameters.gamma = gammaOld;
        parameters.P = interp1(solOld.x, PNew, solNew.x);
        %     parameters.Eb = interp1(solOld.x, EbNew, solNew.x);
        
        % Stop the solution if net growth drops below unity or the curve
        % self-intersects
        %     if ( (trapz(solNew.x, gammaOld) < 1)||(~isempty(InterX([solNew.y(2,:); solNew.y(3,:)]))) )
        if (HasRodHitSelfContact(solNew, parameters)||(times(i) > 2.0))
            
            Sols = Sols(1:(i - 1));
            %    gammaSols = gammaSols(1:(i - 1));
            times = times(1:(i - 1));
            stressSols = stressSols(1:(i - 1));
            
            break
        end
        
        solOld = solNew;
        
        Sols{i} = [solOld.x; solOld.y; parameters.P];
        %     gammaSols{i} = [L.*solNew.x; gammaOld];
        %     stressSols{i} = [L.*solNew.x; parameters.P];
        
        solMesh = solMeshNew;
        
    end
    
    toc
    
    % %%
    % Sols = Sols(1:21);
    % %    gammaSols = gammaSols(1:(i - 1));
    % times = times(1:21);
    % stressSols = stressSols(1:21);
    
    % Save the solutions
    outputDirectory = '../../../Solutions/LinearViscoelasticFoundation/StandardLinearSolid/';
    outputValues = ['Eb_1_beta_1_nu_', nuValueNames{j}, '_kf_0p01_L0_0p125_homoggrowth'];
    save([outputDirectory, 'sols_', outputValues, '.mat'], 'Sols') % Solutions
    % save([outputDirectory, 'gamma_', outputValues,'.mat'], 'gammaSols') % Gamma
    save([outputDirectory, 'stresses_', outputValues,'.mat'], 'stressSols') % Foundation stresses
    save([outputDirectory, 'times_', outputValues, '.mat'], 'times') % Times
    save([outputDirectory, 'parameters_', outputValues, '.mat'], 'parameters') % Times
    
end

