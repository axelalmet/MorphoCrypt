function ExplorePossibleBuckledHomeostasisSolutions
% Initialise the solutions and parameters.
% clear all

outputDirectory = '../../../Solutions/RemodellingFoundation/';
outputValues = 'Eb_1_nu_10_kf_0p01_L0_0p125_current_sigma_2w_w0_0';

Sols = load([outputDirectory, 'sols_',outputValues,'.mat'], 'Sols');
Sols = Sols.Sols;
parameters = load([outputDirectory, 'parameters_',outputValues,'.mat'], 'parameters');
parameters = parameters.parameters;

% %%
% 
% 
% K = 5;
% nu = 0;
% Es = 1;
% mu = 10;
% ns = -0.5;
% sigma = 0.16;
% W = @(s) exp(-(s./sigma).^2);

% parameters.K = K; % Foundation stiffness
% parameters.nu = nu; % Foundation remodelling rate
% parameters.Es = Es; % Stretching stiffness
% parameters.mu = mu; % Relative strength of mechanical feedback
% parameters.ns = ns; % Target stress
% parameters.W = W; % Wnt signal

solOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Vectorized', 'On');

% Set the relevant parameters
parameters.mu = 0.5;
parameters.ns = -15;

n0 = -20;

parameters.n0 = n0;

parameters.Es = (15)^2/(12*125^2);

% Take all the solutions at s = epsilon, rather than s = 0, as this causes
% a singularity in the numerics.
currentArcLengthOld = parameters.currentArcLength;

solMesh = currentArcLengthOld(2:end);

parameters.currentArcLength = solMesh;

PxOld = parameters.Px;
PyOld = parameters.Py;

parameters.Px = PxOld(2:end);
parameters.Py = PyOld(2:end);

homeostaticSol.x = solMesh;
homeostaticSol.y = Sols{end}(2:end,2:end);

% Should give us small, non-zero starting value
epsilon = solMesh(1);

% Initialise the values based on Taylor expansions about s = epsilon
initialValues = [n0; 0; ...
    epsilon*parameters.nu/(parameters.nu + 1 + parameters.mu*(n0 - parameters.ns)); ...
    parameters.Py(2); -epsilon*(1 + parameters.mu*tanh((n0 - parameters.ns)))];

% Define the slenderness

%%

% Iterate through different values of the foundation remodelling rate rho (nu in simulations)
% and the foundation stiffness k

KValues = [0.1 1 10 100];
KValueNames = {'0p1', '1', '10', '100'};
rhoValues = 0:10;

tic
for i = 1:length(rhoValues)
    rho = rhoValues(i);
    parameters.nu = rho;
    for j = 1:length(KValues)
        K = KValues(j);
        parameters.K = K;
        disp([parameters.nu, parameters.K])
        
        
        currentSol = GetBuckledExtensibleHomeostasisSolution(initialValues, solMesh, parameters, homeostaticSol, solOptions);
        
        savedSol = [0 currentSol.x; [n0; 0; 0; homeostaticSol.y(3,1); 0], currentSol.y; Sols{end}([3, 4, 7],:)];
        outputDirectory = '../../../Solutions/BuckledRod/';
        outputValues = ['homeostasis_ext_Es_0p0012_L0_0p125_sigma_current_2w_localmechano_mu_0p5_ns_-15_n0_-20_nu_',num2str(rho),'_k_',num2str(KValueNames{j})];
        save([outputDirectory, 'sols_', outputValues, '.mat'], 'savedSol') % Solutions
    end
end
toc


% %% Plot the different growth profiles
% rhovals = [0, 25, 50];
% K = 800;
% 
% for i = 1:length(rhovals)
%     rho = rhovals(i);
%     outputValues = ['homeostasis_Es_1_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_n0_-0p55_nu_',num2str(rho),'_k_',num2str(K)];
%     outputDirectory = '../../../../Solutions/UnbuckledRod/';
%     load([outputDirectory, 'sols_', outputValues, '.mat'], 'savedSol') % Solutions
%     
%     
%     % %% Plot the resulting growth profile
%     W = @(s) exp(-(s./0.16).^2);
%     
%     figure(1)
%     hold on
%     plot(savedSol(1,:), W(savedSol(1,:)) +10.*(savedSol(2,:) + 0.5).*(savedSol(2,:) < -0.5))
%     plot(-savedSol(1,:), W(savedSol(1,:)) +10.*(savedSol(2,:) + 0.5).*(savedSol(2,:) < -0.5))
%     
%     figure(2)
%     hold on
%     plot(savedSol(1,:), savedSol(2,:))
%     %     ylim([-0.6 1])
%     xlim([0 1])
%     
%     figure(3)
%     hold on
%     plot(savedSol(1,:), -savedSol(3,:))
%     
% end

% %% Plot the different growth profiles
% rho = 35;
% Kvals = 0:200:1000;
% 
% for i = 1:length(Kvals)
%     K = Kvals(i);
%     outputValues = ['homeostasis_Es_1_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_n0_-0p55_nu_',num2str(rho),'_k_',num2str(K)];
%     outputDirectory = '../../../../Solutions/UnbuckledRod/';
%     load([outputDirectory, 'sols_', outputValues, '.mat'], 'savedSol') % Solutions
%     
%     
%     % %% Plot the resulting growth profile
%     W = @(s) exp(-(s./0.16).^2);
%     
%     figure(1)
%     hold on
%     plot(savedSol(1,:), W(savedSol(1,:)) +10.*(savedSol(2,:) + 0.5).*(savedSol(2,:) < -0.5))
%     plot(-savedSol(1,:), W(savedSol(1,:)) +10.*(savedSol(2,:) + 0.5).*(savedSol(2,:) < -0.5))
%     
%     figure(2)
%     hold on
%     plot(savedSol(1,:), savedSol(2,:))
%     ylim([-0.6 2])
%     
%     figure(3)
%     hold on
%     plot(savedSol(1,:), savedSol(4,:))
%     ylim([-5 0.5])
%     
% end
% %
% % %% Plot the phase space that we predicted
% %
% % rhoValues = 0:50;
% % KValues = (2/(mu*sigma^2)).*(rhoValues./(1 + mu*(n0 - ns)) + 1);
% % plot(rhoValues, KValues)
% % ylim([0 1000])
% 
%%
% kvals = [0, 10, 100, 1000];
% rhovals = [0, 0.1, 1, 10, 50];
% stressvals = zeros(length(kvals), length(rhovals));
% 
% for j = 1:length(rhovals)
%     rho = rhovals(j);
%     for i = 1:length(kvals)
%         K = kvals(i);
%         outputValues = ['homeostasis_Es_0p0012_L0_0p125_sigma_current_2w_localmechano_mu_0p6_ns_-15_n0_-17p5_nu_',num2str(rho),'_k_',num2str(K)];
%         outputDirectory = '../../../../Solutions/UnbuckledRod/';
%         load([outputDirectory, 'sols_', outputValues, '.mat'], 'savedSol') % Solutions
%         
%         stressvals(i, j) = savedSol(2, end);
% %     end
% end
% 
% figure(1)
% hold on
% contourf(rhovals, kvals, stressvals)
% shading interp
% colorbar
% colorbar('Ticks', [-0.55 100 200 300 400 500], 'FontSize', 21)
% ncols = 100;
% startCol = [48,128,181]./256;
% endCol = [247,192,45]./256;
% colormap([linspace(startCol(1), endCol(1), ncols)', ...
%         linspace(startCol(2), endCol(2), ncols)', ...
%         linspace(startCol(3), endCol(3), ncols)'])
% % alpha(0.5)
% view(2)
% caxis([-0.55 500])
% hold on
% mu = 10; sigma = 0.16; n0 = -0.55; ns = -0.5;
% rhovalline= 0:50;
% kvalline = (2/(mu*sigma^2)).*(rhovalline./(1 + mu*(n0 - ns)) + 1);
% line(rhovalline, kvalline, 500*ones(1, length(kvalline)), 'color', 'black', 'linewidth', 1.5)
% set(gca, 'linewidth', 1.5, 'fontsize', 21)
% 
% %% 
% kvals = [0, 100, 200, 250, 300, 400, 500, 600, 750, 900, 1000];
% rhovals = 0:10:50;
% velocityvals = zeros(length(kvals), length(rhovals));
% 
% for j = 1:length(rhovals)
%     rho = rhovals(j);
%     for i = 1:length(kvals)
%         K = kvals(i);
%         outputValues = ['homeostasis_Es_1_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_n0_-0p55_nu_',num2str(rho),'_k_',num2str(K)];
%         outputDirectory = '../../../../Solutions/UnbuckledRod/';
%         load([outputDirectory, 'sols_', outputValues, '.mat'], 'savedSol') % Solutions
%         
%         velocityvals(i, j) = savedSol(end, end);
%     end
% end
% 
% figure(2)
% hold on
% contourf(rhovals, kvals, -velocityvals)
% shading interp
% colorbar
% colorbar('Ticks', [0 10 20 30 40], 'FontSize', 21)
% ncols = 100;
% startCol = [48,128,181]./256;
% endCol = [247,192,45]./256;
% colormap([linspace(startCol(1), endCol(1), ncols)', ...
%         linspace(startCol(2), endCol(2), ncols)', ...
%         linspace(startCol(3), endCol(3), ncols)'])
% % alpha(0.5)
% view(2)
% caxis([0 40])
% hold on
% mu = 10; sigma = 0.16; n0 = -0.55; ns = -0.5;
% rhovalline= 0:50;
% kvalline = (2/(mu*sigma^2)).*(rhovalline./(1 + mu*(n0 - ns)) + 1);
% line(rhovalline, kvalline, ones(1, length(kvalline)), 'color', 'black', 'linewidth', 1.5)
% set(gca, 'linewidth', 1.5, 'fontsize', 21)
% 
% % kvals = [0, 100, 200, 250, 300, 400, 500, 600, 750, 900, 1000];
% % rhovals = 0:10:50;
% % stressvals = zeros(length(kvals), length(rhovals));
% % 
% % for j = 1:length(rhovals)
% %     rho = rhovals(j);
% %     for i = 1:length(kvals)
% %         K = kvals(i);
% %         outputValues = ['homeostasis_Es_1_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_n0_-0p55_nu_',num2str(rho),'_k_',num2str(K)];
% %         outputDirectory = '../../../../Solutions/UnbuckledRod/';
% %         load([outputDirectory, 'sols_', outputValues, '.mat'], 'savedSol') % Solutions
% %         
% %         stressvals(i, j) = savedSol(2, end);
% %     end
% % end
