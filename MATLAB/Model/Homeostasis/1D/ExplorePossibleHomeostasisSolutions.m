function ExplorePossibleHomeostasisSolutions
% Initialise the parameters.
K = 5;
nu = 0;
Es = 1;
mu = 10;
ns = -0.5;
sigma = 0.16;
W = @(s) exp(-(s./sigma).^2);

parameters.K = K; % Foundation stiffness
parameters.nu = nu; % Foundation remodelling rate
parameters.Es = Es; % Stretching stiffness
parameters.mu = mu; % Relative strength of mechanical feedback
parameters.ns = ns; % Target stress
parameters.W = W; % Wnt signal

solMesh = 0:(5*1e-3):1;

solOptions = odeset('RelTol', 1e-8,'AbsTol', 1e-8, 'Vectorized', 'On');

n0 = -0.55;
initialValues = [n0; 0; 0];

% Iterate through different values of the foundation remodelling rate rho (nu in simulations)
% and the foundation stiffness k

rhoValues = 0:10:50;
KValues = [250, 300, 500, 750, 900, 1050, 1100];

tic
for i = 1:length(rhoValues)
    parameters.nu = rhoValues(i);
    for j = 1:length(KValues)
        parameters.K = KValues(j);
        disp([parameters.nu, parameters.K])
        
        % currentResiduals = EvaluateHomeostasisSolution(initialValues, parameters, solOptions);
        
        currentSol = GetBuckledHomeostasisSolution(initialValues, solMesh, parameters, homeostaticSolShape, solOptions);
        % numSol = ShootForHomeostasisSolution(initialValues, parameters, solOptions,solTol);
        
        savedSol = [currentSol.x; currentSol.y];
        outputDirectory = '../../../../Solutions/UnbuckledRod/';
        outputValues = ['homeostasis_Es_1_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_n0_-0p55_nu_',num2str(rhoValues(i)),'_k_',num2str(KValues(j))];
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
kvals = [0, 100, 200, 250, 300, 400, 500, 600, 750, 900, 1000];
rhovals = 0:10:50;
stressvals = zeros(length(kvals), length(rhovals));

for j = 1:length(rhovals)
    rho = rhovals(j);
    for i = 1:length(kvals)
        K = kvals(i);
        outputValues = ['homeostasis_Es_1_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_n0_-0p55_nu_',num2str(rho),'_k_',num2str(K)];
        outputDirectory = '../../../../Solutions/UnbuckledRod/';
        load([outputDirectory, 'sols_', outputValues, '.mat'], 'savedSol') % Solutions
        
        stressvals(i, j) = savedSol(2, end);
    end
end

figure(1)
hold on
contourf(rhovals, kvals, stressvals)
shading interp
colorbar
colorbar('Ticks', [-0.55 100 200 300 400 500], 'FontSize', 21)
ncols = 100;
startCol = [48,128,181]./256;
endCol = [247,192,45]./256;
colormap([linspace(startCol(1), endCol(1), ncols)', ...
        linspace(startCol(2), endCol(2), ncols)', ...
        linspace(startCol(3), endCol(3), ncols)'])
% alpha(0.5)
view(2)
caxis([-0.55 500])
hold on
mu = 10; sigma = 0.16; n0 = -0.55; ns = -0.5;
rhovalline= 0:50;
kvalline = (2/(mu*sigma^2)).*(rhovalline./(1 + mu*(n0 - ns)) + 1);
line(rhovalline, kvalline, 500*ones(1, length(kvalline)), 'color', 'black', 'linewidth', 1.5)
set(gca, 'linewidth', 1.5, 'fontsize', 21)

%% 
kvals = [0, 100, 200, 250, 300, 400, 500, 600, 750, 900, 1000];
rhovals = 0:10:50;
velocityvals = zeros(length(kvals), length(rhovals));

for j = 1:length(rhovals)
    rho = rhovals(j);
    for i = 1:length(kvals)
        K = kvals(i);
        outputValues = ['homeostasis_Es_1_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_n0_-0p55_nu_',num2str(rho),'_k_',num2str(K)];
        outputDirectory = '../../../../Solutions/UnbuckledRod/';
        load([outputDirectory, 'sols_', outputValues, '.mat'], 'savedSol') % Solutions
        
        velocityvals(i, j) = savedSol(end, end);
    end
end

figure(2)
hold on
contourf(rhovals, kvals, -velocityvals)
shading interp
colorbar
colorbar('Ticks', [0 10 20 30 40], 'FontSize', 21)
ncols = 100;
startCol = [48,128,181]./256;
endCol = [247,192,45]./256;
colormap([linspace(startCol(1), endCol(1), ncols)', ...
        linspace(startCol(2), endCol(2), ncols)', ...
        linspace(startCol(3), endCol(3), ncols)'])
% alpha(0.5)
view(2)
caxis([0 40])
hold on
mu = 10; sigma = 0.16; n0 = -0.55; ns = -0.5;
rhovalline= 0:50;
kvalline = (2/(mu*sigma^2)).*(rhovalline./(1 + mu*(n0 - ns)) + 1);
line(rhovalline, kvalline, ones(1, length(kvalline)), 'color', 'black', 'linewidth', 1.5)
set(gca, 'linewidth', 1.5, 'fontsize', 21)

% kvals = [0, 100, 200, 250, 300, 400, 500, 600, 750, 900, 1000];
% rhovals = 0:10:50;
% stressvals = zeros(length(kvals), length(rhovals));
% 
% for j = 1:length(rhovals)
%     rho = rhovals(j);
%     for i = 1:length(kvals)
%         K = kvals(i);
%         outputValues = ['homeostasis_Es_1_L0_0p125_sigma_current_2w_localmechano_mu_10_ns_-0p5_n0_-0p55_nu_',num2str(rho),'_k_',num2str(K)];
%         outputDirectory = '../../../../Solutions/UnbuckledRod/';
%         load([outputDirectory, 'sols_', outputValues, '.mat'], 'savedSol') % Solutions
%         
%         stressvals(i, j) = savedSol(2, end);
%     end
% end
