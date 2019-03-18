function ExploreAgeingDistributions
% Script to explore how the age within a crypt will involve under different
% but simplistic growth laws. 

W = @(S, width) exp(-((S - 0.5*S(end))./width).^2);

sigma = 0.16;
dt = 0.05;

S = 0:(1e-2):1;
g = 1;

AOld = zeros(1, length(S));
gammaOld = 1;

%% Uniform growth. Nothing exciting, the age just linearly increases everywhere
figure(1)
hold on
plot(S, AOld)

for i = 1:100
    
    gammaNew = gammaOld + g*dt;
    
    ANew = AOld + dt - AOld.*(gammaNew - gammaOld)./gammaOld;
    
    figure(1)
    hold on
    plot(S, ANew)
    
    AOld = ANew;
    gammaOld = gammaNew;
end

figure(1)
hold on
ylim([0, 2*max(ANew)])


%% Wnt-based growth; parametrised in the current configuration
AOld = zeros(1, length(S));
gammaOld = 1;

figure(2)
hold on
plot(S, AOld)
currentS = S; % Initialise current arc length
a = 0.5;

b = zeros(1, 401);


for i = 1:400
%     
    gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
    currentS = cumtrapz(S, gammaNew);
    
    ANew = AOld + dt - 0.9*dt*W(currentS, sigma);
        
    figure(2)
    hold on
    plot(currentS, ANew)
    
    AOld = ANew;
    gammaOld = gammaNew;
end
%% Wnt-based growth; parametrised in the initial configuration
AOld = zeros(1, length(S));
gammaOld = 1;

b = zeros(1, 501);

figure(3)
hold on
plot(S, AOld)

for i = 1:500
    
    gammaNew = gammaOld.*(1 + dt*W(S, sigma));
    
    ANew = AOld + dt - AOld.*(gammaNew - gammaOld)./gammaOld;
    
    b(i + 1) = 0.525*(max(ANew) - min(ANew))./(max(ANew));
    
    figure(3)
    hold on
    plot(S, ANew)
    
    AOld = ANew;
    gammaOld = gammaNew;
end

    figure(4)
    hold on
    plot(0:500, b)
