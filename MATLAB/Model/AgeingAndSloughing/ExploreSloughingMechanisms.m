function ExploreSloughingMechanisms
% Script to explore how the sloughing region changes under different
% 'sloughing rate' mechanisms
% 
% clear all
% %% Parameter set-up
% L0 = 1;
% L = L0;
% g = 1;
% dt = 0.05;
% T = 3.0;
% tMax = 3.0;
% AStar = 5;
% sigma = 10;
% times = 0:dt:tMax;
% 
% %% Homogeneous, time-dependent sloughing
% 
% LOld = 1;
% newLengths = zeros(1, length(times));
% newLengths(1) = 1;
% AOld = 0;
% gammaOld = 1;
% g = 1;
% sigma = 0.16;
% 
% W = @(S, width) exp(-((S)./width).^2);
% 
% S = 0:1e-3:1;
% currentS = S;
% 
% options = optimset('Display','off');
% 
% % Iterate using derived numerical scheme for derivatives
% for i = 2:(length(times))
%     
%     gammaNew = gammaOld + g*dt;
%     
%     ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
%     
%     t = times(i);
%     
%     sloughedAmount = trapz(0:1e-2:t, 0.5*(1 + tanh(10*((0:1e-2:t) - T))) );
%     
%     currentArcLength = gammaNew.*S;
%     
%     if (sloughedAmount > 1e-4)
%         % Define new functionx
%         sloughLfunct = @(l) sloughedAmount - (currentArcLength(end) - interp1(S, currentArcLength, l));
%         
%         sloughL = min(fsolve(sloughLfunct, [0 1], options));
%         newLengths(i) = sloughL;
%         
%     else
%         newLengths(i) = 1;
%     end
%     
%     gammaOld = gammaNew;
%     AOld = ANew;
%     
% end
% 
% figure
% plot(times, newLengths)
% 
% %% Chemical signalling, time-dependent sloughing
% LOld = 1;
% newLengths = zeros(1, length(times));
% newLengths(1) = 1;
% AOld = 0;
% gammaOld = 1;
% sigma = 0.16;
% 
% W = @(S, width) exp(-((S)./width).^2);
% 
% S = 0:1e-3:1;
% currentS = S;
% 
% options = optimset('Display','off');
% 
% % Iterate using derived numerical scheme for derivatives
% for i = 2:(length(times))
%     
%     gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
%     
%     ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
%     
%     t = times(i);
%     
%     netGrowth = trapz(S, gammaOld.*W(currentS, sigma));
%     
%     sloughedAmount = trapz(0:1e-2:t, 0.5*netGrowth*(1 + tanh((10*((0:1e-2:t) - T)))));
%     
%     currentS = cumtrapz(S, gammaNew);
%     
%     if (sloughedAmount > 1e-4)
%         % Define new functionx
%         sloughLfunct = @(l) sloughedAmount - (currentS(end) - interp1(S, currentS, l));
%         
%         sloughL = min(fsolve(sloughLfunct, [0 1], options));
%         newLengths(i) = sloughL;
%         
%     else
%         newLengths(i) = 1;
%     end
%     
%     gammaOld = gammaNew;
%     AOld = ANew;
%     
% end
% 
% %%
% figure
% plot(times, newLengths, '-^')
% 
% %% Chemical signalling, age-dependent sloughing
% LOld = 1;
% newLengths = zeros(1, length(times));
% newLengths(1) = LOld;
% AOld = 0;
% gammaOld = 1;
% sigma = 0.16;
% 
% W = @(S, width) exp(-((S)./width).^2);
% 
% S = 0:1e-3:LOld;
% currentS = S;
% 
% sloughedAmount = zeros(1, length(times));
% 
% options = optimset('Display','off');
% 
% % Iterate using derived numerical scheme for derivatives
% for i = 2:(length(times))
%     
%     gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
%     
%     ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
%     
%     sloughedAmount(i) = trapz(S, ones(1, length(S)).*(ANew > AStar).*(AOld.*W(currentS, sigma) < 0.01));
%     
%     gammaOld = gammaNew;
%     AOld = ANew;
%     
%     if (sloughedAmount(i) > 1e-4)
%         
%         break
%         
%     else
%         newLengths(i) = 1;
%     end
%     
%     figure(1)
%     hold on
%     plot(S, ANew)
%     
% end
% 
% %%
% 
% % Get the new sloughing boundary
% 
% % Define new functionx
% sloughLfunct = @(l) sloughedAmount(i) - (currentS(end) - interp1(S, currentS, l));
% 
% sloughL = fsolve(sloughLfunct, 0, options);
% newLengths(i) = sloughL;
% 
% %%
% 
% L = sloughL;
% 
% SOld = S;
% gammaPrev = gammaOld;
% APrev = AOld;
% 
% S = [SOld(1:(end - 1))./(SOld(end - 1)).*L, L, 1];
% 
% gammaOld = [gammaPrev(1:(end - 1)), gammaPrev((end - 1):end)];
% currentS = cumtrapz(S, gammaOld);
% 
% AOld = [APrev(1:(end - 1)), APrev((end - 1):end)];
% 
% LOld = L;
% 
% %%
% 
% % Iterate now with the evolving sloughing boundary
% for j = (i + 1):length(times)
%     gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
%     
%     ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
%     
%     LIndex = find(S == LOld, 1);
%     
%     sloughedAmount(j) = trapz(S(1:LIndex), ones(1, length(S(1:LIndex))).* ...
%         (ANew(1:LIndex) > AStar).*(AOld(1:LIndex).*W(currentS(1:LIndex), sigma) < 0.01));
%     
%     % Define function to determine new boundary
%     sloughLfunct = @(l) sloughedAmount(j) - (currentS(end) ...
%         - interp1(S([1:LIndex, (LIndex + 2):end]), ...
%         currentS([1:LIndex, (LIndex + 2):end]), l));
%     
%     sloughL = fsolve(sloughLfunct, 0, options);
%     newLengths(j) = sloughL;
%     
%     S = [SOld(1:(LIndex))./(SOld(LIndex))*sloughL, sloughL, 1];
%     
%     LOld = sloughL;
%     
%     gammaOld = gammaNew;
%     AOld = ANew;
%     
%     currentS = cumtrapz(S, gammaOld);
%     
%     figure(1)
%     hold on
%     plot(S, ANew)
%     
% end
% 
% figure(2)
% hold on
% plot(times, newLengths)
% plot(times, sloughedAmount)
% %% Chemical signal-driven growth, time-dependent sloughing
% % We're going to work out how much needs to be sloughed at each time to
% % balance growth.
% 
% LOld = 1;
% newLengths = zeros(1, length(times));
% newLengths(1) = LOld;
% AOld = 0;
% gammaOld = 1;
% sigma = 0.16;
% 
% W = @(S, width) exp(-((S)./width).^2);
% 
% S = 0:1e-3:LOld;
% currentS = S;
% 
% sloughedAmount = zeros(1, length(times));
% 
% options = optimset('Display','off');
% 
% % Iterate using derived numerical scheme for derivatives
% for i = 2:(length(times))
%     
%     gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
%     
%     %     ANew = 2*AOld + dt - gammaNew./gammaOld.*AOld;
%     
%     currentS = cumtrapz(S, gammaNew);
%     
%     gammaOld = gammaNew;
%     %     AOld = ANew;
%     
%     if (times(i) >= T)
%         
%         break
%         
%     else
%         newLengths(i) = 1;
%     end
%     
%     %     figure(1)
%     %     hold on
%     %     plot(S, ANew)
%     
% end
% 
% %%
% 
% count = i - 1;
% 
% sloughingRateMagnitudes = zeros(1, length(times((count + 1):end)));
% amountGrownInNonSloughedRegion = zeros(1, length(times((count + 1):end)));
% 
% while(newLengths(count) > 0)
%     
%     count = count + 1;
%     
%     gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
%     
%     currentS = cumtrapz(S, gammaNew);
%     
%     gammaOld = gammaNew;
%     
%     sloughedAmountScaled = trapz(times(1:count), 0.5*(1 + tanh(10*(times(1:count) - T))));
%     
%     newLengths(count) = newLengths(count - 1) - 0.005;
%     
%     amountGrown = currentS(end) - interp1(S, currentS, newLengths(count));
%     
%     sloughingRateMagnitudes(count - i + 1) = amountGrown./sloughedAmountScaled;
%     
%     amountGrownInNonSloughedRegion(count - i + 1) = currentS(end) - interp1(S, currentS, newLengths(count))./(1 - newLengths(count));
% end
% 
% %% Explore sloughing mechanisms
% 
% LOld = 1;
% newLengths = zeros(1, length(times));
% AOld = 0;
% gammaOld = 1;
% sigma = 0.16;
% 
% W = @(S, width) exp(-((S)./width).^2);
% 
% S = 0:1e-3:LOld;
% currentS = S;
% 
% As = 2.5;
% 
% ns = -0.5;
% mu = 10.0;
% 
% dt = 0.01;
% 
% Es = 1;
% tMax = 5.0;
% times = 0:dt:tMax;
% Ls = ones(1, length(times));
% LHats = ones(1, length(times));
% LMins = ones(1, length(times));
% eta = 1;
% 
% netStresses = zeros(1, length(times));
% 
% for i = 2:length(times)
%     
%     
%     nOld = Es.*(1 - gammaOld)./gammaOld;
%     
%     gammaNew = gammaOld.*(1 + dt*(W(S, sigma) + mu.*(nOld - ns).*(nOld < ns)));
%     
%     nNew = Es.*(1 - gammaNew)./(gammaNew);
%     
%     gammaOld = gammaNew;
%     
%     netStresses(i) = trapz(S, nNew);
% end
% 
% plot(times, netStresses)
% 
% %%
% 
% for i = 2:length(times)
%     
%     nOld = Es.*(1 - gammaOld)./gammaOld;
%     
%     gammaNew = gammaOld.*(1 + dt*(W(S, sigma) + mu.*(nOld - ns).*(nOld < ns)));
%     
%     ANew = 2*AOld + dt - (gammaNew./gammaOld).*AOld;
%     
%     %     figure(1)
%     %     hold on
%     %     plot(S, (gammaNew - gammaOld)./(dt*gammaOld))
%     
%     gammaOld = gammaNew;
%     
%     
%     AOld = ANew;
%     
%     AgeFunct = @(x) interp1(S, AOld, x) - As;
%     
%     if ( (AgeFunct(0) < 0)&&(AgeFunct(1) > 0) )
%         
%         LMins(i) = fsolve(AgeFunct, 0.5*LMins(i - 1));
%         
%         LHats(i) = LHats(i - 1) + dt*eta*(LMins(i) - LHats(i - 1));
%         Ls(i) = Ls(i - 1) + dt*eta*(LHats(i - 1) - Ls(i - 1));
%         
%     else
%         
%         LMins(i) = LMins(i - 1);
%         LHats(i) = LHats(i - 1);
%         Ls(i) = Ls(i - 1);
%         
%     end
%     
%     
%     SOld = S;
%     SNew = SOld./SOld(end).*Ls(i);
%     gammaOld = interp1(SOld, gammaOld, SNew);
%     AOld = interp1(SOld, AOld, SNew);
%     S = SNew;
%     
% end
% 
% figure(2)
% plot(times, Ls)
% 
% 
% %%

%%

LOld = 1;
gammaOld = 1;
sigma = 0.16;

W = @(S, width) exp(-((S)./width).^2);

S = 0:1e-3:LOld;
currentS = S;

dt = 0.01;

Es = 1;
tMax = 5.0;
t = 0:dt:tMax;
Ls = ones(1, length(t));
LHats = ones(1, length(t));

T = 2.5;
for i = 2:length(t)
    
    gammaNew = gammaOld.*(1 + dt*W(currentS, sigma));
    
    % Update sloughing if we are past the critical time
    if (t(i) > T)
        LNew = LOld - (1/LOld).*trapz(S, gammaNew - gammaOld);
        
        SNew = S./S(end).*LNew;
        gammaNew = interp1(S, gammaNew, SNew);
        
        S = SNew;
        
    else
        LNew = LOld;
    end
    
    currentS = cumtrapz(S, gammaNew);
    
    Ls(i) = currentS(end);
    LHats(i) = LNew;
    
    gammaOld = gammaNew;
    
    LOld = LNew;
end