function Odes = BuckledExtensibleRodHomeostasisOdes(s, M, solOld, parameters)
% Equations governing homeostasis for a growing 2D rod, where
% mechanochemical growth is assumed.

% Get the relevant parameters
K = parameters.K;
nu = parameters.nu;
W = parameters.W;
mu = parameters.mu;
ns = parameters.ns;
Es = parameters.Es;
sigma = parameters.sigma;

% Get the homeostatic shape
sOld = parameters.currentArcLength;
xOld = solOld.y(2,:);
yOld = solOld.y(3,:);
thetaOld = solOld.y(6,:);
thetaPOld = solOld.y(7,:);

l = sOld(end);

    function dMdS = HomeostasisEqns(s, M)
        
        % Define the state variables
        Nx = M(1,:); % Horizontal stress component
        Ny = M(2,:); % Vertical stress component
        Px = M(3,:); % Horizontal foundation attachment component
        Py = M(4,:); % Vertical foundation attachment component
        V = M(5,:); % Flow velocity
        
        X = interp1(sOld, xOld, s);
        Y = interp1(sOld, yOld, s);
        Theta = interp1(sOld, thetaOld, s);
        ThetaP = interp1(sOld, thetaPOld, s);
        
        % Have to be careful with behaviour of Px' and Py' at s = 0. Values
        % are found by application of L'Hopital's rule.
        
        if (s == 0)
            
            dPxdS = nu./(nu + W(0, sigma) + mu*(Nx - ns).*0.5*(1 + tanh(50.*(ns - Nx))));
            dPydS = 0;
            
        else
            
            dPxdS = nu*(Px - X)./V;
            dPydS = nu*(Py - Y)./V;
            
        end
        
        dNxdS = K.*(X - Px);
        dNydS = K.*(Y - Py);
        
        dVdS = (K.*(X - Px).*cos(Theta) + K.*(Y - Py).*sin(Theta) + ...
            ThetaP./(1 + Es*(Nx.*cos(Theta) + Ny.*sin(Theta))).*(Ny.*cos(Theta) - Nx.*sin(Theta))).*V ...
                ./(1 + Es.*(Nx.*cos(Theta) + Ny.*sin(Theta))) ...
                - W(s, sigma).*(1 + mu*tanh(Nx.*cos(Theta) + Ny.*sin(Theta) - ns).*(Nx.*cos(Theta) + Ny.*sin(Theta) < ns));
        
        dMdS = l.*[dNxdS; dNydS; dPxdS; dPydS; dVdS];
        
    end

Odes = HomeostasisEqns(s, M);

end
