function Odes = BuckledInextensibleRodHomeostasisOdes(s, M, sol, parameters)
% Equations governing homeostasis for a growing 2D rod, where
% mechanochemical growth is assumed.

% Get the relevant parameters
K = parameters.K;
nu = parameters.nu;
W = parameters.W;
mu = parameters.mu;
ns = parameters.ns;
sigma = parameters.sigma;

% Get the homeostatic shape
sOld = parameters.currentArcLength;
xOld = sol.y(2,:);
yOld = sol.y(3,:);
thetaOld = sol.y(6,:);

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
        
        % Have to be careful with behaviour of Px' and Py' at s = 0. Values
        % are found by application of L'Hopital's rule.
        
        if (s == sol.x(1))
            
            dPxdS = nu./(nu + W(0, sigma) + mu.*tanh(Nx - ns));
            dPydS = 0;
            
        else
            
            dPxdS = nu*(Px - X)./V;
            dPydS = nu*(Py - Y)./V;
            
        end
        
        dNxdS = K.*(X - Px);
        dNydS = K.*(Y - Py);
        
        dVdS = -W(s, sigma).*(1 + mu*tanh(Nx.*cos(Theta) + Ny.*sin(Theta) - ns).*(Nx.*cos(Theta) + Ny.*sin(Theta) < ns));
        
        dMdS = l.*[dNxdS; dNydS; dPxdS; dPydS; dVdS];
        
    end

Odes = HomeostasisEqns(s, M);

end
