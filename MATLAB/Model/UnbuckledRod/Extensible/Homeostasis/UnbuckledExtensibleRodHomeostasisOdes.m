function Odes = UnbuckledExtensibleRodHomeostasisOdes(s, M, parameters)
% Equations governing homeostasis for a growing 1D rod, where
% mechanochemical growth is assumed.

K = parameters.K;
nu = parameters.nu;
W = parameters.W;
mu = parameters.mu;
ns = parameters.ns;
Es = parameters.Es;

    function dMdS = HomeostasisEqns(s, M)
        
        % Define the state variables
        N = M(1,:);
        U = M(2,:);
        V = M(3,:);
        
        if (s == 0)
            dUdS = K.*(W(0) + mu*(N - ns).*0.5*(1 + tanh(50.*(ns - N))))./(nu + W(0) + mu*(N - ns).*0.5*(1 + tanh(50.*(ns - N))));
        else
            dUdS = K + nu*U./V;
        end
        
        dNdS = U;
        dVdS = U.*V./(Es + N) - (W(s) + mu*(N - ns).*0.5*(1 + tanh(50.*(ns - N))));
        
        dMdS = [dNdS; dUdS; dVdS];
        
    end

Odes = HomeostasisEqns(s, M);

end
