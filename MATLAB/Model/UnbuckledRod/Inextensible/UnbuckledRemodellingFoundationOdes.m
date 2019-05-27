function Odes = UnbuckledRemodellingFoundationOdes(x, M, MOld, parameters)

K = parameters.K;
L = parameters.L;
gamma = parameters.gamma;

Px = parameters.Px;

    function dMdS = FoundationEqns(x, M)
        
        % Define the state variables
        S = M(1,:);
        X = M(2,:);
        Y = M(3,:);
        F = M(4,:);
        
        if (length(gamma) > 1)
            gamma = interp1(MOld.x, gamma, x);
        end
        
        if (length(K) > 1)
            K = interp1(MOld.x, K, x);
        end
        
        if (length(Px) > 1)
            Px = interp1(MOld.x, Px, x);
        end
        
        dSdS = L.*ones(1, length(S));
        dxdS = L.*gamma.*ones(1, length(S));
        dydS = zeros(1, length(X));
        dFdS = L.*K.*gamma.*(X - Px);

        dMdS = [dSdS; dxdS; dydS; dFdS];
        
    end

Odes = FoundationEqns(x, M);

end
