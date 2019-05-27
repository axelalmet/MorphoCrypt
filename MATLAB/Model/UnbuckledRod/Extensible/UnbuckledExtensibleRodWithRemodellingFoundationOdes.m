function Odes = UnbuckledExtensibleRodWithRemodellingFoundationOdes(x, M, MOld, parameters)

K = parameters.K;
L = parameters.L;
gamma = parameters.gamma;
Es = parameters.Es;

P = parameters.P;

    function dMdS = FoundationEqns(x, M)
        
        % Define the state variables
        X = M(1,:);
        F = M(2,:);
        
        if (length(gamma) > 1)
            gamma = interp1(MOld.x, gamma, x);
        end
        
        if (length(K) > 1)
            K = interp1(MOld.x, K, x);
        end
        
        if (length(Es) > 1)
            Es = interp1(MOld.x, Es, x);
        end
        
        if (length(P) > 1)
            P = interp1(MOld.x, P, x);
        end
        
        alpha = 1 + F./Es;
        
        dxdS = L.*gamma.*alpha;
        dFdS = L.*K.*gamma.*alpha.*(X - P);

        dMdS = [dxdS; dFdS];
        
    end

Odes = FoundationEqns(x, M);

end
