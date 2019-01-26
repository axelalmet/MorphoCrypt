function Odes = KellerModelContactOdes(x, M, modelParams)

P = modelParams.P;
L = modelParams.L;

    function dMdS = ContactEqns(x, M)
        
        % Define the state variables
        S = M(1,:);
        x = M(2,:);
        y = M(3,:);
        Q = M(4,:);
        N = M(5,:);
        theta = M(6,:);
        k = M(7,:);
        
        dSdS = (L).*ones(1, length(S));
        dxdS = (L).*cos(theta);
        dydS = (L).*sin(theta);
        dQdS = -(L).*k.*N;
        dNdS = (L).*(k.*Q + P);
        dthetadS = (L).*k;
        dkdS = (L).*N;
                
        dMdS = [dSdS; dxdS; dydS; dQdS; dNdS; dthetadS; dkdS];
    end

Odes = ContactEqns(x, M);

end
