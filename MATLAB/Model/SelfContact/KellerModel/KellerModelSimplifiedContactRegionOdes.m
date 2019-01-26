function Odes = KellerModelSimplifiedContactRegionOdes(x, M, modelParams)

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
        sc = M(8,:);
        
        dsCdS = zeros(1, length(S));
                
        dSdS = (sc).*ones(1, length(S));
        dxdS = (sc).*cos(theta);
        dydS = (sc).*sin(theta);
        dQdS = -(sc).*k.*N;
        dNdS = (sc).*(k.*Q + P);
        dthetadS = (sc).*k;
        dkdS = (sc).*N;
                
        dMdS = [dSdS; dxdS; dydS; dQdS; dNdS; dthetadS; dkdS; dsCdS];
    end

Odes = ContactEqns(x, M);

end
