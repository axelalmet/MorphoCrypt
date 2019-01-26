function Odes = KellerModelDuringRegionContactOdes(x, M, sc, modelParams)

P = modelParams.P;

    function dMdS = ContactEqns(x, M)
        
        % Define the state variables
        S = M(1,:);
        x = M(2,:);
        y = M(3,:);
        Q = M(4,:);
        N = M(5,:);
        theta = M(6,:);
        k = M(7,:);
                
        dSdS = (2*sc).*ones(1, length(S));
        dxdS = zeros(1, length(S));
        dydS = (2*sc).*ones(1, length(S));
        dQdS = zeros(1, length(S));
        dNdS = (2*sc).*(P);
        dthetadS = zeros(1, length(S));
        dkdS = zeros(1, length(S));
        
                
        dMdS = [dSdS; dxdS; dydS; dQdS; dNdS; dthetadS; dkdS];

    end

Odes = ContactEqns(x, M);

end
