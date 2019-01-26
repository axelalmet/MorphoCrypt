function Odes = KellerModelHalfIntervalContactRegionOdes(x, M, region, modelParams)

P = modelParams.P;
L = modelParams.L;

    function dMdS = ContactEqns(x, M, region)
        
        % Define the state variables
        S = M(1,:);
        x = M(2,:);
        y = M(3,:);
        Q = M(4,:);
        N = M(5,:);
        theta = M(6,:);
        k = M(7,:);
        sc = M(8,:);
        
        dscdS = zeros(1, length(S));
        
        switch region
            case 1 % ODEs for [0, sc]
                
                dSdS = (0.5*L - sc).*ones(1, length(S));
                dxdS = (0.5*L - sc).*cos(theta);
                dydS = (0.5*L - sc).*sin(theta);
                dQdS = -(0.5*L - sc).*k.*N;
                dNdS = (0.5*L - sc).*(k.*Q + P);
                dthetadS = (0.5*L - sc).*k;
                dkdS = (0.5*L - sc).*N;
                
            case 2 % ODEs for [sc, L]
                
                dSdS = (0.5*L - sc).*ones(1, length(S));
                dxdS = (0.5*L - sc).*cos(theta);
                dydS = (0.5*L - sc).*sin(theta);
                dQdS = -(0.5*L - sc).*k.*N;
                dNdS = (0.5*L - sc).*(k.*Q + P);
                dthetadS = (0.5*L - sc).*k;
                dkdS = (0.5*L - sc).*N;
                
            otherwise
                
                error('MATLAB:contactOdes:BadRegionIndex','Incorrect region index: %d',region);
                
        end
                
        dMdS = [dSdS; dxdS; dydS; dQdS; dNdS; dthetadS; dkdS; dscdS];

    end

Odes = ContactEqns(x, M, region);

end
