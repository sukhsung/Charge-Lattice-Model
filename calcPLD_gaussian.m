function posLT_PLD = calcPLD_gaussian(posLT, posSL, r_c, sg, pldScal)
    % Calculate Associate Lattice Distortion from Charge Lattice
    % Assume Gaussian Cutoff
    
    numSiteLT = size(posLT,1);
    posLT_PLD = posLT;
    
    for indLT = 1:numSiteLT

        pos_cur = posLT(indLT,:);
        
        dR_nns = pos_cur - posSL;

        r_nns = sqrt(dR_nns(:,1).^2 + dR_nns(:,2).^2);
        
        dR_nns( r_nns > r_c, : ) = [];
        r_nns = sqrt(dR_nns(:,1).^2 + dR_nns(:,2).^2);
        
        
        pldDis = pldScal*sum(dR_nns.*exp(-0.5*r_nns.^2/(sg^2)),1);

        posLT_PLD(indLT,:) = pos_cur - pldDis;

        
    end
        

end