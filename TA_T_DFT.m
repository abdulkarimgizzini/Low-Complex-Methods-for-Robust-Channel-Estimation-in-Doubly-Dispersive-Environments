function [H_STA_DFT] = TA_T_DFT(H_DFT_Interpolation_ZF, Kset, nUSC, nSym, alpha, w, lambda)
H_STA_DFT = zeros(nUSC, nSym);
for i = 1:nSym
    if(i == 1)
        Initial_Channel_Estimate  = H_DFT_Interpolation_ZF(Kset,i);
        Abed_TA = Initial_Channel_Estimate;
        H_STA_DFT(:,i) = Abed_TA;        
    elseif (i > 1)
        Initial_Channel_Estimate = H_DFT_Interpolation_ZF(Kset,i);       
        % Time domain Averging
        Abed_TA = (1 - (1/alpha)) .*  H_STA_DFT(:,i-1) + (1/alpha).* Initial_Channel_Estimate;
        % Update H_Abed Matrix
        H_STA_DFT(:,i) = Abed_TA; 
    end
end

end

