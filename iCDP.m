function [H_iCDP, Equalized_OFDM_Symbols] = iCDP(he_LS_Preamble ,y_r, Kset, rp2, pp, mod, nUSC, nSym, ppositions, alpha, w, lambda, Constellation, Constellation_Vector)

H_iCDP = zeros(nUSC, nSym);
H_STA = zeros(nUSC, nSym);
% Constellation



Equalized_OFDM_Symbols = zeros(nUSC, nSym);
for i = 1:nSym
    if(i == 1)
        % Step 1: Equalization. First, the ith received symbol, Received_OFDM_Symbols(i,k) is equalized by the channel estimate, 
        % of the previous symbol. For i = 1, CDP(0) is the LS estimate at
        % the preamble.
        Equalized_OFDM_Symbol = y_r(Kset,i)./ he_LS_Preamble;
        Equalized_OFDM_Symbols(:,i) = Equalized_OFDM_Symbol;
       
        % Step 2: Constructing Data Pilot. Equalized_OFDM_Symbol is demodulated to obtain De_Equalized_OFDM_Symbol.
        De_Equalized_OFDM_Symbol = wlanClosestReferenceSymbol(Equalized_OFDM_Symbol,mod);
        De_Equalized_OFDM_Symbol(ppositions,:) = [1;1;1;-1]; 
        % Step 3: LS Estimation. the ith received symbol is divided by  De_Equalized_OFDM_Symbol to obtain the initial channel estimate
        Initial_Channel_Estimate = y_r(Kset,i)./ De_Equalized_OFDM_Symbol;
        
        
        
        % STA
        STA_FA = zeros(nUSC,1);
        STA_FA(1,1) = Initial_Channel_Estimate(1,1);
        STA_FA(2,1) = Initial_Channel_Estimate(2,1);
        STA_FA(51,1) = Initial_Channel_Estimate(51,1);
        STA_FA(52,1) = Initial_Channel_Estimate(52,1);
        for j = 3:nUSC -2
           for l = 1: size(lambda,2)
               STA_FA(j,1) = STA_FA(j,1) + (w * Initial_Channel_Estimate(j + lambda(1,l), 1)); 
           end 
        end
        
        STA_TA = (1 - (1/alpha)) .* he_LS_Preamble + (1/alpha).* STA_FA;
        H_STA(:,i) = STA_TA; 
        
        
        
        % Step 4: Equalization and Demapping
        E1_0 = real(rp2 ./ Initial_Channel_Estimate);
        E2_0 = pp;
        X1_0 = wlanClosestReferenceSymbol(E1_0,'BPSK');
        X2_0 =  E2_0;
        
        % Step 5: Comparison
        equal_indices = find(X1_0 == X2_0);
        unequal_indices = find(X1_0 ~= X2_0);
        
        H_iCDP(equal_indices,1) =  (STA_TA(equal_indices,1)+ Initial_Channel_Estimate(equal_indices,1)) ./2;
        H_iCDP(unequal_indices,1) = he_LS_Preamble(unequal_indices,1);
        
 
    elseif (i > 1)
        % Step 1: Equalization 
        Equalized_OFDM_Symbol = y_r(Kset,i)./ H_iCDP(:,i-1);
        Equalized_OFDM_Symbols(:,i) = Equalized_OFDM_Symbol;
        % Step 2: Constructing Data Pilot.
        De_Equalized_OFDM_Symbol = wlanClosestReferenceSymbol(Equalized_OFDM_Symbol,mod);
        De_Equalized_OFDM_Symbol(ppositions,:) = [1;1;1;-1]; 
        % Step 3: LS Estimation.
        Initial_Channel_Estimate = y_r(Kset,i)./ De_Equalized_OFDM_Symbol;
        
        
        
        
        
        STA_FA = zeros(nUSC,1);
        STA_FA(1,1) = Initial_Channel_Estimate(1,1);
        STA_FA(2,1) = Initial_Channel_Estimate(2,1);
        STA_FA(51,1) = Initial_Channel_Estimate(51,1);
        STA_FA(52,1) = Initial_Channel_Estimate(52,1);
%         alpha = 2;
%         Beta = 2;
%         w = 1 / (2*Beta +1);
%         lambda = -Beta:Beta;
        for j = 3:nUSC -2
           for l = 1: size(lambda,2)
               STA_FA(j,1) = STA_FA(j,1) + (w * Initial_Channel_Estimate(j + lambda(1,l), 1)); 
           end 
        end
        

        STA_TA = (1 - (1/alpha)) .*  H_STA(:,i-1) + (1/alpha).* STA_FA;
        H_STA(:,i) = STA_TA; 
        
        
        
        % Step 4: Equalization and Demapping
        E1 =  y_r(Kset,i - 1)./ Initial_Channel_Estimate;
        E2 =  y_r(Kset,i - 1)./ H_iCDP(:,i-1);
        
        X1 =  wlanClosestReferenceSymbol(E1,mod);
        X2 =  wlanClosestReferenceSymbol(E2,mod);
        
%         [~, X1_Distances_IDX ] = min(sqrt((real(E1 - Constellation)).^2 + (imag(E1 - Constellation)).^2),[],2);
%         [~, X2_Distances_IDX] = min(sqrt((real(E2 - Constellation)).^2 + (imag(E2 - Constellation)).^2),[],2);
%         
%         
%         X1 = zeros(52,1);
%         X2 = zeros(52,1);
%         for lk = 1:52
%             X1(lk,1) = Constellation_Vector(X1_Distances_IDX(lk,1));
%             X2(lk,1) = Constellation_Vector(X2_Distances_IDX(lk,1));
%         end
%         %X1 = Constellation(:,X1_Distances_IDX);
%         %X2 = Constellation(:,X2_Distances_IDX);
        
        % Step 5: Comparison
        equal_indices = find(X1 == X2);
        unequal_indices = find(X1 ~= X2);
        
        H_iCDP(equal_indices,i) = (STA_TA(equal_indices,1) + Initial_Channel_Estimate(equal_indices,1)) ./2;
        H_iCDP(unequal_indices,i) = H_iCDP(unequal_indices,i - 1);
    end
end
end

