clc;clearvars; close all; warning('off','all');
ch_func = Channel_functions();
est_func = Estimation_functions();
%% --------OFDM Parameters - Given in IEEE 802.11p Spec--
ofdmBW                 = 10 * 10^6 ;                                % OFDM bandwidth (Hz)
nFFT                   = 64;                                        % FFT size 
nDSC                   = 48;                                        % Number of data subcarriers
nPSC                   = 4;                                         % Number of pilot subcarriers
nZSC                   = 12;                                        % Number of zeros subcarriers
nUSC                   = nDSC + nPSC;                               % Number of total used subcarriers
K                      = nUSC + nZSC;                               % Number of total subcarriers
nSym                   = 100;                                       % Number of OFDM symbols within one frame
deltaF                 = ofdmBW/nFFT;                               % Bandwidth for each subcarrier - include all used and unused subcarriers 
Tfft                   = 1/deltaF;                                  % IFFT or FFT period = 6.4us
Tgi                    = Tfft/4;                                    % Guard interval duration - duration of cyclic prefix - 1/4th portion of OFDM symbols = 1.6us
Tsignal                = Tgi+Tfft;                                  % Total duration of BPSK-OFDM symbol = Guard time + FFT period = 8us
K_cp                   = nFFT*Tgi/Tfft;                             % Number of symbols allocated to cyclic prefix 
pilots_locations       = [8,22,44,58].';                            % Pilot subcarriers positions
pilots                 = [1 1 1 -1].';
data_locations         = [2:7, 9:21, 23:27, 39:43, 45:57, 59:64].'; % Data subcarriers positions
null_locations         = [1, 28:38].';
ppositions             = [7,21, 32,46].';                           % Pilots positions in Kset
dpositions             = [1:6, 8:20, 22:31, 33:45, 47:52].';        % Data positions in Kset
% Pre-defined preamble in frequency domain
dp = [ 0  0 0 0 0 0 +1 +1 -1 -1 +1  +1 -1  +1 -1 +1 +1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1 +1 +1 +1 0 +1 -1 -1 +1 +1 -1 +1 -1 +1 -1 -1 -1 -1 -1 +1 +1 -1 -1 +1 -1 +1 -1 +1 +1 +1 +1 0 0 0 0 0];
Ep                     = 1;                                         % pramble power per sample
dp                     = fftshift(dp);                              % Shift zero-frequency component to center of spectrum    
predefined_preamble    = dp;
Kset                   = find(dp~=0);     
predefined_preamble_Kon = predefined_preamble(Kset).';% set of allocated subcarriers                  
Kon                    = length(Kset);                              % Number of active subcarriers
dp                     = sqrt(Ep)*dp.';
xp                     = sqrt(K)*ifft(dp);
xp_cp                  = [xp(end-K_cp+1:end); xp];                  % Adding CP to the time domain preamble
preamble_80211p        = repmat(xp_cp,1,2);                         % IEEE 802.11p preamble symbols (tow symbols)
DiagP                  = diag(predefined_preamble_Kon);
DiagP_Inv              = inv(DiagP);
%% ------ Bits Modulation Technique------------------------------------------
modu                      = 'QPSK';
Mod_Type                  = 1;              % 0 for BPSK and 1 for QAM 
if(Mod_Type == 0)
    nBitPerSym            = 1;
    Pow                   = 1;
    %BPSK Modulation Objects
    bpskModulator         = comm.BPSKModulator;
    bpskDemodulator       = comm.BPSKDemodulator;
    M                     = 1;
elseif(Mod_Type == 1)
    if(strcmp(modu,'QPSK') == 1)
         nBitPerSym       = 2; 
    elseif (strcmp(modu,'16QAM') == 1)
         nBitPerSym       = 4; 
    elseif (strcmp(modu,'64QAM') == 1)
         nBitPerSym       = 6; 
    end
    M                     = 2 ^ nBitPerSym; % QAM Modulation Order   
    Pow                   = mean(abs(qammod(0:(M-1),M)).^2); % Normalization factor for QAM    
    Constellation         =  1/sqrt(Pow) * qammod(0:(M-1),M); % 
end
 Constellation_Matrix = repmat(Constellation,nUSC,1);
%% ---------Scrambler Parameters---------------------------------------------
scramInit                 = 93; % As specidied in IEEE 802.11p Standard [1011101] in binary representation
%% ---------Convolutional Coder Parameters-----------------------------------
constlen                  = 7;
trellis                   = poly2trellis(constlen,[171 133]);
tbl                       = 34;
rate                      = 1/2;
%% -------Interleaver Parameters---------------------------------------------
% Matrix Interleaver
Interleaver_Rows          = 16;
Interleaver_Columns       = (nBitPerSym * nDSC * nSym) / Interleaver_Rows;
% General Block Interleaver
Random_permutation_Vector = randperm(nBitPerSym*nDSC*nSym); % Permutation vector
%% -----------------Vehicular Channel Model Parameters--------------------------
ChType                    = 'VTV_UC';             % Channel model
mobility                  = 'Low';
ch_l                      = [6, 7, 8, 9];
fs                        = K*deltaF;               % Sampling frequency in Hz, here case of 802.11p with 64 subcarriers and 156250 Hz subcarrier spacing
fc                        = 5.9e9;                  % Carrier Frequecy in Hz.
vel                       = 45;                    % Moving speed of user in km
c                         = 3e8;                    % Speed of Light in m/s
fD                        = 250;%(vel/3.6)/c*fc;         % Doppler freq in Hz
rchan                     = ch_func.GenFadingChannel(ChType, fD, fs);
release(rchan);
rchan.Seed                = 22;
%% ---------Bit to Noise Ratio------------------%
EbN0dB                    = 0:5:40;         % bit to noise ratio
SNR_p                     = EbN0dB + 10*log10(K/nDSC) + 10*log10(K/(K + K_cp)) + 10*log10(nBitPerSym) + 10*log10(rate);
SNR_p                     = SNR_p.';
N0                        = Ep*10.^(-SNR_p/10);


ChType2                    = 'VTV_UC';  
fD2                        = 250;
rchan2                     = ch_func.GenFadingChannel(ChType2, fD2, fs);
release(rchan2);
rchan2.Seed                = 22;
Rhphp = est_func.Estimat_Rhphp (rchan2, K_cp, K, 1, pilots_locations);
release(rchan2);
rchan2.Seed                = 22;
Rsubhp = est_func.Estimat_Rsubhp (rchan2, K_cp, K, 1, pilots_locations, data_locations);
release(rchan2);
rchan2.Seed                = 22;
%% Simulation Parameters 
N_CH                         = 10; 
N_SNR                        = length(SNR_p);
Ber_Ideal                    = zeros(N_SNR,1);
Ber_LSP                      = zeros(N_SNR,1);
Ber_STA                      = zeros(N_SNR,1);
Ber_CDP                      = zeros(N_SNR,1);
Ber_TRFI                     = zeros(N_SNR,1);
Ber_iCDP                     = zeros(N_SNR,1);
Ber_T_DFT                    = zeros(N_SNR,1);
Ber_TA_TDFT                  = zeros(N_SNR,1);
Ber_LMMSE                    = zeros(N_SNR,1);

Err_LSP                      = zeros(N_SNR,1);
Err_STA                      = zeros(N_SNR,1);
Err_CDP                      = zeros(N_SNR,1);
Err_TRFI                     = zeros(N_SNR,1);
Err_iCDP                     = zeros(N_SNR,1);
Err_T_DFT                    = zeros(N_SNR,1);
Err_TA_TDFT                  = zeros(N_SNR,1);
Err_LMMSE                    = zeros(N_SNR,1);
%% 1D DFT Interpolation
D                         = dftmtx (nFFT);
Dt                        = D(:,ch_l+1);
Dp                        = D(pilots_locations,ch_l+1);
temp                      = ((Dp' * Dp)^-1) * Dp';
H_Interpolation           = Dt * temp;
Phf_H_Total                  = zeros(N_SNR,1);

alpha = 2;
Beta = 2;
w = 1 / (2*Beta +1);
lambda = -Beta:Beta;
%% Simulation Loop
for n_snr = 1:N_SNR
    disp(['Running Simulation, SNR = ', num2str(EbN0dB(n_snr))]);
     tic;      
     release(rchan);
     rchan.Seed                = 22;
     
      R_Estimated = ((Rhphp) + N0(n_snr) * eye(nPSC))^(-1);
     
    for n_ch = 1:N_CH % loop over channel realizations
        % Bits Stream Generation 
        Bits_Stream_Coded = randi(2, nDSC * nSym  * nBitPerSym * rate,1)-1;
        % Data Scrambler 
        scrambledData = wlanScramble(Bits_Stream_Coded,scramInit);
        % Convolutional Encoder
        dataEnc = convenc(scrambledData,trellis);
        % Interleaving
        % Matrix Interleaving
        codedata = dataEnc.';
        Matrix_Interleaved_Data = matintrlv(codedata,Interleaver_Rows,Interleaver_Columns).';
        % General Block Interleaving
        General_Block_Interleaved_Data = intrlv(Matrix_Interleaved_Data,Random_permutation_Vector);
        % Bits Mapping: M-QAM Modulation
        TxBits_Coded = reshape(General_Block_Interleaved_Data,nDSC , nSym  , nBitPerSym);
        % Gray coding goes here
        TxData_Coded = zeros(nDSC ,nSym);
        for m = 1 : nBitPerSym
           TxData_Coded = TxData_Coded + TxBits_Coded(:,:,m)*2^(m-1);
        end
        % M-QAM Modulation
         Modulated_Bits_Coded  =1/sqrt(Pow) * qammod(TxData_Coded,M);
         % Check the power of Modulated bits. it must be equal to 1
         avgPower = mean (abs (Modulated_Bits_Coded).^ 2,2);
         % OFDM Frame Generation
         OFDM_Frame_Coded = zeros(K,nSym);
         OFDM_Frame_Coded(data_locations,:) = Modulated_Bits_Coded;
         OFDM_Frame_Coded(pilots_locations,:) = repmat(pilots,1,nSym);
         % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
         IFFT_Data_Coded = sqrt(K)*ifft(OFDM_Frame_Coded);
         % checking the power of the transmit signal (it has to be 1 after normalization)
         power_Coded = var(IFFT_Data_Coded(:)) + abs(mean(IFFT_Data_Coded(:)))^2; 
         % Appending cylic prefix
         CP_Coded = IFFT_Data_Coded((K - K_cp +1):K,:);
         IFFT_Data_CP_Coded = [CP_Coded; IFFT_Data_Coded];
         % Appending preamble symbol 
         %IFFT_Data_CP_Preamble_Coded = [xp_cp IFFT_Data_CP_Coded];
         IFFT_Data_CP_Preamble_Coded = [ preamble_80211p IFFT_Data_CP_Coded];
         power_transmitted =  var(IFFT_Data_CP_Preamble_Coded(:)) + abs(mean(IFFT_Data_CP_Preamble_Coded(:)))^2; 

        % ideal estimation
        [ h, y ] = ch_func.ApplyChannel( rchan, IFFT_Data_CP_Preamble_Coded, K_cp);
        release(rchan);
        rchan.Seed = rchan.Seed+1;
        yp = y((K_cp+1):end,1:2);
        y  = y((K_cp+1):end,3:end);
        
        yFD = sqrt(1/K)*fft(y);
        yfp = sqrt(1/K)*fft(yp); % FD preamble
        
        
        h = h((K_cp+1):end,:);
        hf = fft(h); % Fd channel
        hfp1 = hf(:,1);
        hfp2 = hf(:,2);
        hfp = (hfp1 + hfp2) ./2;
        hf  = hf(:,3:end);        
      
        Phf_H_Total(n_snr,1) = Phf_H_Total(n_snr,1) + mean(sum(abs(hf(Kset,:)).^2));
       
        %add noise
        noise_preamble = sqrt(N0(n_snr))*ch_func.GenRandomNoise([K,2], 1);
        yfp_r = yfp +  noise_preamble;
        noise_OFDM_Symbols = sqrt(N0(n_snr))*ch_func.GenRandomNoise([K,size(yFD,2)], 1);
        y_r   = yFD + noise_OFDM_Symbols;
        HLS_P = y_r(pilots_locations,:) ./pilots;
       %% Channel Estimation
       %% IEEE 802.11p LS Estimation at Preambles
        he_LS_Preamble = ((yfp_r(Kset,1) + yfp_r(Kset,2))./(2.*predefined_preamble(Kset).'));   
        hfe_LSP_Frame = repmat(he_LS_Preamble,1,nSym); 
        Err_LSP (n_snr,1) = Err_LSP (n_snr,1) + mean(sum(abs( hfe_LSP_Frame - hf(Kset,:)).^2));
         %% STA Channel Estimation
         [H_STA, Equalized_OFDM_Symbols_STA] = STA(he_LS_Preamble ,y_r, Kset,modu, nUSC, nSym, ppositions, alpha, w, lambda);        
         Err_STA(n_snr) = Err_STA(n_snr) + mean(sum(abs(H_STA - hf(Kset,:)).^2));        
         %% CDP Channel Estimation
         [H_CDP, Equalized_OFDM_Symbols_CDP] = CDP(he_LS_Preamble([1:6, 8:20, 22:31, 33:45, 47:52].',1) ,y_r, data_locations, yfp_r(data_locations,2), predefined_preamble(1,data_locations).', modu, nDSC, nSym);
         Err_CDP(n_snr) = Err_CDP(n_snr) + mean(sum(abs(H_CDP - hf(data_locations,:)).^2));                 
         %% iCDP Channel Estimation
         [H_iCDP, Equalized_OFDM_Symbols_iCDP] = iCDP(he_LS_Preamble ,y_r, Kset, yfp_r(Kset,2), predefined_preamble(1,Kset).', modu, nUSC, nSym, ppositions, alpha, w, lambda, Constellation_Matrix, Constellation);
         Err_iCDP(n_snr) = Err_iCDP(n_snr) + mean(sum(abs(H_iCDP - hf(Kset,:)).^2));    
         %% TRFI Channel Estimation 
          [H_TRFI, Equalized_OFDM_Symbols_TRFI] = TRFI(he_LS_Preamble ,y_r, Kset, yfp_r(Kset,2), predefined_preamble(1,Kset).', ppositions, modu, nUSC, nSym);
          Err_TRFI(n_snr) = Err_TRFI(n_snr) + mean(sum(abs(H_TRFI - hf(Kset,:)).^2));        
         %% Proposed T-DFT & TA-TDFT
         H_T_DFT = H_Interpolation * HLS_P; 
         Err_T_DFT (n_snr,1) = Err_T_DFT (n_snr,1) + mean(sum(abs(H_T_DFT(Kset,:) - hf(Kset,:)).^2));         
         H_TA_TDFT = TA_T_DFT(H_T_DFT, Kset, nUSC, nSym, alpha, w, lambda);
         Err_TA_TDFT(n_snr,1) = Err_TA_TDFT(n_snr,1) + mean(sum(abs(H_TA_TDFT - hf(Kset,:)).^2)) ; 
        %% 1D-LMMSE
        H_LMMSE = (Rsubhp * R_Estimated) * HLS_P;
        Err_LMMSE(n_snr,1) = Err_LMMSE(n_snr,1) + mean(sum(abs(H_LMMSE - hf(data_locations,:)).^2)) ;
                   
        %%    IEEE 802.11p Rx     
        Bits_Ideal                                          = de2bi(qamdemod(sqrt(Pow) * (y_r(data_locations ,:) ./ hf(data_locations,:)),M)); 
        Bits_LSP                                            = de2bi(qamdemod(sqrt(Pow) * (y_r(data_locations ,:) ./ hfe_LSP_Frame(dpositions,:)),M)); 
        Bits_STA                                            = de2bi(qamdemod(sqrt(Pow) * Equalized_OFDM_Symbols_STA(dpositions,:),M)); 
        Bits_CDP                                            = de2bi(qamdemod(sqrt(Pow) * Equalized_OFDM_Symbols_CDP,M)); 
        Bits_iCDP                                           = de2bi(qamdemod(sqrt(Pow) * Equalized_OFDM_Symbols_iCDP(dpositions,:),M));  
        Bits_T_DFT                                          = de2bi(qamdemod(sqrt(Pow) * (y_r(data_locations ,:) ./ H_T_DFT(data_locations,:)),M));  
        Bits_TA_TDFT                                        = de2bi(qamdemod(sqrt(Pow) * (y_r(data_locations ,:) ./ H_TA_TDFT(dpositions,:)),M));  
        Bits_TRFI                                           = de2bi(qamdemod(sqrt(Pow) * Equalized_OFDM_Symbols_TRFI(dpositions,:),M));                            
        Bits_LMMSE                                          = de2bi(qamdemod(sqrt(Pow) * (y_r(data_locations ,:) ./ H_LMMSE),M));
                
        Ber_Ideal (n_snr)                                   = Ber_Ideal (n_snr) + biterr(wlanScramble((vitdec((matintrlv((deintrlv(Bits_Ideal(:),Random_permutation_Vector)).',Interleaver_Columns,Interleaver_Rows).'),trellis,tbl,'trunc','hard')),scramInit),Bits_Stream_Coded);
        Ber_LSP (n_snr)                                     = Ber_LSP (n_snr) + biterr(wlanScramble((vitdec((matintrlv((deintrlv(Bits_LSP(:),Random_permutation_Vector)).',Interleaver_Columns,Interleaver_Rows).'),trellis,tbl,'trunc','hard')),scramInit),Bits_Stream_Coded);
        Ber_STA(n_snr)                                      = Ber_STA(n_snr) + biterr( wlanScramble((vitdec((matintrlv((deintrlv(Bits_STA(:),Random_permutation_Vector)).',Interleaver_Columns,Interleaver_Rows).'),trellis,tbl,'trunc','hard')),scramInit),Bits_Stream_Coded);
        Ber_CDP(n_snr)                                      = Ber_CDP(n_snr) + biterr(wlanScramble((vitdec((matintrlv((deintrlv(Bits_CDP(:),Random_permutation_Vector)).',Interleaver_Columns,Interleaver_Rows).'),trellis,tbl,'trunc','hard')),scramInit),Bits_Stream_Coded);
        Ber_iCDP(n_snr)                                     = Ber_iCDP(n_snr) + biterr(wlanScramble((vitdec((matintrlv((deintrlv(Bits_iCDP(:),Random_permutation_Vector)).',Interleaver_Columns,Interleaver_Rows).'),trellis,tbl,'trunc','hard')),scramInit),Bits_Stream_Coded);
        Ber_T_DFT(n_snr)                                    = Ber_T_DFT(n_snr) + biterr(wlanScramble((vitdec((matintrlv((deintrlv(Bits_T_DFT(:),Random_permutation_Vector)).',Interleaver_Columns,Interleaver_Rows).'),trellis,tbl,'trunc','hard')),scramInit),Bits_Stream_Coded);
        Ber_TA_TDFT(n_snr)                                  = Ber_TA_TDFT(n_snr) + biterr(wlanScramble((vitdec((matintrlv((deintrlv(Bits_TA_TDFT(:),Random_permutation_Vector)).',Interleaver_Columns,Interleaver_Rows).'),trellis,tbl,'trunc','hard')),scramInit),Bits_Stream_Coded);  
        Ber_LMMSE (n_snr)                                   = Ber_LMMSE (n_snr) + biterr(wlanScramble((vitdec(Bits_LMMSE(:),trellis,tbl,'trunc','hard')),scramInit),Bits_Stream_Coded);        
        Ber_TRFI(n_snr)                                     = Ber_TRFI(n_snr) + biterr(wlanScramble((vitdec((matintrlv((deintrlv(Bits_TRFI(:),Random_permutation_Vector)).',Interleaver_Columns,Interleaver_Rows).'),trellis,tbl,'trunc','hard')),scramInit),Bits_Stream_Coded);
        
    end            
toc;    
end 

 
BER_Ideal                     = Ber_Ideal /(N_CH * nSym * nDSC * nBitPerSym);
BER_LSP                       = Ber_LSP / (N_CH * nSym * nDSC * nBitPerSym);
BER_STA                       = Ber_STA / (N_CH * nSym * nDSC * nBitPerSym);
BER_CDP                       = Ber_CDP / (N_CH * nSym * nDSC * nBitPerSym);
BER_iCDP                      = Ber_iCDP / (N_CH * nSym * nDSC * nBitPerSym);
BER_TRFI                      = Ber_TRFI / (N_CH * nSym * nDSC * nBitPerSym);
BER_T_DFT                     = Ber_T_DFT / (N_CH * nSym * nDSC * nBitPerSym);
BER_TA_TDFT                   = Ber_TA_TDFT / (N_CH * nSym * nDSC * nBitPerSym);
BER_LMMSE                     = Ber_LMMSE  / (N_CH * nSym * nDSC * nBitPerSym);

% NMSE Plot 
Phf_H                                 = Phf_H_Total/(N_CH);
ERR_LSP                               = Err_LSP ./ (Phf_H * N_CH);
ERR_STA                               = Err_STA ./ (Phf_H * N_CH);
ERR_CDP                               = Err_CDP ./ (Phf_H * N_CH);
ERR_iCDP                              = Err_iCDP ./ (Phf_H * N_CH);
ERR_TRFI                              = Err_TRFI ./ (Phf_H * N_CH);
ERR_T_DFT                             = Err_T_DFT ./ (Phf_H * N_CH);
ERR_TA_TDFT                           = Err_TA_TDFT ./ (Phf_H * N_CH);
ERR_LMMSE                             = Err_LMMSE ./ (Phf_H * N_CH);

%% 

figure(1)
semilogy(EbN0dB, BER_Ideal,'k-o','LineWidth',2);
hold on;
semilogy(EbN0dB, BER_LSP,'k--','LineWidth',2);
hold on;
semilogy(EbN0dB, BER_STA,'g--*','LineWidth',2);
hold on;
semilogy(EbN0dB, BER_CDP,'b--*','LineWidth',2);
hold on;
semilogy(EbN0dB, BER_TRFI,'m--*','LineWidth',2);
hold on;
semilogy(EbN0dB, BER_iCDP,'k--*','LineWidth',2);
hold on;
semilogy(EbN0dB, BER_T_DFT,'r--o','LineWidth',2);
hold on;
semilogy(EbN0dB, BER_TA_TDFT,'c--o','LineWidth',2);
hold on;
semilogy(EbN0dB, BER_LMMSE,'k--^','LineWidth',2);
hold on;
grid on;
legend('Perfect Channel','LS-Preamble','STA','CDP','TRFI','iCDP','Proposed T-DFT','Proposed TA-TDFT','1D-LMMSE');
xlabel('SNR(dB)');
ylabel('Bit Error Rate (BER)');
title(['fD = ',num2str(fD), ' Hz, Modulation = QPSK']);

figure(2),
semilogy(EbN0dB, ERR_LSP,'k--','LineWidth',2);
hold on;
semilogy(EbN0dB, ERR_STA,'g--*','LineWidth',2);
hold on;
semilogy(EbN0dB, ERR_CDP,'b--*','LineWidth',2);
hold on;
semilogy(EbN0dB, ERR_TRFI,'m--*','LineWidth',2);
hold on;
semilogy(EbN0dB, ERR_iCDP,'k--*','LineWidth',2);
hold on;
semilogy(EbN0dB, ERR_T_DFT,'r--*','LineWidth',2);
hold on;
semilogy(EbN0dB, ERR_TA_TDFT,'c--o','LineWidth',2);
hold on;
semilogy(EbN0dB, ERR_LMMSE,'k--^','LineWidth',2);
hold on;
grid on;
legend('LS-Preamble','STA','CDP','TRFI','iCDP','Proposed T-DFT','Proposed TA-TDFT','1D-LMMSE');
xlabel('SNR(dB)');
ylabel('Normalized Mean Sqaure Error (NMSE)');
title(['fD = ',num2str(fD), ' Hz, Modulation = QPSK']);


