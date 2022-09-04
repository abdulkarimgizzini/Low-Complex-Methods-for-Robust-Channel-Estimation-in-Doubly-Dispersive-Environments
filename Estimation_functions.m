function func = Estimation_functions()
func.LS = @LS;
func.MMSE = @MMSE;
func.MMSE_matrix = @MMSE_matrix;
func.Estimat_Rh = @Estimat_Rh;
func.Estimat_Rhphp = @Estimat_Rhphp;
func.Estimat_Rhphp_Frame = @Estimat_Rhphp_Frame;
func.Estimat_Rsubhp = @Estimat_Rsubhp;
func.Estimat_Rsubhp_Frame = @Estimat_Rsubhp_Frame;
func.Estimat_Rh1 = @Estimat_Rh1;
func.Estimat_Rh12 = @Estimat_Rh12;
end
% LS estimation
function [he, err] = LS(yp, xp, h)
he = yp./xp;
err = norm(he-h)^2;
end

function [he, err] = MMSE(yp, W, h)
% W the MMSE filter generated from the preambles
he = W*yp;
err = norm (he-h)^2;
% correct bias
end

function W = MMSE_matrix (xp, Rh, SNR)
K = size(Rh,1);
% to perform A*diag(xb)^{-1}
Lamp_p_inv = repmat(1./xp', K,1);
W = (Rh/(Rh + 1/SNR *eye(K))).*Lamp_p_inv;
end

function Rh = Estimat_Rh (rchan, K_cp, K, ch_l, nSym)
NR_CH = 1000;
Rh = zeros(8,8);
% dummy signal
xp_cp = rand(K_cp+K,nSym);
for n_ch = 1:NR_CH % loop over channel realizations
  
    % ideal estimation
    ch_func = Channel_functions();
    [ h, ~ ] = ch_func.ApplyChannel( rchan, xp_cp, K_cp);
    release(rchan);
    rchan.Seed = rchan.Seed+1;
    
    h = h((K_cp+1):end,:);
    h_eliminated_Tabs = h(ch_l,:);  
   
    Rh = Rh + h_eliminated_Tabs*h_eliminated_Tabs';
end
Rh = Rh/NR_CH;
end

function Rh1 = Estimat_Rh1 (rchan, K_cp, K, nSym,index)
NR_CH = 1000;
Rh1 = zeros(12,12);
% dummy signal
xp_cp = rand(K_cp+K,nSym);
for n_ch = 1:NR_CH % loop over channel realizations
  
    % ideal estimation
    ch_func = Channel_functions();
    [ h, ~ ] = ch_func.ApplyChannel( rchan, xp_cp, K_cp);
    release(rchan);
    rchan.Seed = rchan.Seed+1;
    
    h = h((K_cp+1):end,:);
    h_eliminated_Tabs = h(1:12,index);  
   
    Rh1 = Rh1 + h_eliminated_Tabs*h_eliminated_Tabs';
end
Rh1 = Rh1/NR_CH;
end

function Rh12 = Estimat_Rh12 (rchan, K_cp, K, nSym,index1,index2)
NR_CH = 1000;
Rh12 = zeros(12,12);
% dummy signal
xp_cp = rand(K_cp+K,nSym);
for n_ch = 1:NR_CH % loop over channel realizations
  
    % ideal estimation
    ch_func = Channel_functions();
    [ h, ~ ] = ch_func.ApplyChannel( rchan, xp_cp, K_cp);
    release(rchan);
    rchan.Seed = rchan.Seed+1;
    
    h = h((K_cp+1):end,:);
    h_eliminated_Tabs1 = h(1:12,index1);
    h_eliminated_Tabs2 = h(1:12,index2);
   
    Rh12 = Rh12 + h_eliminated_Tabs1*h_eliminated_Tabs2';
end
Rh12 = Rh12/NR_CH;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Estimation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rh = Estimat_Rhphp (rchan, K_cp, K, nSym, pilots_locations)
NR_CH = 1000;
PathDelays = 1e-9.*[0, 1, 100, 101, 200, 300, 400, 401, 500, 600, 700, 701];
avgPathGains = [0, 0, -11.2,-11.2,-19,-21.9, -25.3, -25.3, -24.4, -28.0, -26.1,-26.1];   
Paths_Index = 1:12;

Rh = zeros(4,4);
% dummy signal
xp_cp = rand(K_cp+K,nSym);
for n_ch = 1:NR_CH % loop over channel realizations
  
    % ideal estimation
    ch_func = Channel_functions();
%     release(rchan);
%     Path_Index_off = randi ([1 12], 1,1);
%     Path_Index_Updated = setdiff(Paths_Index,Path_Index_off);
%     PathDelays_Updated = PathDelays(1, Path_Index_Updated);
%     avgPathGains_Updated = avgPathGains(1, Path_Index_Updated);
%     rchan.PathDelays = PathDelays_Updated;
%     rchan.AveragePathGains = avgPathGains_Updated;
%     rchan.Seed = rchan.Seed+1;

    [ h, ~ ] = ch_func.ApplyChannel( rchan, xp_cp, K_cp);
      release(rchan);
      rchan.Seed = rchan.Seed+1;
%     
    h = h((K_cp+1):end,:);
    hf = fft(h);
    hp = hf(pilots_locations,:);
    
    Rh = Rh + hp*hp';
end
Rh = Rh/NR_CH;
end

function Rh = Estimat_Rsubhp (rchan, K_cp, K, nSym, pilots_locations, data_locations)
NR_CH = 1000;
% PathDelays = 1e-9.*[0, 1, 100, 101, 200, 300, 400, 401, 500, 600, 700, 701];
% avgPathGains = [0, 0, -11.2,-11.2,-19,-21.9, -25.3, -25.3, -24.4, -28.0, -26.1,-26.1];   
% Paths_Index = 1:12;

Rh = zeros(48,4);
% dummy signal
xp_cp = rand(K_cp+K,nSym);
for n_ch = 1:NR_CH % loop over channel realizations
  
    % ideal estimation
    ch_func = Channel_functions();
%     release(rchan);
%     Path_Index_off = randi ([1 12], 1,1);
%     Path_Index_Updated = setdiff(Paths_Index,Path_Index_off);
%     PathDelays_Updated = PathDelays(1, Path_Index_Updated);
%     avgPathGains_Updated = avgPathGains(1, Path_Index_Updated);
%     rchan.PathDelays = PathDelays_Updated;
%     rchan.AveragePathGains = avgPathGains_Updated;
%     rchan.Seed = rchan.Seed+1;
    [ h, ~ ] = ch_func.ApplyChannel( rchan, xp_cp, K_cp);
      release(rchan);
      rchan.Seed = rchan.Seed+1;
    h = h((K_cp+1):end,:);
    hf = fft(h); 
    hp = hf(pilots_locations,:);
    hd = hf(data_locations,:);
    for sid = 1: nSym
        for subid = 1: 48
            Rh(subid,:) = Rh(subid,:) + (hd(subid,sid) *  hp(:,sid)');
        end
    end
end
Rh = Rh/NR_CH;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rh = Estimat_Rhphp_Frame (rchan, K_cp, K, nSym, pilots_locations)
NR_CH = 1000;
Rh = zeros(nSym*4,nSym*4);
% dummy signal
xp_cp = rand(K_cp+K,nSym);
for n_ch = 1:NR_CH % loop over channel realizations
  
    % ideal estimation
    ch_func = Channel_functions();
    [ h, ~ ] = ch_func.ApplyChannel( rchan, xp_cp, K_cp);
    release(rchan);
    rchan.Seed = rchan.Seed+1;
    
    h = h((K_cp+1):end,:);
    hf = fft(h);
    hp = hf(pilots_locations,:);
    hp = hp(:);
    Rh = Rh + hp*hp';
end
Rh = Rh/NR_CH;
end

function [Rhsub,Rhp] = Estimat_Rsubhp_Frame (rchan, K_cp, K, nSym, pilots_locations, data_locations)
NR_CH = 1000;
Rhsub = zeros(48,nSym*4);
Rhp = zeros(nSym*4,nSym*4);
ch_func = Channel_functions();
% dummy signal
xp_cp = rand(K_cp+K,nSym);
for n_ch = 1:NR_CH % loop over channel realizations
    % ideal estimation
    [ h, ~ ] = ch_func.ApplyChannel(rchan, xp_cp, K_cp);
    release(rchan);
    rchan.Seed = rchan.Seed+1;
    h = h((K_cp+1):end,:);
    hf = fft(h); 
    hp = hf(pilots_locations,:);
    hp = hp(:);
    Rhp = Rhp + hp*hp';
    hd = hf(data_locations,:);
    for sid = 1: nSym
        for subid = 1: 48
            Rhsub(subid,:) = Rhsub(subid,:) + (hd(subid,sid) *  hp');
        end  
    end
end
Rhp = Rhp/NR_CH;
Rhsub = Rhsub/(NR_CH);
end


