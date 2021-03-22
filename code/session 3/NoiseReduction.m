function [S_mvdr_mwfL_stft,X_mvdr_mwfL_stft,N_mvdr_mwfL_stft] = NoiseReduction(y_STFT,x_STFT,n_STFT,...
                                                      SPP,TX_Mask,NMR,model,h_steer,gama)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
N_freqs = size(y_STFT,1);
N_frames = size(y_STFT,2);
num_mics = size(y_STFT,3);

% model.lambda = 0.995;                                                       % Forgetting factors for correlation matrices - can change
% model.SPP_thr = 0.95;                                                       % Threshold for SPP - can change
% model.mu = 1;
Rnn = cell(N_freqs,1);  Rnn(:) = {1e-9*randn(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {1e-9*randn(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values
Rxx = cell(N_freqs,1);  Rxx(:) = {1e-9*randn(num_mics,num_mics)};      % esitimated speech corr. matrix. Initialize to small random values                                                                                                                                       

% For MWF filter for left ear
S_mvdr_mwfL_stft = zeros(N_freqs,N_frames,num_mics);
X_mvdr_mwfL_stft = zeros(N_freqs,N_frames,num_mics);
N_mvdr_mwfL_stft = zeros(N_freqs,N_frames,num_mics);
W_mvdr_mwfL = (1/num_mics)*ones(num_mics,N_freqs,num_mics);
weight = (1/num_mics)*ones(num_mics,N_freqs,N_frames);
K = ones(N_freqs,num_mics,N_frames);alpha = 0.09;
% STFT Processing
% Looping through each time frame and each frequency bin

% tic
for l=1:N_frames % Time index
    for k = 2:N_freqs % Freq index
                
        % Create a vector of mic. signals
        Y_kl = squeeze(y_STFT(k,l,1:num_mics));  % M x 1 noisy mic sigs for this freq and time frame
        X_kl = squeeze(x_STFT(k,l,1:num_mics));
        N_kl = squeeze(n_STFT(k,l,1:num_mics));
        % ## Update the correlation matrices using the forgetting factor.
        % Threshold the SPP in order to distinguish between periods of speech and non-speech activity
        if(SPP(k,l)>model.SPP_thr)
             Ryy{k} = model.lambda^2*Ryy{k}+(1-model.lambda^2).* Y_kl* Y_kl';
            % speech + noise
        else
            % noise only
             Rnn{k} = model.lambda^2*Rnn{k}+(1-model.lambda^2).* Y_kl* Y_kl'; % N*N' is not available in real situation 
        end
      
       % calculate mu used in SDW MWF     
       % mu1: fix constant; mu2 spp model; mu3 psycho1;mu4 psycho2
%        model.mu = 1;
       mu = mu_calc(SPP(k,l),TX_Mask(k,l),NMR(k,l),model.mu);
       muxx(k,l) =  mu;  
       % calculate W_update in different model ,func =1:4
%        model.weight = 1; model.power = 1;
       [W_update,W_compare,mvdr_n1(k,l),mvdr_n2(k,l)] = weight_cal(Ryy{k}, Rnn{k},TX_Mask(k,l),gama{k},h_steer(:,k),model.power,model.weight,mu);  
       weight(:,k,l) =W_update;  
       weight_com(:,k,l) =W_compare; 
       % update w
         for m = 1:num_mics    
              % After experiments, we found that Filter parameter is not likely
              % larger than 2 if it converges.
            if(max(abs(W_update)) < 2.2) % This will skip some initial frames  
                W_mvdr_mwfL(:,k,m) =W_update;              
            end
            % Filtering the noisy speech, the speech-only, and the noise-only.
            S_mvdr_mwfL_stft(k,l,m) = W_mvdr_mwfL(:,k,m)'* Y_kl(1:num_mics);
            X_mvdr_mwfL_stft(k,l,m) = W_mvdr_mwfL(:,k,m)'* X_kl(1:num_mics);
            N_mvdr_mwfL_stft(k,l,m) = W_mvdr_mwfL(:,k,m)'* N_kl(1:num_mics);
        end
        
        
    end % end freqs
end % end time frames
% toc;
end

