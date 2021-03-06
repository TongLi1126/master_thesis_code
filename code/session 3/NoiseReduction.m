function [S_mvdr_mwfL_stft,X_mvdr_mwfL_stft,N_mvdr_mwfL_stft] = NoiseReduction(y_STFT,x_STFT,n_STFT,...
                                          xe_STFT,xr_STFT,SPP,TX_Mask,model,h_steer,gama)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
N_freqs = size(y_STFT,1);
N_frames = size(y_STFT,2);
num_mics = size(y_STFT,3);


Rnn = cell(N_freqs,1);  Rnn(:) = {1e-6*randn(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {1e-6*randn(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values
Rxx = cell(N_freqs,1);  Rxx(:) = {1e-6*randn(num_mics,num_mics)};      % esitimated speech corr. matrix. Initialize to small random values                                                                                                                                       

% Defination of model parameter                                            % Threshold for SPP - can change
% model.lambda = 0.90;                                                       % Forgetting factors for correlation matrices - can change
% model.SPP_thr = 0.60; % Threshold for SPP - can change
% model.alpha = 0.33;

% For MWF filter for left ear
S_mvdr_mwfL_stft = zeros(N_freqs,N_frames);
X_mvdr_mwfL_stft = zeros(N_freqs,N_frames);
N_mvdr_mwfL_stft = zeros(N_freqs,N_frames);
noise_mvdr = zeros(N_freqs,N_frames);
W_single_mvdr = zeros(N_freqs,N_frames);
W_single_pw = zeros(N_freqs,N_frames);
% For estimation reverberant power
rho_s = zeros(N_freqs,N_frames);
rho_r = zeros(N_freqs,N_frames);v = zeros(N_freqs,N_frames);
% For estimation of speech power
postSNR = zeros(N_freqs,N_frames);
prioiSNR = zeros(N_freqs,N_frames);
sigPow = zeros(N_freqs,N_frames);
W_mvdr_mwfL = (1/num_mics)*ones(num_mics,N_freqs);

% STFT Processing
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
             Ryy{k} = model.lambda*Ryy{k}+(1-model.lambda).* Y_kl* Y_kl';
            % speech + noise
        else
            % noise only
             Rnn{k} = model.lambda*Rnn{k}+(1-model.lambda).* Y_kl* Y_kl'; % N*N' is not available in real situation 
        end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% estimation of reverberant component%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if(l>1);residual = rho_s (k,l-1); else;residual = 0;end
        rho_s(k,l)  = (1-model.alpha)*abs(xe_STFT(k,l,1)).^2 +model.alpha*residual;
        if(l>1);residual_r = rho_r (k,l-1);else;residual_r = 0;end
        rho_r (k,l) = (1-model.alpha)*abs(xr_STFT(k,l,1)).^2 +model.alpha*residual_r;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% estimation of speech power component%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%           noise_power = abs(n_STFT(k,l,1))^2;
%           noise_power = noisePowMat_y(k,l);
            noise_power = Rnn{k}(1,1);
          postSNR(k,l) = abs(y_STFT(k,l,1))^2/noise_power;
         if l==1 ;resi = 0;else ;resi = sigPow(k,l-1)/noise_power; end
          a = 0.6;  
          prioiSNR(k,l) = a * resi +(1-a)*max(0,postSNR(k,l)-1);
          v(k,l) =  (prioiSNR(k,l)/(1+prioiSNR(k,l)))*postSNR(k,l);
          v(k,l)  = min(500,v(k,l));
          % mmse
         sig_spectral = gamma(1.5)*(sqrt(v(k,l))/postSNR(k,l))*exp(-v(k,l)/2)...
         *((1+v(k,l))*besseli(0,v(k,l)/2) + v(k,l)*besseli(1,v(k,l)/2))...
         *y_STFT(k,l,1);
     
         sigPow(k,l) = sig_spectral*conj(sig_spectral);
     % MWF filter parameter
%          model.rev = 1;
         rho_rev =   model.rev*rho_r(k,l) ;
         rho_ss  =  abs(sigPow(k,l) ); h = h_steer(:,k); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%MWF filter    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
              
           Rrr = rho_rev *gama{k};         
           noise_mvdr(k,l) = 1/(abs(h'*((Rrr + Rnn{k})\h)));
           W_single_mvdr(k,l) = rho_ss/(rho_ss + noise_mvdr(k,l));
           W_single_mvdr(k,l) = max(0, W_single_mvdr(k,l));
           W_mvdr = 2*(1/abs(h'*((Rnn{k}+Rrr)\h)))*((Rnn{k}+Rrr)\h);
           if ((l<24))
              W_mvdr = 0;
           end 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%perceptual weighting   %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
           xi = 0.005;           
           psd_TX = (10^((TX_Mask(k,l)-90)/10))* 512*512/2;
           W_single_pw(k,l) = (xi + sqrt(psd_TX*0.1/noise_mvdr(k,l)));           
           if (abs(W_single_pw(k,l) ) >5 && abs(W_single_pw(k,l) ) < 60)
              W_single_pw(k,l) = 1 ;
           elseif (abs(W_single_pw(k,l) ) > 60)
              W_single_pw(k,l) = 0;
           end          
          
                                 
          if (model.weight==1)
               W_mvdr_mwfL(:,k) = W_single_mvdr(k,l)* W_mvdr;   
          else
               W_mvdr_mwfL(:,k) = W_single_pw(k,l) * W_mvdr;
          end
           
           % Filtering the noisy speech, the speech-only, and the noise-only.
           S_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* Y_kl(1:num_mics);
           X_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* X_kl(1:num_mics);
           N_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* N_kl(1:num_mics);
       


    end % end freqs
end % end time frames
% toc;
