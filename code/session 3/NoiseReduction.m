function [S_mvdr_mwfL_stft,X_mvdr_mwfL_stft,N_mvdr_mwfL_stft] = NoiseReduction(y_STFT,x_STFT,n_STFT,...
                                          xe_STFT,xr_STFT,SPP,TX_Mask,NMR,model,h_steer,gama)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
N_freqs = size(y_STFT,1);
N_frames = size(y_STFT,2);
num_mics = size(y_STFT,3);

Rnn = cell(N_freqs,1);  Rnn(:) = {1e-6*randn(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {1e-6*randn(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values
Rxx = cell(N_freqs,1);  Rxx(:) = {1e-6*randn(num_mics,num_mics)};      % esitimated speech corr. matrix. Initialize to small random values                                                                                                                                       
Rxx_real = cell(N_freqs,1);  Rxx_real(:) = {1e-9*randn(num_mics,num_mics)};      % esitimated speech corr. matrix. Initialize to small random values                                                                                                                                       
Rxx_esti = cell(N_freqs,1);  Rxx_esti(:) = {1e-9*randn(num_mics,num_mics)};      % esitimated speech corr. matrix. Initialize to small random values                                                                                                                                       
Rnn_real = cell(N_freqs,1);  Rnn_real(:) = {1e-9*randn(num_mics,num_mics)};      % esitimated speech corr. matrix. Initialize to small random values                                                                                                                                       

% model.lambda = 0.8;                                                       % Forgetting factors for correlation matrices - can change
% model.SPP_thr = 0.5;                                                       % Threshold for SPP - can change

model.alpha = 0.33;
% For MWF filter for left ear
S_mvdr_mwfL_stft = zeros(N_freqs,N_frames,num_mics);
X_mvdr_mwfL_stft = zeros(N_freqs,N_frames,num_mics);
N_mvdr_mwfL_stft = zeros(N_freqs,N_frames,num_mics);
W_mvdr_mwfL = (1/num_mics)*ones(num_mics,N_freqs,num_mics);
weight = (1/num_mics)*ones(num_mics,N_freqs,N_frames);
K = ones(N_freqs,num_mics,N_frames);alpha = 0.09;
rho_s = zeros(N_freqs,N_frames);
rho_r = zeros(N_freqs,N_frames);v = zeros(N_freqs,N_frames);
% STFT Processing
% Looping through each time frame and each frequency bin

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
%              Ryy{k} = Rnn{k};
        end           
            Rxx_real{k} = (model.alpha*Rxx_real{k}+(1-model.alpha).* X_kl* X_kl');
            Rnn_real{k} = (model.alpha*Rnn_real{k}+(1-model.alpha).* N_kl* N_kl');
            Rn1(k,l) = Rnn_real{k}(1,1);Rn2(k,l) = Rnn{k}(1,1);
        % model.Rx choose which model to estimate Rx, 1: real update 
        % 2: direct estimate 3: decesion directed 
           if (abs(Ryy{k}(1,1))> abs(Rnn{k}(1,1)))
           Rxx_esti{k} = Ryy{k};
           else
           Rxx_esti{k} = 1e-6*Rnn{k}(1,1);
           end
           Rxx{k} = Rxx_esti{k};
           [vector,eigen] = eig(Rxx{k});
           eigen(find(eigen<0)) = 0;
           Rxx_remove0{k} = pinv(vector')* eigen* pinv(vector);
           L = chol(gama{k});
           [~,lamada] = eig(inv(L) *  Rxx_remove0{k} * inv(L'),'vector');
          
           if k<81
               rho_rev_realx(k,l) = abs(lamada(2)/300) ;
           else
               rho_rev_realx(k,l) = abs(lamada(2));
           end
            if(l>1);residual = rho_s (k,l-1); else;residual = 0;end
            rho_s(k,l)  = (1-model.alpha)*abs(xe_STFT(k,l,1)).^2 +model.alpha*residual;
            if(l>1);residual_r = rho_r (k,l-1);else;residual_r = 0;end
            rho_r (k,l) = (1-model.alpha)*abs(xr_STFT(k,l,1)).^2 +model.alpha*residual_r;
        
                a = 0.6;  
              noise_power = abs(n_STFT(k,l,1))^2;
        %       noise_power = noisePowMat_y(k,l);
              postSNR(k,l) = abs(y_STFT(k,l,1))^2/noise_power;
             if l==1 ;resi = 0;else ;resi = Rx3(k,l-1)/noise_power; end
              prioiSNR(k,l) = a * resi +(1-a)*max(0,postSNR(k,l)-1);
              v(k,l) =  (prioiSNR(k,l)/(1+prioiSNR(k,l)))*postSNR(k,l);
              v(k,l)  = min(100,v(k,l));
              % mmse
%              sig_spectral = gamma(1.5)*(sqrt(v(k,l))/postSNR(k,l))*exp(-v(k,l)/2)...
%              *((1+v(k,l))*besseli(0,v(k,l)/2) + v(k,l)*besseli(1,v(k,l)/2))...
%              *y_STFT(k,l,1);           
              sig_spectral = (prioiSNR(k,l)/(1+prioiSNR(k,l))) *y_STFT(k,l,1);
              
              Rx3(k,l) = sig_spectral*conj(sig_spectral);
             rho_rev =   rho_r(k,l) ;
%                rho_rev =  0;
              rho_ss  =  abs(Rx3(k,l) ); h = h_steer(:,k);Rn =  Rnn{k};
%        calculate mu used in SDW MWF     
%        mu1: fix constant; mu2 spp model; mu3 psycho1;mu4 psycho2
           model.mu = 1;
           mu = mu_calc(SPP(k,l),TX_Mask(k,l),NMR(k,l),model.mu);
           muxx(k,l) =  mu;  
           Rrr = rho_rev *gama{k};           
           W_single = rho_ss/(rho_ss + mu /(abs(h'*pinv(Rrr + Rn)*h)));
           W_single = max(0,W_single );
           W_mvdr = (1/abs(h'*pinv(Rn+Rrr)*h))*pinv(Rn+Rrr)*h;
%            W_update =  W_single * W_mvdr; 
% [W_update] = weight_cal(Ryy{k}, Rnn{k},Rxx{k},rho_rev,rho_ss,TX_Mask,gama,h,1,1);
 [W_update] = weight_cal(Ryy{k}, Rnn{k},Rxx{k},rho_rev,rho_ss,TX_Mask(k,l),gama{k},...
     h,model.weight,model.mu);

  
       weight(:,k,l) =W_update;  
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
function mu = mu_calc(SPP,TX_Mask,NMR,func)

switch func
    
    case 1 % mu is constanct
        mu = 1;
    case 2   
        alpha = 0.75;     % spp weighted  mu_spp = 1/(alpha/mu +(1-alpha)*SPP(k,l));
        mu_spp = 2;
        mu = 1/(alpha/mu_spp +(1-alpha)*SPP);
    case 3 % mu based on psychoacoustic model T1s(1st model)
        a = 4.374; beta = -0.0282;v = 40;
        mu = (a*exp(beta*TX_Mask))*(TX_Mask<=v);
    case 4 % mu based on 2nd psychoacoustic model
        gama = 0.3225; delta = 0.9; epsilon = 0.97;
        if NMR>0
        mu = (gama*NMR^delta+epsilon);
        else 
            mu = 0;
        end
    case 5
        if NMR>0
        mu = 4*exp(0.04*NMR);
        else 
            mu = 0;
        end       

end      
end
