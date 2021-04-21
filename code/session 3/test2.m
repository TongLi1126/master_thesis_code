clear;
close all
mex  ./pesq/*.c -output ./bin/PESQ_MEX
if ispc
    addpath('..\..\audio_files');
    addpath('..\..\sim_environment');
    addpath('..\..\sim_environment\RIR_generator_RIM');
    addpath('..\..\sim_environment\babblegen');
    addpath('..\evaluation');
    addpath('..\auxiliary function');
    addpath('.\pesq');
    addpath('.\bin');
else
    addpath('../../audio_files');
    addpath('../../sim_environment');
    addpath('../../sim_environment/RIR_generator_RIM');
    addpath('../../sim_environment/babblegen');
    addpath('../evaluation');
    addpath('../auxiliary function');
    addpath('./pesq');
    addpath('./bin');
end

%% Setup Configuration
% Here we set up the initial parameters for the simulations

%%%%% Initialize Simulation Parameters %%%%% 

num_mics = 2;                      % # mics
micsp = 0.03;               % mic spacing
room_dim = [4.3 6.9 2.6];   % room dimensions - typical living room - see https://www.nahb.org/en/research/housing-economics/special-studies/spaces-in-new-homes-2013.aspx
Rd = 0.13;                  % random displacement for randomized image method 
refmicpos = [2.1 3.3 1.7];  % ref mic (right most mic in room);
fs = 16000;                 % sampling frequency (Hz)
RIR_len = 8000;             % Length of RIR
cair = 340;                 % speed of sound in air (m/s)
siglen = 15;                % Length for generated signal (seconds)
siglen_samp = siglen*fs;    % Length for generated signal (samples)
source_angle = 0;           % All angles for which the source is to be defined (degrees, it is converted to radians in functions)
source_dist = 1;            % define source as 1m away
noise_angles = 90;          % angle for which there is a localised noise source.
noise_dist = 1*ones(size(noise_angles));            % noises at some distance away from source. will be scaled later to change SNR
% RT = 0.25;                   % reverberation time (s)
% SNR_in_amb = 2;             % desired SNR input - for all noises
SNR_in_uncorr = 30;         % desired SNR input for uncorrelated noise

% sim_param = struct('num_mics',num_mics,'micsp',micsp,'room_dim',room_dim,'Rd',Rd,'refmicpos',refmicpos,'fs',fs,'RIR_len',RIR_len,...
%                     'cair',cair,'siglen',siglen,'source_angle',source_angle,'noise_angles',noise_angles,'noise_dist',noise_dist,'RT',RT,...
%                     'SNR_in_amb',SNR_in_amb,'SNR_in_uncorr',SNR_in_uncorr);
% STFT parameter
nfft = 512;    % number of DFT points
window = sqrt(hann(nfft,'periodic')); % analysis window
noverlap = 2;   % factor for overlap. noverlap = 2 corresponds to 50% overlap
time = 0:1/fs:((size(siglen_samp,1)-1)/fs);
fs_RIR = fs;                
%% Load audio and noise signal
% Load audio signal
audiosignal{1} = 'speech1.wav';

% Load noise signal
noisesignal{1} = 'Babble_noise1.wav';
% Generating diffuse noise field instead - if desired
type_nf = 'spherical';    % Type of noise field: % 'spherical' or 'cylindrical'
% diff_noise = gen_babblesp_fcn(fs, cair, nfft, M, micsp, type_nf, siglen_samp);

RT_set = [0.3 0.6 0.9 1.5 2.5];
SNR_set = [-10 -5 ];
delta_SI_SNR = zeros(length(RT_set),length(SNR_set));
delta_SD = zeros(length(RT_set),length(SNR_set));
delta_PESQ = zeros(length(RT_set),length(SNR_set));
delta_STOI = zeros(length(RT_set),length(SNR_set));
%% Start loop
revi = 1;
for RT = RT_set
% Senario Generation
% For the scenario generation, we will use the Randomized Image Method
% (RIM). Here we specify the details pertaining to number/types of sources,
% noises, etc. 
% Generate convolved speech, noises, etc.
[speechIR, noiseIR, speech,speech_clean,n_local,m_pos,pc,s_pos,v_pos] = genmics_array_rim(audiosignal, noisesignal, source_angle, source_dist, noise_angles, noise_dist, room_dim, Rd, refmicpos, num_mics, micsp, fs, cair, RT, RIR_len, siglen);
speech_r = fftfilt(speechIR(800:end,:),speech_clean,siglen*fs);
% Setting and scaling the powers of all signals
     snrj= 1;
  for SNR_in_amb = SNR_set
sp_pow = var(speech(:,1));                      % Speech power at 1st mic
sf_uncorr = sqrt(sp_pow/10^(SNR_in_uncorr/10)); % Scale factor for desired SNR for uncorr noise, assume initial std. dev of uncorr noise = 1
n_uncorr = sf_uncorr*(randn(length(speech),num_mics));

n_amb = n_local;                                            % choose here ambient noise field - n_local and/OR diff_noise 
sf_amb = sqrt(sp_pow/(var(n_amb(:,1))*10^(SNR_in_amb/10))); % scale noise for desired SNR
n_amb_sc = sf_amb*n_amb;

% ## Obtain the noisy microphone signals 
noise = n_uncorr + n_amb_sc;
% Create noisy mic signals in the time domain:
noisy_sig = speech + noise ;  %  stacked microphone signals
snrb = snr(speech(:,1),n_amb_sc(:,1));
snrw = snr(speech(:,1),n_uncorr(:,1));
%% Apply WOLA analysis to observe signals in the STFT domain, Apply the SPP.

% ## Apply the WOLA analysis to the noisy mic. signals, the speech, and the noise.
%  reveberation
[y_STFT,f] = WOLA_analysis(noisy_sig,fs,window,nfft,noverlap) ;% To complete
[n_STFT,~] = WOLA_analysis(noise,fs,window,nfft,noverlap) ;% To complete
[x_STFT,~] = WOLA_analysis(speech,fs,window,nfft,noverlap) ;% To complete
[xe_STFT,~] = WOLA_analysis(speech_clean,fs,window,nfft,noverlap) ;% To complete
[xr_STFT,~] = WOLA_analysis(speech_r,fs,window,nfft,noverlap) ;% To complete
% Observe the STFT
clow = -60; chigh = 10; % lower and upper limits for signal power in spectrogram (can change for different resolutions)
[N_freqs, N_frames] = size(y_STFT(:,:,1));
% 
% figure;subplot(2,1,2);
% imagesc(time, f/1000, mag2db(abs(y_STFT(:,:,1))),[clow chigh]); colorbar;
% axis xy; set(gca,'fontsize', 14);
% set(gcf,'color','w'); xlabel('Time Frame'); ylabel('Frequency (kHz)');
% title('Spectrogram of reverberant noisy signal ')
% subplot(2,1,1);
% imagesc(time, f/1000, mag2db(abs(x_STFT(:,:,1))),[clow chigh]); colorbar;
% axis xy; set(gca,'fontsize', 14);
% set(gcf,'color','w'); xlabel('Time Frame'); ylabel('Frequency (kHz)')
% title('Spectrogram of speech ')
%% ## Compute the Speech Presence Probability on the reference microphone
[noisePowMat, SPP] =  spp_calc(speech(:,1),nfft,nfft/noverlap);
[noisePowMat_y, SPP_y] =  spp_calc(noisy_sig(:,1),nfft,nfft/noverlap);
% Create perfect VAD - use x_STFT (cheating...)
VAD_th = -60;   % VAD threshold
pVAD = ((mag2db(abs(x_STFT(5,:,1))))>VAD_th).';

% Observe the SPP
% figure; subplot(2,1,1);
% imagesc(1:N_frames, f,SPP); colorbar; axis xy; set(gcf,'color','w');
% set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('Speech Presence Probability for ref mic');
% subplot(2,1,2);
% imagesc(1:N_frames, f,SPP_y); colorbar; axis xy; set(gcf,'color','w');
% set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('Speech Presence Probability for ref non reverberant mic');
%% ## Calculate and observe TX_mask, TY_mask, TN_mask, NMR
TY_Mask = zeros(N_freqs,N_frames);
TX_Mask = zeros(N_freqs,N_frames);
TN_Mask = zeros(N_freqs,N_frames);
noise_spl = zeros(N_freqs,N_frames);
NMR = zeros(N_freqs,N_frames);
% 
for t = 1:N_frames
%  use 1st mic speech signal as masker
%      [TY_Mask(:,t),bark] = mask_cal(y_STFT(:,t,1),f,fs,nfft,t,1);
   [TX_Mask(:,t),bark] = mask_cal(x_STFT(:,t,1),f,fs,nfft,t,10);
   noise_spl(:,t) = psd2spl(n_STFT(:,t,1),nfft);
   NMR(:,t) = noise_spl(:,t)-TX_Mask(:,t);          
end

% for t = 1:N_frames
% %  use 1st mic noisy signal as masker
%    [TY_Mask(:,t),bark] = mask_cal(y_STFT(:,t,1),f,fs,nfft,t,10);
%    noise_spl(:,t) = psd2spl(n_STFT(:,t,1),nfft);
%    NMR(:,t) = noise_spl(:,t)-TX_Mask(:,t);          
% end

% % Observe the TM
% figure; subplot(2,1,1);
% imagesc(1:N_frames, f,TX_Mask); colorbar; axis xy; set(gcf,'color','w');
% set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('Masking threshold for ref mic');
% subplot(2,1,2);
% imagesc(1:N_frames, f, TY_Mask); colorbar;
% axis xy; set(gca,'fontsize', 14);
% set(gcf,'color','w'); xlabel('Time Frame'); ylabel('Frequency (Hz)'),title('Masking threshold for NON REV ref mic');
%% ##Implementation of the MWF

% % calculate RTF d and diffuse reverberation gama
[h_steer,gama] = RTF_cal(fs_RIR,m_pos,s_pos,num_mics,N_freqs) ;

Rnn = cell(N_freqs,1);  Rnn(:) = {1e-6*randn(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {1e-6*randn(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values
Rxx = cell(N_freqs,1);  Rxx(:) = {1e-6*randn(num_mics,num_mics)};      % esitimated speech corr. matrix. Initialize to small random values                                                                                                                                       

% Defination of model parameter                                            % Threshold for SPP - can change
model.lambda = 0.9;                                                       % Forgetting factors for correlation matrices - can change
model.SPP_thr = 0.60; % Threshold for SPP - can change
model.alpha = 0.33;

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
%              Rnn{k} = model.lambda*Rnn{k}+(1-model.lambda).* N_kl* N_kl'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% estimation of reverberant component%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        if(l>1);residual = rho_s (k,l-1); else;residual = 0;end
        rho_s(k,l)  = (1-model.alpha)*abs(xe_STFT(k,l,1)).^2 +model.alpha*residual;
        if(l>1);residual_r = rho_r (k,l-1);else;residual_r = 0;end
        rho_r (k,l) = (1-model.alpha)*abs(xr_STFT(k,l,1)).^2 +model.alpha*residual_r;
        
          
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
%             rho_s_esti(k,l) = Rxx_esti{k}(1,1)-rho_rev_realx(k,l);
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
          % mmse amplitude estimator
         sig_spectral = gamma(1.5)*(sqrt(v(k,l))/postSNR(k,l))*exp(-v(k,l)/2)...
         *((1+v(k,l))*besseli(0,v(k,l)/2) + v(k,l)*besseli(1,v(k,l)/2))*y_STFT(k,l,1);    
         sigPow(k,l) = sig_spectral*conj(sig_spectral);
     % MWF filter parameter
         model.rev =1;
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
%            W_mvdr_mwfL(:,k) = W_single_mvdr(k,l)* W_mvdr;   

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
           W_mvdr_mwfL(:,k) = W_single_pw(k,l) * W_mvdr;
                       
           % Filtering the noisy speech, the speech-only, and the noise-only.
           S_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* Y_kl(1:num_mics);
           X_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* X_kl(1:num_mics);
           N_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* N_kl(1:num_mics);
       


    end % end freqs
end % end time frames
% toc;

% Observe processed STFTst
% err = abs(x_STFT(:,:,1))-0.5*abs(S_mvdr_mwfL_stft);
% figure; subplot(2,1,1);
% imagesc(1:N_frames,f/1000,mag2db(abs(rho_r)),[-60 10]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), ...
%     title('residual noise');
% subplot(2,1,2); imagesc(1:N_frames,f/1000,mag2db(abs( rho_rev_realx)),[-60 10]); colorbar; axis xy; set(gcf,'color','w');
% set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('W_single_mvdr');
% 
% figure; subplot(2,1,1);
% imagesc(1:N_frames,f/1000,mag2db(abs(sigPow)),[-60 10]); colorbar; axis xy; set(gcf,'color','w');
% set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('microphne signal, 1st mic');
% S_L_enhanced = S_mvdr_mwfL_stft(:,:,1);
% subplot(2,1,2); imagesc(1:N_frames,f/1000,mag2db(abs(S_mvdr_mwfL_stft)),[-60 10]); colorbar; axis xy; 
% set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Enhanced Signal L - MWF');

% Apply the synthesis stage of the WOLA framework to obtain the time domain equivalents:

s_mwfL = WOLA_synthesis(S_mvdr_mwfL_stft,window,nfft,noverlap);% To complete (time-domain version of S_mvdr_mwfL_stft)
x_mwfL = WOLA_synthesis(X_mvdr_mwfL_stft,window,nfft,noverlap);% To complete (time-domain version of X_mvdr_mwfL_stft)
n_mwfL = WOLA_synthesis(N_mvdr_mwfL_stft,window,nfft,noverlap);% To complete (time-domain version of N_mvdr_mwfL_stft)


% PLOT SIGNALS
% figure;
% subplot(2,1,1);
% plot(speech); title('clean signal');
% subplot(2,1,2);
% plot(real(noisy_sig(:,1)));hold on; plot(real(s_mwfL(:,1)));title('signal after enhancement');

% LISTEN TO SIGNALS!
% soundsc(real(s_mwfL(:,1)),fs);
% soundsc(signal,fs);
%  soundsc(noisy_sig,fs_RIR);
% EVALUATION

% transform vad to time vad
pvad = repmat(pVAD',nfft,1); L = size(pvad,2);
vad = zeros((nfft/noverlap)*(L-1) + nfft, 1);
for m = 0:size(pvad,2)-1
    vad((floor(m * nfft/noverlap) + 1):floor(m * nfft/noverlap + nfft),:) = pvad(:,m+1)...
        +vad((floor(m * nfft/noverlap) + 1):floor(m *nfft/noverlap + nfft),:);
end
vad((vad >=1)) = 1;

enhance_sig = real(s_mwfL(:,1));
noisy_sig_ref = noisy_sig(:,1);
speech_ref = speech(:,1);
%  snr
delta_SI_SNR(revi,snrj) = SI_SNRatio(s_mwfL(:,1),n_mwfL(:,1),1, fs, vad)-SI_SNRatio(speech(:,1),noise(:,1),1, fs, vad);

% spectral distance SD
delta_SD(revi,snrj) =  distor_cal(x_STFT,S_mvdr_mwfL_stft,nfft,TX_Mask)-distor_cal(x_STFT,y_STFT,nfft,TX_Mask);
% PESQ
delta_PESQ(revi,snrj) = pesq_mex(speech_ref, enhance_sig, fs, 'wideband')-pesq_mex(speech_ref, noisy_sig_ref, fs, 'wideband')  ;
% STOI short time objective Intelligibility Measure
delta_STOI(revi,snrj) = stoi(speech_ref, enhance_sig, fs)-stoi(speech_ref,noisy_sig_ref, fs) ;

snrj = snrj +1;
 end
revi = revi +1;
end
%% Plot results
figure;
for line = 1:length(SNR_set)
    
plot(RT_set,delta_SI_SNR(:,line),'-s'); hold on
end
legend('SNR =-10','SNR =-5','SNR =0','SNR =5','SNR =10','SNR =20');grid on;ylabel('\Delta SI-SNR');xlabel('Reverberation time T_{60}')

figure;
for line = 1:length(SNR_set)
plot(RT_set,delta_SD(:,line),'-s'); hold on
end
legend('SNR =-10','SNR =-5','SNR =0','SNR =5','SNR =10','SNR =20');grid on;ylabel('\Delta SD');xlabel('Reverberation time T_{60}')

figure;
for line = 1:length(SNR_set)
plot(RT_set,delta_PESQ(:,line),'-s'); hold on
end
legend('SNR =-10','SNR =-5','SNR =0','SNR =5','SNR =10','SNR =20');grid on;ylabel('\Delta PESQ');xlabel('Reverberation time T_{60}')

figure;
for line = 1:length(SNR_set)
plot(RT_set,delta_STOI(:,line),'-s'); hold on
end
legend('SNR =-10','SNR =-5','SNR =0','SNR =5','SNR =10','SNR =20');grid on;ylabel('\Delta STOI');xlabel('Reverberation time T_{60}')




