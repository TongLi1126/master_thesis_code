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
source_angle = 0;           % All angles for which the source is to be defined (degrees, it is converted to radians in functions)
source_dist = 1;            % define source as 1m away
noise_angles = 90;          % angle for which there is a localised noise source.
noise_dist = 1*ones(size(noise_angles));            % noises at some distance away from source. will be scaled later to change SNR
% RT = 0.25;                   % reverberation time (s) Loop and compare
% SNR_in_amb = 2;             % desired SNR input - for all noises Loop and compare
SNR_in_uncorr = 30;         % desired SNR input for uncorrelated noise

% sim_param = struct('num_mics',num_mics,'micsp',micsp,'room_dim',room_dim,'Rd',Rd,'refmicpos',refmicpos,'fs',fs,'RIR_len',RIR_len,...
%                     'cair',cair,'siglen',siglen,'source_angle',source_angle,'noise_angles',noise_angles,'noise_dist',noise_dist,'RT',RT,...
%                     'SNR_in_amb',SNR_in_amb,'SNR_in_uncorr',SNR_in_uncorr);
      
%% Load audio and noise signal
% Load audio signal
audiosignal{1} = 'speech1.wav';

% Load noise signal
noisesignal{1} = 'Babble_noise1.wav';
% Generating diffuse noise field instead - if desired
type_nf = 'spherical';    % Type of noise field: % 'spherical' or 'cylindrical'
% diff_noise = gen_babblesp_fcn(fs, cair, nfft, M, micsp, type_nf, siglen_samp);

%% LOOP RT AND SNR
RT_set = [0.3 0.6 0.9 1.5 2.5];
SNR_set = [10 15 20];
delta_SI_SNR = zeros(length(RT_set),length(SNR_set)+1);
delta_SD = zeros(length(RT_set),length(SNR_set)+1);
delta_PESQ = zeros(length(RT_set),length(SNR_set)+1);
delta_STOI = zeros(length(RT_set),length(SNR_set)+1);
% Generate TX standard
model.rev = 1;
model.weight = 1;
model.mask = 1;
% Start loop
revi = 1;
for RT = 0.61
  snrj= 1;
  for SNR_in_amb = SNR_set
      [delta_SI_SNR(revi,snrj),delta_SD(revi,snrj),delta_PESQ(revi,snrj),delta_STOI(revi,snrj),SI_SNR(revi,snrj),...
        SD(revi,snrj),PESQ(revi,snrj),STOI(revi,snrj)] = ...
      SpeechEnhance (model,audiosignal, noisesignal, source_angle, source_dist, noise_angles, noise_dist,num_mics, micsp, RT,SNR_in_uncorr,SNR_in_amb );
      snrj = snrj +1;
  end
  revi = revi +1;
end

%% compare performance with TS_enhanced weight 3; 
% % compare performance with MWF weight 2; 
% % compare performance with TY,mask ==0 
model.rev = 1;
model.weight = 3;
% % Start loop
revi = 1;
for RT = RT_set
  snrj= 2;
  for SNR_in_amb = 20
      [delta_SI_SNR(revi,snrj),delta_SD(revi,snrj),delta_PESQ(revi,snrj),delta_STOI(revi,snrj),SI_SNR(revi,snrj),...
        SD(revi,snrj),PESQ(revi,snrj),STOI(revi,snrj)] = ...
      SpeechEnhance (model,audiosignal, noisesignal, source_angle, source_dist, noise_angles, noise_dist,num_mics, micsp, RT,SNR_in_uncorr,SNR_in_amb );
      snrj = snrj +1;
  end
  revi = revi +1;
end


%% Plot results

%DELTA

illustrate(RT_set,SNR_set,delta_SI_SNR,'\Delta SI-SNR');
illustrate(RT_set,SNR_set,delta_SD,'\Delta SD');
illustrate(RT_set,SNR_set,delta_PESQ,'\Delta PESQ');
illustrate(RT_set,SNR_set,delta_STOI,'\Delta STOI');
% Absolute OUTPUT
illustrate(RT_set,SNR_set,SI_SNR,'SI-SNR');
illustrate(RT_set,SNR_set,SD,'SD');
illustrate(RT_set,SNR_set,PESQ,'PESQ');
illustrate(RT_set,SNR_set,STOI,'STOI');



