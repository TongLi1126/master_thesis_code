function [speechIR, noiseIR, speech,speech_clean, n_uncorr,n_local,pm,pc,ps_const,s_pos,v_pos, VAD] = gen_micsignals_rim_array(dc_angle, dc_dist, source_angle, source_dist, noise_angles, noise_dist, room_dim, Rd, refmicpos, M, micsp, fs_RIR, c,RT)
%28/03/2016 %
% OLDER - see genmics_array_rim.
%
%   function [speechIR, noiseIR, speech, speech_clean, n_uncorr,n_local_init,pm,pc,ps_const,s_pos,v_pos, VAD] = gen_micsignals_rim_array(dc_angle, dc_dist, source_angle(z), source_dist, noise_angles, noise_dist, room_dim,Rd, refmicpos, M, micsp, fs_RIR, cair,RT);
%
% Inputs:
% dc_angle, dc_dist  - angle and distances from the constraint direction to
%                      the centre of array (Angles in degrees!)
% source angle, dist - angle and distances from the source to centre of
%                      array. - Angles in degrees!
% noise angle, dist  - angle and distances from the noise to centre (can
%                      specify multiple angles, distances for noises)
%                      Angles in degrees!
% room_dim - dimensions of room - 1x3 = [x y z]
% Rd - random displacement for randomized IM
% refmicpos - position of first mic 1x3 = [x y z]
% M - number of mics
% micsp - mic spacing
% fs_RIR - sampling frequency
% c - speed of sound (m/s)
% RT - Reverberation Time (seconds)
%
% Outputs:
% SpeechIR - impulse response from speech signal
% NoiseIR - impulse response from noise signal
% speech - convolved speech with RIR
% speech_clean - raw speech with no convolution, actual clean speech
% n_uncorr - uncorrelated noise (mic noise), spatially white
% n_local - correlated noise (babble talk, etc)
% pm - dimensions of the mic array
% pc - centre of array
% ps_const - position of the constraint direction
% s_pos - position of source
% v_pos - position of noises
% vad - vad used in detection (set of 1s and 0's) - see below for threshold
%
% NOTE : DID NOT THROW ERROR IF COORDINATES DO NOT MATCH WITH ROOM
% DIMENSIONS, DOUBLE CHECK THIS ON YOUR OWN!!!
%
%UPDATE - This uses Randomized IM - see DeSena et al..

%clear all

% In generating the microphone signals, make sure to think about the room
% scenario and that there are no conflicting coordinates.
% Confirmed that this method works with the delay and sum BF
% It also works with MUSIC DOA estimation - however it is affected by
% reverberation and signal to noise ratios.

v_postr = zeros(length(noise_angles),3);

%convert angles to radians!
source_angle = (pi/180)*source_angle;
noise_angles = (pi/180)*noise_angles;
dc_angle = (pi/180)*dc_angle;


%For sources (we assume one):
[xt, yt] = pol2cart(source_angle, source_dist);     % distance to add to the centre of the array.
[pm, pc] = gen_micarray(refmicpos, micsp, M);        % gen array and get centre
zt = 0;                                             % z-coord offset from centre mic (=0 since we are in 2D)
ptemp = [xt yt zt].';                                
s_pos = ptemp + pc;                                 % translate the coordinates of the source
s_postr = s_pos'; %transposes
pmtr = pm';         

%For constraint:
[xtc ytc] = pol2cart(dc_angle, dc_dist);    
ptempc = [xtc ytc zt].';
ps_const = ptempc + pc;     % Direction of the constraint

%For noises (no need to regenerate array, only coords for noise):
for g=1:length(noise_angles)
    [xtn, ytn] = pol2cart(noise_angles(g), noise_dist(g));
    ptempn = [xtn ytn zt].';                                
    v_pos = ptempn + pc;                                  %translate the coordinates of the source
    v_postr(g,:) = v_pos.';
end 

rev_time = RT;                         % Reverberation time (s)
Nsamp = 8000;                           % Number of samples
%mtype = 'omnidirectional';              % Type of microphone - other options: `subcardioid',`cardioid', `hypercardioid', `bidirectional','omnidirectional'
%order = -1;                             % -1 equals maximum reflection order!
%dim = 3;                                % Room dimension
%orientation = 0;                        % Microphone orientation (rad) - default is to point in dir of +ve x axis, so leave as is.
%hp_filter = 1;                          % Disable high-pass filter
Tw = 20;                                % samples of Low pass filter 
Fc = 0.9;                               % cut-off frequency
Nrefl =  [ 0;0;0];                      % Reflection order

% Generate RIR_sources - i.e. the RIRs for the source signals
numsource = length(s_postr(:,1));         % get the number of source signals we have
RIR_sources = zeros(Nsamp, M, numsource);
for n=1:numsource
    % Randomized IM
    [htemps, Sr]  = ISM(pmtr',s_postr(n,:),room_dim,rev_time,Nrefl,Nsamp, Rd, [],Tw,Fc,fs_RIR,c);
    RIR_sources(:,:,n) = htemps;
    
    % IM Habets
    %htempsource = rir_generator(c, fs_RIR, pmtr, s_postr(n,:), room_dim, rev_time, Nsamp, mtype, order, dim, orientation, hp_filter);
    %RIR_sources(:,:,n) = htempsource.';
end
speechIR = RIR_sources;

% Generate RIR_noise - i.e. the RIRs for the noise signals
numnoise = length(v_postr(:,1));         % get the number of noise signals we have
RIR_noise = zeros(Nsamp, M, numnoise);
for p=1:numnoise
    
    [htempn, Sr]  = ISM(pmtr',v_postr(p,:),room_dim,rev_time,Nrefl,Nsamp, Rd, [],Tw,Fc,fs_RIR,c);
    RIR_noise(:,:,p) = htempn; 
    
    %htempn = rir_generator(c, fs_RIR, pmtr, v_pos(p,:), room_dim, rev_time, Nsamp, mtype, order, dim, orientation, hp_filter);
    %RIR_noise(:,:,p) = htempn.';
end
noiseIR = RIR_noise;

%setup names of source and noise files.
%these will be convolved with the IR to get the sound
%You can add more source files or noise files as desired.
%audiosignal will contain all the audio signals
%noisesignal will contain all the noise signals

%audiosignal{1} = 'speech1_hint_trun.wav';
%audiosignal{1} = 'cleanspeech_sg.wav';
audiosignal{1} = 'speech5_trun_hint16k.wav';    %longer signals
%audiosignal{2} = 'speech1_hint_trun_16k.wav';
%audiosignal{1} = 'speech2.wav';
%audiosignal{3} = 'White_noise1.wav';
%audiosignal{4} = 'speech1.wav';

noisesignal{1} = 'Multitalker_track3_44k_longer1.wav';  %longer noises
noisesignal{2} = 'Multitalker_track3_44k_longer2.wav';  
noisesignal{3} = 'Multitalker_track3_44k_longer3.wav';  
% noisesignal{1} = 'Multitalker_track3_10.wav';
% noisesignal{2} = 'Multitalker_track3_20.wav';
% noisesignal{3} = 'Multitalker_track3_30.wav';
% noisesignal{4} = 'Multitalker_track3_40.wav';
% noisesignal{5} = 'Multitalker_track3_50.wav';
% noisesignal{6} = 'Multitalker_track3_60.wav';
% noisesignal{7} = 'Multitalker_track3_70.wav';
% noisesignal{8} = 'Multitalker_track3_80.wav';
% noisesignal{9} = 'Multitalker_track3_90.wav';
% noisesignal{10} = 'Multitalker_track3_100.wav';
% noisesignal{11} = 'Multitalker_track3_110.wav';
% noisesignal{12} = 'Multitalker_track3_120.wav';
% noisesignal{13} = 'Multitalker_track3_130.wav';


%noisesignal{2} = 'carnoise.wav';
%noisesignal{2} = 'factorynoise.wav';
%noisesignal{3} = 'speechnoise.wav';

%DISCLAIMER: THE LENGTHS OF THE VECTORS IN audiosignal and noisesignal
%should be the same. Can pad with zeros. But do this manually.

%if fs_RIR ~= 44100
%    error ('RIR sample freq should be 44.1k!')
%end

siglen = 20;                %Set here the total time of the signal desired.
Ntotal = siglen*fs_RIR;
dt = 1/fs_RIR;
t = 0:dt:(Ntotal-1)*dt;

%Check to see how many microphone signals, noise sources and audio sources
%exist. Audio sources are the sources as defined from the room GUI, similar
%as noise sources. These are independent of the audio and noise signals.
%The sources just refer to the excitation points. The signals refer to the
%actual signals that will be used for excitation.

%M = length(RIR_sources(1,:,1)); %#microphones
%NAS = length(RIR_sources(1,1,:)); %#audio sources

%---LOADING SIGNALS---%
clear AS;
for h = 1:numsource
    [asignal,fs_as] = audioread(audiosignal{h});    %read in audio source
    if fs_as ~= fs_RIR                              %checking the sampling rate
        asignal = resample(asignal,fs_RIR,fs_as);   %resample the data
    end
    asignal = asignal(1:Ntotal);
    asignal = asignal/max(asignal);
    AS(:,h) = asignal;                              %create a vector whose colums are the audio signals
end
speech_clean = AS;

%Sometimes there may not be a localized noise signal, so check for this:
NSvar = isempty(RIR_noise); %check to see if variable is empty
clear NS;
if NSvar == 1 %if the variable is empty
    numnoise = 0;
else
    %NNS = length(RIR_noise(1,1,:)); %#noise sources
    %Load noise signals:
    for k = 1:numnoise
        [nsignal,fs_ns] = audioread(noisesignal{k}); %read in audio source
        if fs_ns ~= fs_RIR 
            nsignal = resample(nsignal,fs_RIR,fs_ns); 
        end
        nsignal = nsignal(1:Ntotal);
        nsignal = nsignal/max(nsignal);             %normalize for consistent scaling later on
        NS(:,k) = nsignal;                           %create a vector whose colums are the noise signals
    end
end



%---CONVOLUTIONS---%
%Creating the Speech only section of the signals:
%The matrix "speech" contains the clean speech signal as received by the
%microphones.

speech_init = zeros(Ntotal,M);      %initial speech as acquired directly from the scenario. This will serve as a reference later on.
for m=1:M %loop over microphones
   for k=1:numsource %loops over number of audio sources
       speech_init(:,m)=speech_init(:,m)+fftfilt(RIR_sources(:,m,k),AS(:,k));
   end
end

%scaled speech to use for controlling the input SNR.
speech = speech_init;

%Creating Noise only sections
%There are 2 types of noise: 
%1. Localized noise source from the GUI - n_local
%2. Uncorrelated noise added to each microphone - n_uncorr
% Total noise = n_local + n_uncorr
% We change the amplitude of n_local to vary the SNR!

n_local_init = zeros(Ntotal,M);
if numnoise > 0
    for m=1:M 
        for k=1:numnoise 
            n_local_init(:,m)=n_local_init(:,m)+fftfilt(RIR_noise(:,m,k),NS(:,k));
        end
    end
end

%Scale localized noise source here to change SNR
n_local = n_local_init;

%Use ideal VAD below to select speech only sections:
VADthresh = 1e-8;                                   %threshold for VAD
VAD = abs(speech(:,1))>std(speech(:,1))*VADthresh;  %To identify when speech is active

%speechpow_init = var(speech_init(VAD==1,1));        %power of the initial speech ONLY when the VAD is on.
speechpow = var(speech(VAD==1,1));                  %power of the scaled speech ONLY when the VAD is on.

n_uncorr = zeros(Ntotal,M);
n_uncorrpow = 0.01*speechpow;                        %power of the uncorrelated noise in relation to speech (can change fraction, 0.05 = 5% of speech pow)
for m = 1:M
    n_uncorr(:,m) = sqrt(n_uncorrpow)*(randn(Ntotal,1)); % generate separate uncorrelated noise.
end

noise = n_local + n_uncorr;                         %total noise contribution
noisepow_loc = var(n_local(:,1));
noisepow_uncorr = var(n_uncorr(:,1));               %individual noise contributions

noisepow = noisepow_loc + noisepow_uncorr;          %total noise power
%noisepowchk = var(noise(:,1));                      %2nd check for noise power

mic = speech + noise;                               %total microphone signals with speech and noise

SNRin = 10*log10(speechpow/noisepow);


