% 20/11/2017 R. Ali
% This code utilizes the Randomized Image Method to get
% the appropriate signals.

function [speechIR, noiseIR, speech, speech_clean, n_local, pm, pc, s_pos, v_pos] = genmics_array_rim(audiosignal, noisesignal, source_angle, source_dist, noise_angles, noise_dist, room_dim, Rd, refmicpos, M, micsp, fs_RIR, cair, RT, RIR_len, siglen)
%
%
%  function [speechIR, noiseIR, speech,speech_clean,n_local,pm,pc,s_pos,v_pos] = genmics_array_rim(audiosignal, noisesignal, source_angle, source_dist, noise_angles, noise_dist, room_dim, Rd, refmicpos, M, micsp, fs_RIR, cair, RT, Nsamp, siglen)%
%
%
%
% Inputs:
% audiosignal        - struct of the wav file names
% noisesignal        - struct of the wav file names%
% source angle, dist - angle and distances from the source to centre of
%                      array. - Angles in degrees!
% noise angle, dist  - angle and distances from the noise to centre (can
%                      specify multiple angles, distances for noises)
%                      Angles in degrees!
% room_dim          - dimensions of room - 1x3 = [x y z]
% Rd                - random displacement for randomized IM
% refmicpos         - position of first mic 1x3 = [x y z]
% M                 - number of mics
% micsp             - mic spacing
% fs_RIR            - sampling frequency
% cair              - speed of sound in air (m/s)
% RT                - Reverberation Time (seconds)
% RIR_len           - Number of samples for RIR length
% siglen            - Length of the signal (seconds)
%
% Outputs:
% SpeechIR  - impulse response from speech signal
% NoiseIR   - impulse response from noise signal
% speech    - convolved speech with RIR
% speech_clean  - raw speech with no convolution, actual clean speech
% n_uncorr      - uncorrelated noise (mic noise), spatially white
% n_local       - correlated noise (babble talk, etc)
% pm            - positions of mics in array
% pc            - centre of array
% s_pos     - position of source
% v_pos     - position of noises
%
% NOTE : DID NOT THROW ERROR IF COORDINATES DO NOT MATCH WITH ROOM
% DIMENSIONS, DOUBLE CHECK THIS ON YOUR OWN!!!
% In generating the microphone signals, make sure to think about the room
% scenario and that there are no conflicting coordinates.



%convert angles to radians!
source_angle = (pi/180)*source_angle;
noise_angles = (pi/180)*noise_angles;

%For sources (we assume one):
[xt, yt] = pol2cart(source_angle, source_dist);     % distance to add to the centre of the array.
[pm, pc] = gen_micarray(refmicpos, micsp, M);       % gen array and get centre
zt = 0;                                             % z-coord offset from centre mic (=0 since we are in 2D)
ptemp = [xt yt zt].';                                
s_pos = ptemp + pc;                                 % translate the coordinates of the source

%For noises (no need to regenerate array, only coords for noise):
for g=1:length(noise_angles)
    [xtn, ytn] = pol2cart(noise_angles(g), noise_dist(g));
    ptempn = [xtn ytn zt].';                                
    v_pos(:,g) = ptempn + pc;                                  %translate the coordinates of the source
end 

rev_time = RT;                         % Reverberation time (s)
Tw = 20;                                % samples of Low pass filter 
Fc = 0.9;                               % cut-off frequency
Nrefl =  [ 0;0;0];                      % Reflection order

% Generate RIR_sources - i.e. the RIRs for the source signals
numsource = length(s_pos(1,:));         % get the number of source signals we have
RIR_sources = zeros(RIR_len, M, numsource);
for n=1:numsource
    % Randomized IM
    [htemps, Sr]  = ISM(pm,s_pos(:,1)',room_dim,rev_time,Nrefl,RIR_len, Rd, [],Tw,Fc,fs_RIR,cair);
    RIR_sources(:,:,n) = htemps;
end
speechIR = RIR_sources;

% Generate RIR_noise - i.e. the RIRs for the noise signals
numnoise = length(v_pos(1,:));         % get the number of noise signals we have
RIR_noise = zeros(RIR_len, M, numnoise);
for p=1:numnoise
    
    [htempn, Sr]  = ISM(pm,v_pos(:,p)',room_dim,rev_time,Nrefl,RIR_len, Rd, [],Tw,Fc,fs_RIR,cair);
    RIR_noise(:,:,p) = htempn; 
    
end
noiseIR = RIR_noise;

%DISCLAIMER: THE LENGTHS OF THE VECTORS IN audiosignal and noisesignal
%should be the same. Can pad with zeros. But do this manually.

%Check to see how many microphone signals, noise sources and audio sources
%exist. Audio sources are the sources as defined from the room GUI, similar
%as noise sources. These are independent of the audio and noise signals.
%The sources just refer to the excitation points. The signals refer to the
%actual signals that will be used for excitation.

%---LOADING SIGNALS---%
Ntotal = siglen*fs_RIR;

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

n_local = n_local_init;



