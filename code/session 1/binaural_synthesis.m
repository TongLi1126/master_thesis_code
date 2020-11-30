clear;
% close all

if ispc
    addpath('..\..\audio_files');
    addpath('..\..\sim_environment');
else
    addpath('../../audio_files');
    addpath('../../sim_environment');
end
load HRTF.mat;

fs = 8000;
source_length = 10;
samples = fs * source_length;

%% process speech1
speech_filename = 'speech1.wav';
[source_signals_raw{1}, source_signals_raw{2}] = audioread(speech_filename);
x  = resample(source_signals_raw{1},fs,source_signals_raw{2});
x = x(1:samples);


binaural_sig1_1 = [x x];
binaural_sig2_1 = [x 0.5*x];
binaural_sig3_1 = [x delayseq(x,3)];
binaural_sig4_1 = fftfilt(HRTF,binaural_sig1_1);

%% process speech2
speech_filename = 'speech2.wav';
[source_signals_raw{1}, source_signals_raw{2}] = audioread(speech_filename);
x  = resample(source_signals_raw{1},fs,source_signals_raw{2});
x = x(1:samples);


binaural_sig1_2 = [x x];
binaural_sig2_2 = [0.5*x x];
binaural_sig3_2 = [delayseq(x,3) x];
binaural_sig4_2 = fftfilt(HRTF(:,[2 1]),binaural_sig1_2);

%%
binaural_sig1 = binaural_sig1_2+binaural_sig1_1;
binaural_sig2 = binaural_sig2_2+binaural_sig2_1;
binaural_sig3 = binaural_sig3_2+binaural_sig3_1;
binaural_sig4 = binaural_sig4_2+binaural_sig4_1;

%% listen
soundsc(binaural_sig4_1,fs);

% method 2, 3 and 4 make these two source signal distinguishable in
% binaural sense. 
% 2: By suppress the signal's magnitude in another direction.
% 3: By simulating a DOA of the signal.（the delay is the same with
% 90-degree HRTF)
% 4: By a left 90-degree HRTF (HRTF is measured in an anechoic room)

%%
figure
subplot(2,1,1);
plot(HRTF(:,1));
subplot(2,1,2);
plot(HRTF(:,2));

