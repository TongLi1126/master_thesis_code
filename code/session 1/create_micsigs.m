clear;
% close all

if ispc
    addpath('..\..\audio_files');
    addpath('..\..\sim_environment');
else
    addpath('../../audio_files');
    addpath('../../sim_environment');
end


speech_filename{1} = 'speech1.wav';
%speech_filename{2} = 'speech2.wav';

noise_filename{1} = 'Babble_noise1.wav';
source_length = 10;

load Computed_RIRs.mat;

%% load audio file
for i=1:numel(speech_filename)
    [source_signals_raw{i, 1}, source_signals_raw{i, 2}] = audioread(speech_filename{i});
end

for i=1:numel(noise_filename)
    [noise_signals_raw{i, 1}, noise_signals_raw{i, 2}] = audioread(noise_filename{i});
end

%% resample them using fs_RIR, then truncate to source_length
samples = fs_RIR * source_length;
for i=1:numel(speech_filename)
    source_signals{i}  = resample(source_signals_raw{i,1},fs_RIR,source_signals_raw{i,2});
    source_signals{i} = source_signals{i}(1:samples);
end

for i=1:numel(noise_filename)
    noise_signals{i} =  resample(noise_signals_raw{i,1},fs_RIR,noise_signals_raw{i,2});
    noise_signals{i} = noise_signals{i}(1:samples);
end
source_signals = cell2mat(source_signals);
noise_signals = cell2mat(noise_signals);

%%
if numel(speech_filename) > 0
    RIR_sources_perm = permute(RIR_sources, [1, 3, 2]);
end
if numel(noise_filename) > 0
    RIR_noises_perm = permute(RIR_noise, [1, 3, 2]);
end

if(isempty(RIR_noises_perm))
    IRs = RIR_sources_perm;
    signals = source_signals;
else
    IRs = [RIR_sources_perm, RIR_noises_perm];
    signals = [source_signals noise_signals];
end
%% generate mic signals, save with fs_RIR
mics_num = size(RIR_sources,2);
mic = zeros(samples,mics_num);
for i=1:mics_num
    mic(:,i) = sum(fftfilt(squeeze(IRs(:,:,i)),signals),2);
end
save('mic','mic','fs_RIR');
%% plot
figure;
plot(mic(:, 1));
hold on;
plot(mic(:, 2));
%plot(micSignals(:, 3));
legend({'Mic 1', 'Mic 2'});
hold off;
figure(2);
plot(RIR_sources(:, 1),'r');
hold on;
plot(RIR_sources(:, 2),'g');
legend({'RIR_sources 1', 'RIR_sources 2'});
hold off;
%% play

soundsc(mic(:, 1), fs_RIR);
% Q2 Low reviberation better
% Q3 Low reviberation better
