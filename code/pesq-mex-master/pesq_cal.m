function delta_pesq =pesq_cal(speech,noisy_sig,enhance_sig,fs) 
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
mex *.c -output ./bin/PESQ_MEX

addpath('./bin');
noisy_pesq = pesq_mex(speech, noisy_sig, fs, 'wideband');
enhance_pesq = pesq_mex(speech, enhance_sig, fs, 'wideband');
delta_pesq = enhance_pesq - noisy_pesq;
end

