% Lab 2 for Digital Audio Signal Processing Lab Sessions
% Session 2: Binaural synthesis and 3D audio: OLA and WOLA frameworks
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2020
%
% The following is the skeleton code for the analysis stage of the WOLA method, which you need to
% complete


function [X,f] = WOLA_analysis(x,fs,window,nfft,noverlap,varargin)
%WOLA_analysis  short-time fourier transform
% INPUT:
%   x           : input time signal(s) (samples x channels)
%   fs          : sampling rate
%   window      : window function
%   nfft        : FFT size
%   noverlap    : frame overlap; default: 2 (50%)
%   varargin    : (optional) transfer function      Lh X 1
%
% OUTPUT:
%   X           : STFT matrix (bins x frames x channels)
%   f           : frequency vector for bins


% use only half FFT spectrum
N_half = nfft / 2 + 1;

% get frequency vector
f = 0:(fs / 2) / (N_half - 1):fs / 2;

% init
L = floor((length(x) - nfft + (nfft / noverlap)) / (nfft / noverlap));
M = size(x,2);
X = zeros(N_half, L, M);
if(isa(window,'dsp.Window'))
    w = window(ones(nfft,1));
elseif(length(window)==nfft)
    w = window;
end

% for possible filtering
H_half = ones(N_half,1);
if(~isempty(varargin))
    h = varargin{1};
    H = fft(h,nfft);
    H_half = H(1:N_half);
end
for m = 1:M
    for index = 0:L-1 % Frame index
        left = index*nfft/noverlap+1; %%%
        right = left+nfft-1;
        x_window = w.*x(left:right,m); % apply window
        x_fft = fft(x_window);
%        periodogram(x_window,rectwin(length(x_fft)),length(x_fft),fs)
%         s_multi = pmtm(x(left:right,m),3,nfft,fs);
%          [pss,s_multi] = pmtm(x(left:right,m),3,nfft,fs);
%          pss2 = pmtm(x(left:right,m),3,nfft,fs);
%         figure
%         subplot(2,1,1) 
% %         pmtm(x(left:right,m),3,nfft,fs)
%         plot(f,10*log10(s_multi)+14)
%          subplot(2,1,2) 
%          plot(f,10*log10(2*(abs(x_fft(1:257))/nfft).^2))
        X(:,index+1,m) = x_fft(1:N_half).*H_half(:);
        X_multi(:,index+1,m) = 10*ones(N_half,1);
        P_multi(:,index+1,m) = 10*ones(N_half,1);
        
    end
end

end
