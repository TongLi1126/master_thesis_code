% Lab 2 for Digital Audio Signal Processing Lab Sessions
% Session 2: Binaural synthesis and 3D audio: OLA and WOLA frameworks
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2020
%
% The following is the skeleton code for the synthesis stage of the WOLA method, which you need to
% complete


function x_syn = WOLA_synthesis(X_half,window,nfft,noverlap)
%WOLA_synthesis inverse short-time fourier transform.
%
% INPUT:
%   X           : input matrix (bins x frames x channels)
%   window      : window function
%   nfft        : FFT size
%   noverlap    : frame overlap; default: 2 (50%)
%
% OUTPUT:
%   x           : output time signal(s)


L = size(X_half,2);
M = size(X_half,3);

% ## Perform IFFT

X = [X_half ; conj(flip(X_half(2:nfft/2,:,:)))];
x = ifft(X);


% ## Apply synthesis window
if(isa(window,'dsp.Window'))
    w = window(ones(nfft,1));
elseif(length(window)==nfft)
    w = window;
end
v = synthesis_window(w,noverlap); % calculate inverse window
x = x.*v;
x = permute(x,[1 3 2]);

% ## Obtain re-synthesised signals
x_syn = zeros(L*nfft/noverlap+nfft*(1-1/noverlap),M);
for index = 0:L-1
    left = index*nfft/noverlap+1;
    right = left+nfft-1;
    x_syn(left:right,:) = x_syn(left:right,:) + x(:,:,index+1);
end

    

end

function v =  synthesis_window(window,noverlap)
   %WOLA_synthesis inverse window.
%
% INPUT:
%   window      : window function: nfft x 1
%   nfft        : FFT size
%   noverlap    : frame overlap; default: 2 (50%)
%
% OUTPUT:
%   v           : inversed window
nfft = size(window,1);
downsample = nfft/noverlap;
% check if nfft / overlap is integer
assert(downsample == floor(downsample))
v = zeros(nfft,1);
for i = 1: downsample
    v(i:downsample:end) = window(i:downsample:end)./(sum(window(i:downsample:end).^2));  % DSP charpter14 page 14/15
end
end
