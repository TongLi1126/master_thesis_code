%11/02/2016
%[SNR_weighted] = SI_SNRatio(cleansp, filtsp, start_time, fs, vad)
%
%The following is a measure of an intelligibility weighted Signal to Noise
%ratio (in dB). The formula is given by eqn (13) of Spriet et
%al "Robustness Analysis of MWF and GSC for multimicrophone NR in hearing
%aid applications" (IEEE Trans. Speech and Audio Proc, 2005).
%
%Input Arguments:
%filtsp  - filtered speech signal from output of algorithm
%
%filtn  - filtered noise signal from output of algorithm
%
%start_time - time at which you would like to start the comparison of the
%algorithm (included in case we only want to calculate the results after
%some settling time of the algorithm).
%
%vad - ideal vad reference consisting of 1's and 0's indicating where there
%is speech and noise. We need this as we evaluate this algorithm for speech
%only parts.
%
%fs - sampling frequency (Hz)
%
%
%Output arguments:
%SNR_weighted - the weighted SNR according to the importance bands in dB

function [SNR_weighted] = SI_SNRatio(filtsp, filtn, start_time, fs, vad)

stN = round(start_time*fs)+1;       % starting sample for SNR calc
speechsig = filtsp(stN:end,1);      % truncate the signal accordingly
noisesig = filtn(stN:end,1);        % truncate the signal accordingly
VADtrun = vad(stN:end);             % truncated VAD

%truncation according to VAD
S = speechsig(VADtrun==1,1);       
Nn = noisesig(VADtrun==0,1);

%Band importance function - Taken from table 3 of ANSI S3.5 1997. 1/3
%octave bands are used
r=2^(1/6);
I=[83 95 150 289 440 578 653 711 818 844 882 898 868 844 771 527 364 185]*1e-4;         %Importance weightings @ 1/3 octaves
F=[160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000];  %1/3 octave centre freqs
f2=F*r;
n=sum(fs/2>f2); %ensuring that our max freq of interest is not above fs/2
F=F(1:n);
I=I(1:n);
I=I/sum(I);     %normalization in case of the frequency bands of interest changing
%disp([num2str(n) '/18 banden, max: ' num2str(sum(I)*100) '%']);

SNR=zeros(1,n);

%calculation of the SNR values in the different bands
for i=1:n
%    [b, a] = third_octave_filter(F(i), fs);  %Gabrielson's 1/3 band octave filter
   [b,a]	 = oct3dsgn(F(i),fs,3);         %Use matlab's 1/3 band octave filter if don't have the one above.
   
   oct_out_sp = filter(b,a,S);
   Sp(i) = var(oct_out_sp);                 %Filtered speech power
   
   %noise spectrum level -> should include the noncorrelated as well as the 
   %correlated noise (i.e. reverberation). Here, only the uncorrelated noise is
   %taken into account (-> only valid in case of small reverberation)!! 
   
   oct_out_n = filter(b,a,Nn);
   Np(i) = var(oct_out_n);                   %Filtered noise power
   
   SNR(i)=10*log10(Sp(i)/Np(i));

% condition to provide upper limit for SNR
%    if clipping==1 
%       SNR(i)=min(max(0,SNR(i)),30);
%    end

end

SNR_weighted=I*SNR';    %weighted SNR
end

function [B,A] = oct3dsgn(Fc,Fs,N)
% OCT3DSGN  Design of a one-third-octave filter.
%    [B,A] = OCT3DSGN(Fc,Fs,N) designs a digital 1/3-octave filter with 
%    center frequency Fc for sampling frequency Fs. 
%    The filter is designed according to the Order-N specification 
%    of the ANSI S1.1-1986 standard. Default value for N is 3. 
%    Warning: for meaningful design results, center frequency used
%    should preferably be in range Fs/200 < Fc < Fs/5.
%    Usage of the filter: Y = FILTER(B,A,X). 
%
%    Requires the Signal Processing Toolbox. 
%
%    See also OCT3SPEC, OCTDSGN, OCTSPEC.

% Author: Christophe Couvreur, Faculte Polytechnique de Mons (Belgium)
%         couvreur@thor.fpms.ac.be
% Last modification: Aug. 25, 1997, 2:00pm.

% References: 
%    [1] ANSI S1.1-1986 (ASA 65-1986): Specifications for
%        Octave-Band and Fractional-Octave-Band Analog and
%        Digital Filters, 1993.

% Copyright (c) 1997, Christophe COUVREUR
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

if (nargin > 3) | (nargin < 2)
  error('Invalide number of arguments.');
end
if (nargin == 2)
  N = 3; 
end
if (Fc > 0.88*(Fs/2))
  error('Design not possible. Check frequencies.');
end
  
% Design Butterworth 2Nth-order one-third-octave filter 
% Note: BUTTER is based on a bilinear transformation, as suggested in [1]. 
pi = 3.14159265358979;
f1 = Fc/(2^(1/6)); 
f2 = Fc*(2^(1/6)); 
Qr = Fc/(f2-f1); 
Qd = (pi/2/N)/(sin(pi/2/N))*Qr;
alpha = (1 + sqrt(1+4*Qd^2))/2/Qd; 
W1 = Fc/(Fs/2)/alpha; 
W2 = Fc/(Fs/2)*alpha;
[B,A] = butter(N,[W1,W2]); 

end