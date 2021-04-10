function [T_g,bark] = mask_cal(signal,f,fs,nfft,t,index_fig)
%% change frequency domain to bark domain
bark = hz2bark(f);
N_freqs = nfft/2 +1;


%% calculate ATH, auditory threshold
Tq = 3.64*(f/1000).^-0.8 -6.5*exp(-0.6*(f/1000-3.3).^2)+0.001*(f/1000).^4;
Tq(2) = 40;Tq(3) = 30;Tq(4) = 25;
%% spectral analysis
spl = psd2spl(signal,nfft);
% determine tonal masker and power of tonal masker
%% determine neighbour delta_k of a tonal masker
 ST = []; P_TM=[];
for k = 2:N_freqs-1
    if(k< 64)
    delta_k= 2;
    elseif(k>=64 && k<128)
    delta_k = [2,3];
    
    else
    delta_k = [2,3,4,5,6];
    end
  % determine lower bound and upper bound of neighboorhood of tonal masker
  down = (k- delta_k>1).*(k- delta_k)+(k- delta_k<=1).*1 ;  
  up = (k+delta_k >N_freqs).*N_freqs +(k+delta_k <=N_freqs).*(k+delta_k) ; 
 
  %% tonal masker and power of tonal masker
  if(spl(k)>max(spl(k-1),spl(k+1)) && (spl(k)-7)>max([spl(down),spl(up)],[],'all') )
   ST = [ST;spl(k) k]; % power of only tonal masker
    % power of tonal masker , additive power of three adjecent spectral lines
   P_TM = [P_TM;10*log10(10^(0.1*spl(k-1))+10^(0.1*spl(k))+10^(0.1*spl(k+1))) k];  
  end  %
  
end
%%%% noisek[]  store the fre bins in bark domain, 1st column stands for fre,2nd
%%%% column stands for bark domain
i = 1; noisek = [];
for k = 1:N_freqs

    if(bark(k)>i)
        i=i+1; 
    end
     noisek = [noisek;k i];    
end

%%%%%%%%%%%%%%%%%%%%%%%determine noise masker and power of noise masker%%
% remove spectral lines contribute to tonal maskers
if(ST)
spec_contribute2TM =[];
%%%% define neighborhood of tonal masker
    for i=1:size(ST,1) 
        if(ST(i,2)<64 && ST(i,2)>3)
            delta_k = 2;
        elseif(ST(i,2)>=64 && ST(i,2)<128)
            delta_k = 3;
        elseif(ST(i,2)>=128 && ST(i,2)<=257)
            delta_k = 6;
        else 
            delta_k = 1;
        end
    up = (ST(i,2)+delta_k <=N_freqs).*(ST(i,2)+delta_k) +(ST(i,2)+delta_k>N_freqs)*N_freqs;
    down = (ST(i,2)-delta_k>=1).*(ST(i,2)-delta_k) +(ST(i,2)-delta_k<1)*1;
%     up = ST(i,2) +1; 
%     down = ST(i,2) -1;
    spec_k_contri = down:1:up;
    spec_contribute2TM = [spec_contribute2TM;spec_k_contri'];
    end

noisek_residule = noisek;
noisek_residule (spec_contribute2TM,:) = [];
else
    noisek_residule = noisek;
end
%% calculate geomean of each bark as index of noise masker 
%% power additive of spectral lines which don't contribut tonal masker

P_NM = [];
for i=1:max(noisek(:,2)) % number of barks
    temp_index = noisek(find(noisek(:,2)==i),1); % index k SORT based on bark
    temp_npower = spl(noisek_residule(find(noisek_residule(:,2)==i),1));
    k_n = round(geomean(temp_index));
    P_NM = [P_NM; 10*log10(sum(10.^(0.1*temp_npower))) k_n]; % power of non tonal masker
    % composed of noise mask power per critical bank and index of noise
    % masker
end

%%% decimation P_TM,P_NM based on 2 criteria
% %1.ONLY select maskers above Tq_k 2.replace pair of
%%% maskers occuring withen 0.5 bark by the stronger of two
P_TM = mask_decimation(P_TM,bark,Tq);
P_NM = mask_decimation(P_NM,bark,Tq);


% calculate spread function AND indivadual threshold
% T_TM,T_NM based on spread function 
% 1.generate SF matrix (N_freq *number of maskers) S(i,j) correspond to effect
% of masker located at j makes to i(frequency bins)
% 2. generate Threshold matrix(N_freq *number of maskers)
% TM/NM(i,j),correspond to j (masker)influence i(maskee)
% Consider computation cost we only consider 10 neighboring points of
% masker
Tonal_Flag_TM = true;Tonal_Flag_NM = false;
[T_TM,SF_tm] = get_Threshold(P_TM,N_freqs,size(P_TM,1),bark,Tonal_Flag_TM);
[T_NM,SF_nm] = get_Threshold(P_NM,N_freqs,size(P_NM,1),bark,Tonal_Flag_NM);

%%%% calculate globle threshold
% T_g based on power additive of T_TM,T_NM,T_q(ATH)
T_q = Tq';% SUM(:,2) sum by row
T_g = 10.*log10(10.^(0.1*T_q)+sum(10.^(0.1*T_TM),2)+sum(10.^(0.1*T_NM),2)) ;
%%%%%

%%% visualize T_globle OR T_TM,T_NM,T_q

illustrate_mask(f,bark,spl,T_q,T_g,P_TM,P_NM,T_TM,T_NM,t,index_fig);
drawnow;
pause(0.1);

end 



function P_i = mask_decimation(P_k,bark,Tq)
if(P_k)
 % remove masker below  ATH   
P_k(P_k(:,1)<Tq(P_k(:,2))',:) = [];  
% shift 0.5 bark window to remove pair of too close maskers
masker_length = size(P_k,1);
i=1;
while(i<=masker_length-1)
    if(abs(bark(P_k(i,2))-bark(P_k(i+1,2)))<=0.5)
        if(P_k(i,1) >= P_k(i+1,1))
            P_k(i+1,:) = [];            
        else
            P_k(i,:) = [];            
        end
        masker_length = masker_length-1;
    end
    i = i+1;
end
P_i = P_k;
else
   P_i = P_k;
end
end

function [TM,SF] = get_Threshold(P_k,N_freqs,N_masker,bark,tonalflag)
SF = -1./zeros(N_freqs,N_masker);
TM = -1./zeros(N_freqs,N_masker);
if(P_k)
index = P_k(:,2); % index stores location of masker
%index(j) stores location of jth masker in frequency domain k
for j = 1:length(index) % number of masker loop of masker
    % find lower bound and upper bound of spread of masker
    down_bark= ((bark(index(j))-3)<=-0.15)*(-0.15) +...
       ((bark(index(j))-3)>-0.15)*(bark(index(j))-3);    
    [~,down] = min(abs(bark -down_bark));
    if((bark(down) -down_bark)<0)
        down = down +1;
    end
    upper_bark = ((bark(index(j))+7)>=24.9695)*(24.9695) +...
       ((bark(index(j))+7)<24.9695)*(bark(index(j))+7);  
     [~,up] = min(abs(bark -upper_bark)); % not exceed boundary
     if(bark(up)>upper_bark)
         up = up-1;
     end
    for i = down:up %% loop of maskee
        %% distance between i(maskee) and j(masker) in bark domain!!!
        delta_z = bark(i)- bark(index(j));
        if(delta_z>=-3 && delta_z <-1)
            SF(i,j) = 17*delta_z -0.4*P_k(j,1)+11;
        elseif(delta_z>=-1 && delta_z<0)
            SF(i,j) =(0.4*P_k(j,1)+6)*delta_z; 
        elseif(delta_z>=0 && delta_z<1 )
            SF(i,j) = -17*delta_z;
        else
            SF(i,j) = (0.15*P_k(j,1)-17)*delta_z-0.15*P_k(j,1);
        end
    end 
    if(tonalflag)
        TM(:,j) = P_k(j,1)-0.275*bark(index(j))+SF(:,j)-6.025; % tonal masker threshold
    else
        TM(:,j) = P_k(j,1)-0.175*bark(index(j))+SF(:,j)-2.025;
    end
end
    
end
end
%%%%% illustrate spl/bark
function [] =illustrate_mask(f,bark,spl,ATH,T_g,P_TM,P_NM,T_TM,T_NM,t,index_fig)
figure(index_fig);
subplot(2,1,1)
plot(bark,spl,':',bark,ATH,'--',bark,T_g);grid on;
title(['Periodogram of bark scale on time frame ',num2str(t)]);
xlabel('Bark(z)')
ylabel('Power/Frequency (dB/bark)'); 
xlim([0 25]);
ylim([-20 100])
if(P_TM )
text(bark(P_TM(:,2)),P_TM(:,1),'x','color','r')
end
if(P_NM)
text(bark(P_NM(:,2)),P_NM(:,1),'o','color','b')
end

subplot(2,1,2)
plot(f,spl,':',f,ATH,'--',f,T_g);grid on;
title(['Periodogram of f scale on time frame ',num2str(t)]);
xlabel('f(Hz)')
ylabel('Power/Frequency (dB/f)'); 
% xlim([0 258]);
% ylim([-20 100])
%  hold on;
%  for i=1:size(T_TM,2)
%      plot(bark,T_TM(:,i),'g')
%      hold on;
%  end
%  for i=1:size(T_NM,2)
%       plot(bark,T_NM(:,i),'y')
%       hold on;
%  end
%  hold off
if(P_TM )
text(f(P_TM(:,2)),P_TM(:,1),'x','color','r')
end
if(P_NM)
text(f(P_NM(:,2)),P_NM(:,1),'o','color','b')
end

%legend('Signal','ATH','T_g');
end


