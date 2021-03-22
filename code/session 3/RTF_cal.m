function [d,gama] = RTF_cal(fs_RIR,m_pos,s_pos,num_mics,N_freqs)
% RTF_cal calculate RTF d(steering vector) and diffuse reverberation gama
% d=[1 e^(-jw*t2)], can be calculate by DOA or RIR_sources
% gama is calculated on each k
d = ones(num_mics,N_freqs); 
gama = cell(N_freqs,1);
gama_update = ones(num_mics,num_mics);
% w = 0 : pi / (nfft-1) : pi; % digital w (0-pi)
dis_mic = m_pos(2,2)-m_pos(1,2); % distance between mics
DA = mean(m_pos) -s_pos; 
DB = m_pos(1,:) - m_pos(2,:);
DOA = acos(dot(DA,DB)/sqrt(norm(DA)*norm(DB)));
for k = 1:N_freqs
    for p =  1:num_mics
        for q = 1:num_mics
            gama_update(p,q) = sinc(((2*pi*(k-1)/(N_freqs-1))* abs((p-q))...
                *dis_mic*fs_RIR/340));
        end
         d(p,k) = exp(-j * ((2*pi*(k-1)/(N_freqs-1)) * (p-1)* dis_mic*cos(DOA)*fs_RIR/340));
    end
    gama{k} = gama_update;
end
    
end

