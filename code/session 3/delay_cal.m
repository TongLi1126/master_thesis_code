function [K,a,b] = delay_cal(Rnn,Ryy,K,alpha,N_freqs)
% delay_cal 
% calculate delay and decay factor K adaptively, for mic = 2 K1=1 K2 = a2*exp(jw*b2)
% alpha: learning rate, K: delay vector [1; a2*exp(jw*b2)]
% Rss: estimated power spectrum
w = 0: pi/(N_freqs-1) :pi;
a(1) = 1; b(1) = 0;
beta = 1.1;
PS1 = ones(N_freqs,1);
for m = 2:size(K,2)
     a(m) = abs(K(1,m));
     b(m) =  angle(K(2,m))/(pi/(N_freqs-1));
     gredient_a(m) = 0;
     gredient_b(m) = 0;
    for k = 1:N_freqs
       
        if(Ryy{k}(1,1)>beta*Rnn{k}(1,1))
            PS1(k) = abs(Ryy{k}(1,1) -Rnn{k}(1,1)); 
        else
            PS1(k) = (beta-1)*abs(Rnn{k}(1,1)); 
        end
        E = Ryy{k}-Rnn{k}-PS1(k)*(K(k,:)'*K(k,:));
        theta = pi*(k-1)/(N_freqs-1);
        gredient_a(m) = gredient_a(m) + (-4*(PS1(k)* real(exp(1i*theta* b(m))*...
                        (K(k,:)* E(:,m)))));
        gredient_b(m) = gredient_b(m) + (-2 * a(m)*theta*PS1(k)* ...
                        imag(exp(1i*theta* b(m)) * (K(k,:)* E(:,m))));
    
    end
    
    a_update(m) = a(m) - alpha *  gredient_a(m);
    b_update(m) = b(m) - alpha *  gredient_b(m);
    K(:,m) = a_update(m) .* exp(1i.* w *b_update(m));
end
    K(:,1) = 1;
end
