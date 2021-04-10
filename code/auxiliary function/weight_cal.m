function [W_update] = weight_cal(Ryy, Rnn,Rxx,rho_rev,rho_ss,TX_Mask,gama,h,func_model,mu)
    I = eye(2); 
    
 
% signal general eigenvalue decomposition
    [sig_yy,sig_nn,sig_ss,d,Q,Qh,V] = GEVD(Rnn,Ryy);
    Rss = sig_ss(1,1)* Q(:,1)*Qh(1,:);

    mvdr_n1 = 1; 
    mvdr_n2 = 1;
   switch func_model 
       case 1 % basic MWF general eigenvalue decomposition
            % See DASP-CH4-37/40
            % d(1) corresponds to (d_noise/d_s&n)
            W_update = (1-mu/(d(1)+mu-1))*Qh(1,1).*V(:,1);% Final expression for filter.
            W_compare = W_update;
      
       case 2   % basic MWF in direct way
            alpha = 0; m = 1;
            denominator = (Q*(sig_yy+(mu-1)*sig_nn)*Qh) + alpha * I;
            nominator =  Rss + alpha * m * I;
            W_update = pinv(denominator)* nominator *[1;0];
%               W_update = pinv( Q*(sig_yy+(mu-1)*sig_nn)*Qh)* sig_ss *[1;0];
       case 3    % MWF with dereverberation
%            Rxx =   Ryy - Rnn; % estimate speech correlation matrix, consist of target sppech + late reverberate
%            [vector,eigen] = eig(Rxx);
%            eigen(find(eigen<0)) = 0;
%            Rxx = pinv(vector')* eigen* pinv(vector);
%            L = chol(gama);
%            [~,lamada] = eig(inv(L) * Rxx * inv(L'),'vector');
%            rho_rev = lamada(2);
%            Rrr =  rho_rev * gama; 
%            Rnn_rev = Rnn + Rrr;
%            [sig_yy,sig_nn,sig_ss,d,Q,Qh,V] = GEVD(Rnn_rev,Ryy);
%            rho_ss = sig_ss(1,1);
%            h = Q(:,1);
%            W_single = rho_ss/(rho_ss + mu /(h'*inv(Rrr + Rnn)*h));
%            W_mvdr = Q(1,1)*(1/(h'*inv(Rnn+Rrr)*h))*inv(Rnn+Rrr)*h;
          % traditional d bad performance
           Rrr = rho_rev *gama;           
           W_single = rho_ss/(rho_ss + mu /(abs(h'*pinv(Rrr + Rnn)*h)));
           W_single = max(0,W_single );
           W_mvdr = (1/abs(h'*pinv(Rnn+Rrr)*h))*pinv(Rnn+Rrr)*h;
           W_update =  W_single * W_mvdr; 
           W_compare = (1-mu/(d(1)+mu-1))*Qh(1,1).*V(:,1);
       case 4 % MWF Hjd with MVDR beamformer 
           xi = 0.005;          
           psd_TX = (10^((TX_Mask-90)/10))* 512*512/2;
           % USE Q->d, d is ATF
%            h = Q(:,1);                    
           mvdr_n1 = sig_nn(1,1)*1; 
           mvdr_n2 = (1/(h'*inv( Rnn)*h))*1;
%            W_single = (xi + sqrt(psd_TX/(15*mvdr_n1*(Q(1,1)^2)))); % 
%            W_mvdr = Qh(1,1)*(1/(h'*inv( Rnn)*h))*pinv(Rnn)*h;
           % USE traditional d, d is RTF
           Rrr = rho_rev *gama;
        
             W_single = (xi + sqrt(psd_TX*1/10*abs(h'*pinv(Rnn+Rrr)*h))); %  traditional d
             W_mvdr = (1/abs(h'*pinv(Rnn+Rrr)*h))*pinv(Rnn+Rrr)*h;
           W_update =  W_single * W_mvdr; 
           W_compare = Qh(1,1)*(1/(h'*pinv( Rnn)*h))*pinv(Rnn)*h;
       case 5 % MWF with MVDR beamformer
           % d is ATF
           h = Q(:,1);
           rho_ss = sig_ss(1,1);
%          h_tude = h./h(1,1); rho_ss_tude = rho_ss/(h(1,1)^2);
           W_single = rho_ss/(rho_ss + mu /(h'*inv( Rnn)*h));
           W_mvdr = Qh(1,1)*(1/(h'*inv( Rnn)*h))*pinv(Rnn)*h;
           W_update =  W_single * W_mvdr; 
           W_compare = (1-mu/(d(1)+mu-1))*Qh(1,1).*V(:,1);% Final expression for filter.
           
   end
   
end
function [sig_yy,sig_nn,sig_ss,d,Q,Qh,V] = GEVD(Rnn,Ryy)
    % Computing the MWF filter using a generalised eigenvalue
    % decomposition of the correlation matrices.

    [V,d] = eig(Ryy,Rnn,'vector');
    [d,ind] = sort(d,'descend');
    V = V(:,ind);
    Qh = pinv(V);
    sig_yy = V'*Ryy *V;
    sig_yy = [sig_yy(1,1) 0;0 sig_yy(2,2)];
    sig_nn = V'*Rnn*V;
    sig_nn = [sig_nn(1,1) 0;0 sig_nn(2,2)];
    Q= Qh';
    sig_ss = sig_yy(1,1) - sig_nn(1,1) ; 
    sig_ss = [sig_ss 0; 0 0];
end



