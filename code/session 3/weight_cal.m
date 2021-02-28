function [W_update,W_compare] = weight_cal(Ryy_pre, Rnn_pre,func_r,func_model,mu)
    I = eye(2); 
    alpha = 0; m = 0.9;
    switch func_r
        case 1
            Ryy = Ryy_pre;
            Rnn = Rnn_pre;
        case 2
            
            Ryy = (1-alpha)*Ryy_pre + alpha*I;
            Rnn = (1-alpha)*Rnn_pre + alpha*(1-m)*I;
    end

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
    sig_ss = (sig_yy(1,1) - sig_nn(1,1)) * Q(:,1)*Qh(1,:); 
    
    
   switch func_model 
       case 1 % basic MWF general eigenvalue decomposition
            % See DASP-CH4-37/40
            % d(1) corresponds to (d_noise/d_s&n)
            W_update = (1-mu/(d(1)+mu-1))*Qh(1,1).*V(:,1);% Final expression for filter.
            W_compare = W_update;
      
       case 2   % basic MWF in direct way
            denominator = ((1-alpha)*(Q*(sig_yy+(mu-1)*sig_nn)*Qh) + alpha * I);
            nominator = (1-alpha) * sig_ss + alpha * m * I;
            W_update = pinv(denominator)* nominator *[1;0];
            W_compare = (1-1/d(1))*Qh(1,1).*V(:,1);
%               W_update = pinv( Q*(sig_yy+(mu-1)*sig_nn)*Qh)* sig_ss *[1;0];
       case 3    % basic MWF without eigenvalue decomposition,only use diagnal information
              
              W_update = diag(pinv(diag(diag(Ryy + (mu-1)*Rnn),0))*diag(diag(Ryy-Rnn),0));
              W_compare =  W_update;
       case 5
           
           
   end
end