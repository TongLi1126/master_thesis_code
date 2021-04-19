function rho_r_esti = revPower_cal(Rnn,Ryy,gama,k)
  % model.Rx choose which model to estimate Rx, 1: real update 
        % 2: direct estimate 3: decesion directed 
           if (abs(Ryy(1,1))> abs(Rnn(1,1)))
           Rxx = Ryy;
           else
           Rxx= 1e-6*Rnn(1,1);
           end
           [vector,eigen] = eig(Rxx);
           eigen(find(eigen<0)) = 0;
           Rxx_remove0 = pinv(vector')* eigen* pinv(vector);
           L = chol(gama);
           [~,lamada] = eig(inv(L) *  Rxx_remove0 * inv(L'),'vector');
          
           if k<81
               rho_r_esti = abs(lamada(2)/300) ;
           else
               rho_r_esti = abs(lamada(2));
           end
%        
end

