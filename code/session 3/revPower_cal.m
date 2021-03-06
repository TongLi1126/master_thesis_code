function rho_r_esti = revPower_cal(Rnn,Ryy,Rxx_real,gama,k)
  % model.Rx choose which model to estimate Rx, 1: real update 
        % 2: direct estimate 3: decesion directed 
           if (abs(Ryy(1,1))> abs(Rnn(1,1)))
           Rxx = Ryy;
           else
           Rxx= 1e-6*Rnn;
           end
           [vector,eigen] = eig(Rxx);
           eigen(find(eigen<0)) = 0;
           Rxx_remove0 = pinv(vector')* eigen* pinv(vector);
           L = chol(gama);
           [~,lamada] = eig(inv(L) *  Rxx_remove0 * inv(L'),'vector');
          
           if k<42
               rho_r_esti = abs(lamada(2)/10) ;
           else
               rho_r_esti = abs(lamada(2));
           end
%        
end

