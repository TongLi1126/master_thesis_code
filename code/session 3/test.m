ryy = Ryy{250}; rnn = Rnn{250};
[vv,dd] = eig(ryy,rnn,'vector');
sig_yy = vv'*ryy *vv;
sig_nn = vv'*rnn*vv;
