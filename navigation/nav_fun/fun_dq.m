function dq = fun_dq(q, w)

dq = 0.5*[ 0,   -w(1), -w(2), -w(3);
          w(1),   0,    w(3), -w(2);
          w(2), -w(3),   0,    w(1);
          w(3),  w(2), -w(1),   0 ]*q;

end