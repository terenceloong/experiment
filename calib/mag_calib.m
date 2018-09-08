function [c, S, lamda] = mag_calib(x, y, z)

H = [x.^2, y.^2, z.^2, x, y, z];
a = (H'*H)\H'*ones(length(x),1);
Q = diag(a(1:3));
c = (-0.5*a(4:6)'/Q)';
M = Q/(1+c'*Q*c);
S = sqrt(M);
lamda = 1./diag(S);

% H = [x.^2, y.^2, z.^2, x.*y, x.*z, y.*z, x, y, z];
% a = (H'*H)\H'*ones(length(x),1);
% Q = [a(1),a(4)/2,a(5)/2; a(4)/2,a(2),a(6)/2; a(5)/2,a(6)/2,a(3)];
% c = (-0.5*[a(7),a(8),a(9)]/Q)';
% M = Q/(1+c'*Q*c);
% [V, A] = eig(M);
% [~, index] = max(abs(V));
% V(:,index) = V;
% Ad = diag(A);
% Ad(index) = Ad;
% A = diag(Ad);
% lamda = sqrt((1./diag(A)));
% S = sqrt(A)*V';

end