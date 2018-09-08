function [X, P, E, Xc] = kalman_filter(Phi, X, P, Q, H, Z, R)

switch nargin
    case 4 %only time update
        X = Phi*X;
        P = Phi*P*Phi' + Q;
    case 7 %time update and measure update
        X = Phi*X;
        P = Phi*P*Phi' + Q;
        K = P*H' / (H*P*H'+R);
        E = Z - H*X;
        Xc = K*E;
        X = X + Xc;
        P = (eye(length(X))-K*H)*P;
end
P = (P+P')/2;

end