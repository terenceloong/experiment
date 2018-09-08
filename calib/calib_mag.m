clear;
% clc;

%% 1.Calibrate bias
folder = 'data1';
D0 = fun_parse_imu(['./',folder,'/D0.DAT']);
D0 = D0(:,8:10);
n = size(D0,1);

% figure
% plot3(D0(:,1), D0(:,2), D0(:,3))
% axis equal

[c, S, lamda] = mag_calib(D0(:,1), D0(:,2), D0(:,3));
disp('c:');
disp(c');
disp('lamda:');
disp(lamda');

length = zeros(n,1);
for k=1:n
    D0(k,:) = (S*(D0(k,:)'-c))';
    length(k) = norm(D0(k,:));
end

figure
[x,y,z] = sphere;
h = surf(x,y,z, 'FaceAlpha',0.5);
h.EdgeColor = 'none';
axis equal
hold on
plot3(D0(:,1), D0(:,2), D0(:,3), 'LineWidth',2)
quiver3(0,0,0, 1,0,0, 'Color','b', 'LineWidth',3)
quiver3(0,0,0, 0,1,0, 'Color','m', 'LineWidth',3)
quiver3(0,0,0, 0,0,1, 'Color','r', 'LineWidth',3)

figure
plot(length)

%% 2.Calibrate axis
D1 = fun_parse_imu(['./',folder,'/D1.DAT']);
D2 = fun_parse_imu(['./',folder,'/D2.DAT']);
D3 = fun_parse_imu(['./',folder,'/D3.DAT']);
D1 = D1(:,8:10);
D2 = D2(:,8:10);
D3 = D3(:,8:10);
n1 = size(D1,1);
n2 = size(D2,1);
n3 = size(D3,1);

for k=1:n1
    D1(k,:) = (S*(D1(k,:)'-c))';
%     D1(k,:) = D1(k,:)/norm(D1(k,:));
end
for k=1:n2
    D2(k,:) = -(S*(D2(k,:)'-c))';
%     D2(k,:) = D2(k,:)/norm(D2(k,:));
end
for k=1:n3
    D3(k,:) = (S*(D3(k,:)'-c))';
%     D3(k,:) = D3(k,:)/norm(D3(k,:));
end

figure
[x,y,z] = sphere;
h = surf(x,y,z, 'FaceAlpha',0.5);
h.EdgeColor = 'none';
axis equal
hold on
quiver3(0,0,0, 1,0,0, 'Color','b', 'LineWidth',3)
quiver3(0,0,0, 0,1,0, 'Color','m', 'LineWidth',3)
quiver3(0,0,0, 0,0,1, 'Color','r', 'LineWidth',3)
plot3(D1(:,1), D1(:,2), D1(:,3), '*', 'Color','b')
plot3(D2(:,1), D2(:,2), D2(:,3), '*', 'Color','m')
plot3(D3(:,1), D3(:,2), D3(:,3), '*', 'Color','r')

%-------------------------------------------------------------------------%
meanx = mean(D1(:,1));
meany = mean(D1(:,2));
meanz = mean(D1(:,3));
F11 = mean((meanx-D1(:,1)).^2);
F12 = mean((meanx-D1(:,1)).*(meany-D1(:,2)));
F13 = mean((meanx-D1(:,1)).*(meanz-D1(:,3)));
F21 = mean((meany-D1(:,2)).*(meanx-D1(:,1)));
F22 = mean((meany-D1(:,2)).^2);
F23 = mean((meany-D1(:,2)).*(meanz-D1(:,3)));
F31 = mean((meanz-D1(:,3)).*(meanx-D1(:,1)));
F32 = mean((meanz-D1(:,3)).*(meany-D1(:,2)));
F33 = mean((meanz-D1(:,3)).^2);

fun = @(x) [F11*x(1) + F12*x(2) + F13*x(3) + x(1)*x(4);...
            F21*x(1) + F22*x(2) + F23*x(3) + x(2)*x(4);...
            F31*x(1) + F32*x(2) + F33*x(3) + x(3)*x(4);...
            x(1)^2 + x(2)^2 + x(3)^2 - 1];

u1 = fsolve(fun, [1;0;0;0]);
quiver3(0,0,0, u1(1),u1(2),u1(3), 'Color','b', 'LineWidth',2)
%---------------------------------
meanx = mean(D2(:,1));
meany = mean(D2(:,2));
meanz = mean(D2(:,3));
F11 = mean((meanx-D2(:,1)).^2);
F12 = mean((meanx-D2(:,1)).*(meany-D2(:,2)));
F13 = mean((meanx-D2(:,1)).*(meanz-D2(:,3)));
F21 = mean((meany-D2(:,2)).*(meanx-D2(:,1)));
F22 = mean((meany-D2(:,2)).^2);
F23 = mean((meany-D2(:,2)).*(meanz-D2(:,3)));
F31 = mean((meanz-D2(:,3)).*(meanx-D2(:,1)));
F32 = mean((meanz-D2(:,3)).*(meany-D2(:,2)));
F33 = mean((meanz-D2(:,3)).^2);

fun = @(x) [F11*x(1) + F12*x(2) + F13*x(3) + x(1)*x(4);...
            F21*x(1) + F22*x(2) + F23*x(3) + x(2)*x(4);...
            F31*x(1) + F32*x(2) + F33*x(3) + x(3)*x(4);...
            x(1)^2 + x(2)^2 + x(3)^2 - 1];

u2 = fsolve(fun, [0;1;0;0]);
quiver3(0,0,0, u2(1),u2(2),u2(3), 'Color','m', 'LineWidth',2)
%---------------------------------
meanx = mean(D3(:,1));
meany = mean(D3(:,2));
meanz = mean(D3(:,3));
F11 = mean((meanx-D3(:,1)).^2);
F12 = mean((meanx-D3(:,1)).*(meany-D3(:,2)));
F13 = mean((meanx-D3(:,1)).*(meanz-D3(:,3)));
F21 = mean((meany-D3(:,2)).*(meanx-D3(:,1)));
F22 = mean((meany-D3(:,2)).^2);
F23 = mean((meany-D3(:,2)).*(meanz-D3(:,3)));
F31 = mean((meanz-D3(:,3)).*(meanx-D3(:,1)));
F32 = mean((meanz-D3(:,3)).*(meany-D3(:,2)));
F33 = mean((meanz-D3(:,3)).^2);

fun = @(x) [F11*x(1) + F12*x(2) + F13*x(3) + x(1)*x(4);...
            F21*x(1) + F22*x(2) + F23*x(3) + x(2)*x(4);...
            F31*x(1) + F32*x(2) + F33*x(3) + x(3)*x(4);...
            x(1)^2 + x(2)^2 + x(3)^2 - 1];

u3 = fsolve(fun, [0;0;1;0]);
quiver3(0,0,0, u3(1),u3(2),u3(3), 'Color','r', 'LineWidth',2)

U = [u1(1:3),u2(1:3),u3(1:3)];
S = U'*S;

% fun = @(x) [- U(1,1)*cos(x(2))*sin(x(1)) + U(2,1)*cos(x(2))*cos(x(1))...
%             + U(1,2)*(-sin(x(2))*sin(x(3))*sin(x(1))-cos(x(1))*cos(x(3))) + U(2,2)*(sin(x(2))*sin(x(3))*cos(x(1))-sin(x(1))*cos(x(3)))...
%             + U(1,3)*(-sin(x(2))*cos(x(3))*sin(x(1))+cos(x(1))*sin(x(3))) + U(2,3)*(sin(x(2))*cos(x(3))*cos(x(1))+sin(x(1))*sin(x(3)));...
%             - U(1,1)*sin(x(2))*cos(x(1)) - U(2,1)*sin(x(2))*sin(x(1)) - U(3,1)*cos(x(2))...
%             + U(1,2)*cos(x(2))*sin(x(3))*cos(x(1)) + U(2,2)*cos(x(2))*sin(x(3))*sin(x(1)) - U(3,2)*sin(x(2))*sin(x(3))...
%             + U(1,3)*cos(x(2))*cos(x(3))*cos(x(1)) + U(2,3)*cos(x(2))*cos(x(3))*sin(x(1)) - U(3,3)*sin(x(2))*cos(x(3));...
%               U(1,2)*(sin(x(2))*cos(x(3))*cos(x(1))+sin(x(1))*sin(x(3))) + U(2,2)*(sin(x(2))*cos(x(3))*sin(x(1))-cos(x(1))*sin(x(3))) + U(3,2)*cos(x(2))*cos(x(3))...
%             + U(1,3)*(-sin(x(2))*sin(x(3))*cos(x(1))+sin(x(1))*cos(x(3))) + U(2,3)*(-sin(x(2))*sin(x(3))*sin(x(1))-cos(x(1))*cos(x(3))) - U(3,3)*cos(x(2))*sin(x(3))];
% angle = fsolve(fun, [0;0;0]);
% T = angle2dcm(angle(1), angle(2), angle(3));
% 
% quiver3(0,0,0, T(1,1),T(1,2),T(1,3), 'Color','b', 'LineWidth',1)
% quiver3(0,0,0, T(2,1),T(2,2),T(2,3), 'Color','m', 'LineWidth',1)
% quiver3(0,0,0, T(3,1),T(3,2),T(3,3), 'Color','r', 'LineWidth',1)
% 
% S = T*S;
%-------------------------------------------------------------------------%

D1 = fun_parse_imu(['./',folder,'/D1.DAT']);
D2 = fun_parse_imu(['./',folder,'/D2.DAT']);
D3 = fun_parse_imu(['./',folder,'/D3.DAT']);
D1 = D1(:,8:10);
D2 = D2(:,8:10);
D3 = D3(:,8:10);
n1 = size(D1,1);
n2 = size(D2,1);
n3 = size(D3,1);

for k=1:n1
    D1(k,:) = (S*(D1(k,:)'-c))';
%     D1(k,:) = D1(k,:)/norm(D1(k,:));
end
for k=1:n2
    D2(k,:) = -(S*(D2(k,:)'-c))';
%     D2(k,:) = D2(k,:)/norm(D2(k,:));
end
for k=1:n3
    D3(k,:) = (S*(D3(k,:)'-c))';
%     D3(k,:) = D3(k,:)/norm(D3(k,:));
end

figure
[x,y,z] = sphere;
h = surf(x,y,z, 'FaceAlpha',0.5);
h.EdgeColor = 'none';
axis equal
hold on
quiver3(0,0,0, 1,0,0, 'Color','b', 'LineWidth',3)
quiver3(0,0,0, 0,1,0, 'Color','m', 'LineWidth',3)
quiver3(0,0,0, 0,0,1, 'Color','r', 'LineWidth',3)
plot3(D1(:,1), D1(:,2), D1(:,3), '*', 'Color','b')
plot3(D2(:,1), D2(:,2), D2(:,3), '*', 'Color','m')
plot3(D3(:,1), D3(:,2), D3(:,3), '*', 'Color','r')

figure
plot(D1(:,1))
hold on
grid on
plot(D1(:,2))
plot(D1(:,3))

figure
plot(D2(:,1))
hold on
grid on
plot(D2(:,2))
plot(D2(:,3))

figure
plot(D3(:,1))
hold on
grid on
plot(D3(:,2))
plot(D3(:,3))

S_mag = S;
c_mag = c;
save S_mag.mat S_mag
save c_mag.mat c_mag