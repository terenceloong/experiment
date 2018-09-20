% Pre-process.

%% 1.Correct imu
imu = imu0;

load('./calib/bias_gyro.mat');
% load('./calib/bias_acc.mat');
% load('./calib/c_mag.mat');
% load('./calib/S_mag.mat');
c_mag = [4.6; 3.0; -6.5];

for k=1:size(imu0,1)
    imu(k,2:4) = imu0(k,2:4) - bias_gyro + [0.5*0,0.5*0,0.5*0];
%     imu(k,5:7) = imu0(k,5:7) - bias_acc;
%     imu(k,8:10) = (S_mag*(imu0(k,8:10)'-c_mag))';
    imu(k,8:10) = imu0(k,8:10)-c_mag';
end

%% 2.Delete redundant imu data
k = 1;
while 1
    if imu(k,1)>info(1,2)
        imu(1:k-1,:) = [];
        break;
    end
    k = k+1;
end

k = 0;
while 1
    if imu(end-k,1)<info(end,2)
        imu(end-k+1:end,:) = [];
        break;
    end
    k = k+1;
end

%% 3.Convert velocity
% a = 6378137;
% f = 1/298.257223563;
% vel = zeros(size(pos,1),3); %[vn, ve, vd]
% for k=1:size(pos,1)
%     lat = pos(k,2)/180*pi;
%     h = pos(k,4);
%     Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5;
%     Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5;
%     vel(k,1) = pos(k,8)*(Rm+h);
%     vel(k,2) = pos(k,9)*(Rn+h)/sec(lat);
%     vel(k,3) = -pos(k,10);
% end

%% 4.Clear variables
clearvars -except ephemeris* ion info pos vel rou doppler ele azi imu imu0
disp('Finish!');