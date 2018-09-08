clear; clc;

folder = 'data1';

gyro = fun_parse_imu(['./',folder,'/gyro.DAT']);
accx = fun_parse_imu(['./',folder,'/accx.DAT']);
accy = fun_parse_imu(['./',folder,'/accy.DAT']);

accx = accx(:,5);
accy = accy(:,6);
accz = gyro(:,7);
gyro = gyro(:,2:4);

bias_gyro = mean(gyro,1);

bias_accx = mean(accx);
if bias_accx>0
    bias_accx = bias_accx - 1;
else
    bias_accx = bias_accx + 1;
end
bias_accy = mean(accy);
if bias_accy>0
    bias_accy = bias_accy - 1;
else
    bias_accy = bias_accy + 1;
end
bias_accz = mean(accz);
if bias_accz>0
    bias_accz = bias_accz - 1;
else
    bias_accz = bias_accz + 1;
end
bias_acc = [bias_accx, bias_accy, bias_accz];

disp('Gyro bias (deg/s):');
disp(bias_gyro);
disp('Acc bias (g):');
disp(bias_acc);

save bias_gyro.mat bias_gyro
save bias_acc.mat bias_acc