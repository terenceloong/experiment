%% 1.Magnetometer vector length
mag_norm = zeros(size(imu,1),1);
for k=1:size(imu,1)
    mag_norm(k) = norm(imu(k,8:10));
end
figure
plot(mag_norm)

%% 2.Heading
figure
plot(imu(1:2:end,1),nav(:,7))
hold on
plot(imu(3:2:end,1),heading(:,1))