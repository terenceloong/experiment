%% Gyro bias
t_imu = (imu(3:2:end,1)-imu(1,1))/1000;

figure
plot(t_imu, imu(3:2:end,2))
hold on
plot(t_imu, bias_esti(:,3))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\epsilon\it_x\rm(\circ/s)')

figure
plot(t_imu, imu(3:2:end,3))
hold on
plot(t_imu, bias_esti(:,4))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\epsilon\it_y\rm(\circ/s)')

figure
plot(t_imu, imu(3:2:end,4))
hold on
plot(t_imu, bias_esti(:,5))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\epsilon\it_z\rm(\circ/s)')

%% Acc bias
t_imu = (imu(3:2:end,1)-imu(1,1))/1000;

figure
plot(t_imu, imu(3:2:end,5))
hold on
plot(t_imu, bias_esti(:,6))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\nabla\it_x\rm(g)')

figure
plot(t_imu, imu(3:2:end,6))
hold on
plot(t_imu, bias_esti(:,7))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\nabla\it_y\rm(g)')

figure
plot(t_imu, imu(3:2:end,7)+1)
hold on
plot(t_imu, bias_esti(:,8))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\nabla\it_z\rm(g)')

%% Position
t_gps = (nav_gps(:,1)-imu(1,1))/1000;
t_imu = (imu(1:2:end,1)-imu(1,1))/1000;

figure
plot(t_gps, nav_gps(:,2))
hold on
plot(t_imu, nav(:,1))
plot((pos(:,11)-imu(1,1))/1000, pos(:,2))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\itL\rm(\circ)')

figure
plot(t_gps, nav_gps(:,3))
hold on
plot(t_imu, nav(:,2))
plot((pos(:,11)-imu(1,1))/1000, pos(:,3))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\lambda(\circ)')

figure
plot(t_gps, nav_gps(:,4))
hold on
plot(t_imu, nav(:,3))
plot((pos(:,11)-imu(1,1))/1000, pos(:,4))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\ith\rm(m)')

%% Velocity
t_gps = (nav_gps(:,1)-imu(1,1))/1000;
t_imu = (imu(1:2:end,1)-imu(1,1))/1000;

figure
plot(t_gps, nav_gps(:,5))
hold on
plot(t_imu, nav(:,4))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\itv_n\rm(m/s)')

figure
plot(t_gps, nav_gps(:,6))
hold on
plot(t_imu, nav(:,5))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\itv_e\rm(m/s)')

figure
plot(t_gps, nav_gps(:,7))
hold on
plot(t_imu, nav(:,6))
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\itv_d\rm(m/s)')