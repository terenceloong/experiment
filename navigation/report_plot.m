%% Gyro bias
t_imu = (imu(3:2:end,1)-imu(1,1))/1000;

figure
plot(t_imu, imu(3:2:end,2))
hold on
grid on
plot(t_imu, bias_esti(:,3), 'LineWidth',1.2)
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\epsilon\it_x\rm(\circ/s)')
title('x����������ƫ��������')
legend('���ٶ�','��ƫ����')

figure
plot(t_imu, imu(3:2:end,3))
hold on
grid on
plot(t_imu, bias_esti(:,4), 'LineWidth',1.2)
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\epsilon\it_y\rm(\circ/s)')
title('y����������ƫ��������')
legend('���ٶ�','��ƫ����')

figure
plot(t_imu, imu(3:2:end,4))
hold on
grid on
plot(t_imu, bias_esti(:,5), 'LineWidth',1.2)
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\epsilon\it_z\rm(\circ/s)')
title('z����������ƫ��������')
legend('���ٶ�','��ƫ����')

%% Gyro bias
t_imu = (imu(3:2:end,1)-imu(1,1))/1000;

figure
plot(t_imu, bias_esti(:,3), 'LineWidth',1.2)
grid on
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\epsilon\it_x\rm(\circ/s)')
title('x����������ƫ��������')

figure
plot(t_imu, bias_esti(:,4), 'LineWidth',1.2)
grid on
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\epsilon\it_y\rm(\circ/s)')
title('y����������ƫ��������')

figure
plot(t_imu, bias_esti(:,5), 'LineWidth',1.2)
grid on
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\epsilon\it_z\rm(\circ/s)')
title('z����������ƫ��������')

%% Acc bias
t_imu = (imu(3:2:end,1)-imu(1,1))/1000;

figure
plot(t_imu, bias_esti(:,6), 'LineWidth',1.2)
grid on
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\nabla\it_x\rm(g)')
title('x����ٶȼ���ƫ��������')

figure
plot(t_imu, bias_esti(:,7), 'LineWidth',1.2)
grid on
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\nabla\it_y\rm(g)')
title('y����ٶȼ���ƫ��������')

figure
plot(t_imu, bias_esti(:,8), 'LineWidth',1.2)
grid on
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\nabla\it_z\rm(g)')
title('z����ٶȼ���ƫ��������')

%% Velocity
t_gps = (nav_gps(:,1)-imu(1,1))/1000;
t_imu = (imu(1:2:end,1)-imu(1,1))/1000;

figure
plot(t_gps, nav_gps(:,5), 'LineWidth',1.2)
hold on
grid on
plot(t_imu, nav(:,4), 'LineWidth',1.2)
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\itv_n\rm(m/s)')
title('�����ٶ�����')
legend('GPS����','��ϵ������')

figure
plot(t_gps, nav_gps(:,6), 'LineWidth',1.2)
hold on
grid on
plot(t_imu, nav(:,5), 'LineWidth',1.2)
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\itv_e\rm(m/s)')
title('�����ٶ�����')
legend('GPS����','��ϵ������')

figure
plot(t_gps, nav_gps(:,7), 'LineWidth',1.2)
hold on
grid on
plot(t_imu, nav(:,6), 'LineWidth',1.2)
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\itv_d\rm(m/s)')
title('�����ٶ�����')
legend('GPS����','��ϵ������')

%% Attitude
t_imu = (imu(1:2:end,1)-imu(1,1))/1000;

figure
plot(t_imu, nav(:,7), 'LineWidth',1.2)
grid on
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\psi(\circ)')
title('��ϵ������������')

figure
plot(t_imu, nav(:,8), 'LineWidth',1.2)
grid on
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\psi(\circ)')
title('��ϵ�������������')

figure
plot(t_imu, nav(:,9), 'LineWidth',1.2)
grid on
set(gca, 'xlim', [t_imu(1),t_imu(end)])
xlabel('\itt\rm(s)')
ylabel('\psi(\circ)')
title('��ϵ�����ת������')

%% SV number
t_gps = (info(:,2)-info(1,2))/1000;

figure
plot(t_gps, info(:,4), 'LineStyle','none', 'Marker','.', 'MarkerSize',10)
set(gca, 'xlim', [t_imu(1),t_imu(end)])
set(gca, 'ylim', [min(info(:,4))-1,max(info(:,4))+1])
xlabel('\itt\rm(s)')
ylabel('(��)')
title('�ɼ�������Ŀ')