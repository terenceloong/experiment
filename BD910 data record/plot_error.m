t = (info(:,2)-info(1,2))/1000;

%% 1.Plot drange
figure
hold on
str = {};
for k=1:32
    if isempty(find(drange(:,k)>0, 1))==0
        plot(t, drange(:,k), 'LineWidth',1.5)
        str = [str, ['SV',num2str(k)]];
    end
end
legend(str)
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\rho(m)')
title('可见卫星伪距误差')

%% 2.Plot elevation
figure
hold on
str = {};
for k=1:32
    if isempty(find(ele(:,k)>0, 1))==0
        plot(t, ele(:,k), 'LineWidth',1.5)
        str = [str, ['SV',num2str(k)]];
    end
end
legend(str)
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\beta(\circ)')
title('可见卫星高度角')

%% 3.Plot drange_rate
figure
hold on
str = {};
for k=1:32
    if isempty(find(drange_rate(:,k)>0, 1))==0
        plot(t, drange_rate(:,k), 'LineWidth',1.5)
        str = [str, ['SV',num2str(k)]];
    end
end
legend(str)
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\itd\rm\rho(m/s)')
title('可见卫星伪距率误差')