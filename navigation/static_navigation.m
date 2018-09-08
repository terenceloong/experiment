%% 1.Earth constant
earth_constant;

%% 2.Navigation initial value
kgps = find(info(:,2)>imu(1,1), 1) - 1;
p = pos(kgps,2:4); %deg, [lat,lon,h]
v = [0, 0, 0]; %m/s, [vn,ve,vd]
att = [0, 0, 0]; %deg, [yaw,pitch,roll]
pva0 = [p, v, att]; %record initial value

%% 3.Initialize state vector
p(1:2) = p(1:2)/180*pi; %rad
att = att/180*pi; %rad
q = angle2quat(att(1), att(2), att(3));
avp = [q, v, p]'; %column vector

%% 4.Step
dt = 2 * (imu(end,1)-imu(1,1)) / (size(imu,1)-1) / 1000;
n = floor((size(imu,1)-1)/2);

%% 5.Variables
dgyro = [0; 0; 0]; %gyro compensation, deg/s
dacc = [0; 0; 0]; %accelerometer compensation, g
g = gravity(p(1), p(3));

%% 6.Output
nav = zeros(n,9); %[lat, lon, h, vn, ve, vd, yaw, pitch, roll]
nav_gps = ones(size(info,1),9)*NaN; %[t, lat, lon, h, vn, ve, vd, clock_offset, clock_rate]

%% Filter
N = 15;
X = zeros(N,1); %[phix,phiy,phiz, dvn,dve,dvd, dlat,dlon,dh, dtr,dtv, ex,ey,ez, ax,ay,az]
X(10) = pos(kgps,5); %m
X(11) = pos(kgps,6)/1575.42e6*c; %m/s
if N == 14
%     P = diag([[1,1,1]*(2/180*pi), [1,1,1]*1, [1/a,1/a,1]*5, 10,1, [1,1,1]*(0.1/180*pi)].^2);
%     Q = diag([[1,1,1]*(0.1/180*pi)*0.707 *1, [1,1,1]*0.01*0.707 *1,...
%               [1/a,1/a,1]*0.01*0.707*dt/2 *1, 0.2*dt/2, 0.2, [1,1,1]*(0.005/180*pi)].^2) *dt^2;
elseif N == 15
    P = diag([[1,1,1]*(2/180*pi), [1,1,1]*1, [1/a,1/a,1]*5, 10,1, [1,1,1]*(0.1/180*pi), 0.02].^2);
    Q = diag([[1,1,1]*(0.1/180*pi)*0.707 *1, [1,1,1]*0.01*0.707 *1,...
              [1/a,1/a,1]*0.01*0.707*dt/2 *1, 0.2*dt/2, 0.2, [1,1,1]*(0.005/180*pi), 0.002].^2) *dt^2;
elseif N == 17
%     P = diag([[1,1,1]*(2/180*pi), [1,1,1]*1, [1/a,1/a,1]*5, 10,1, [1,1,1]*(0.1/180*pi), [1,1,1]*0.02].^2);
%     Q = diag([[1,1,1]*(0.1/180*pi)*0.707 *1, [1,1,1]*0.01*0.707 *1,...
%               [1/a,1/a,1]*0.01*0.707*dt/2 *1, 0.2*dt/2, 0.2, [1,1,1]*(0.005/180*pi)*5, [1,1,1]*0.002*5].^2) *dt^2;
end
R_rou = (0.8)^2;
R_drou = (0.1)^2;
bias_esti = zeros(n,8); %[dtr,dtv, ex,ey,ez, ax,ay,az]
filter_P = zeros(n,N); %state variable standard deviation
filter_Xc = zeros(n,N); %state variable correction

%% Ephemeris
ion_c = ion(1,2:end);
ephemeris_sv = []; %the SVs that have ephemeris
for k=1:32
    id = num2str(k);
    if exist(['ephemeris_',id], 'var')
        eval(['ephemeris_c_',id,' = ephemeris_',id,'(1,2:end)'';']);
        eval('ephemeris_sv = [ephemeris_sv, k];');
    end
end

%% 7.Calculate
kgps = kgps+1;
for k=1:n
    %--IMU data--%
    kj = 2*k+1;
    gyro0 = (imu(kj-2, 2:4)'-dgyro) /180*pi; %rad/s
    gyro1 = (imu(kj-1, 2:4)'-dgyro) /180*pi;
    gyro2 = (imu(kj  , 2:4)'-dgyro) /180*pi;
    acc0  = (imu(kj-2, 5:7)'-dacc) * g; %m/s^2
    acc1  = (imu(kj-1, 5:7)'-dacc) * g;
    acc2  = (imu(kj  , 5:7)'-dacc) * g;
    
    %--inertial navigation solution--%
%     avp = RK4(@ins_avp_qn, avp, dt, [gyro0;acc0],[gyro1;acc1],[gyro2;acc2]);
%     avp(1:4) = quatnormalize(avp(1:4)')'; %quaternion normalization
    avp = ins_solve_3(avp, dt, [gyro0;acc0],[gyro1;acc1],[gyro2;acc2]);
%     avp(1:4) = RK4(@fun_dq, avp(1:4), dt, gyro0,gyro1,gyro2);

%-----------------------------------------------------------------------------------------------%
    Phi = state_deep(avp, acc1, dt, N);
    if imu(kj,1)<info(kgps,2)
        [X, P] = kalman_filter(Phi, X, P, Q);
    elseif (info(kgps,1)>=0)==0
        [X, P] = kalman_filter(Phi, X, P, Q);
        kgps = kgps+1;
    else
        t = info(kgps,2); %time, ms
        %---update ephemeris
        for ki=1:length(ephemeris_sv)
            id = num2str(ephemeris_sv(ki));
            eval(['index = find(ephemeris_',id,'(:,1)==t, 1);']);
            if ~isempty(index)
                eval(['ephemeris_c_',id,' = ephemeris_',id,'(index,2:end)'';']);
            end
        end
        %---update ionosphere parameter
        index = find(ion(:,1)==t, 1);
        if ~isempty(index)
            ion_c = ion(index,2:end);
        end
        %---delete the SVs that don't have ephemeris
        visible_sv = find(rou(kgps,:)>0);
        for ki=1:length(visible_sv)
            id = num2str(visible_sv(ki));
            if ~exist(['ephemeris_c_',id], 'var')
                visible_sv(ki) = 0;
            end
        end
        visible_sv(visible_sv==0) = [];
        %---calulate SVs
        sv = zeros(0,8); %[x,y,z, rou, vx,vy,vz, drou]
%         clock_offset = info(kgps,3)/1000; %s
%         lat = pos(kgps,2); %deg
%         lon = pos(kgps,3); %deg
        clock_offset = X(10)/c; %s
        lat = avp(8)/pi*180; %deg
        lon = avp(9)/pi*180; %deg
        for ki=1:length(visible_sv)
            idn = visible_sv(ki); %id number
            id = num2str(idn); %id string
            eval(['sv = [sv; [sv_ecef_ephemeris(ephemeris_c_',id,', t/1000, rou(kgps,',id,'), clock_offset), -doppler(kgps,',id,')/1575.42e6*c]];']);
            eval(['tiono = Klobuchar_iono(ion_c, ele(kgps,',id,'), azi(kgps,',id,'), lat, lon, t/1000);']);
            sv(end,4) = sv(end,4) - tiono*c + clock_offset*c;
        end
        nav_gps(kgps,:) = [t,gps_pos(sv)];
        kgps = kgps+1;
        
        [H0, Z0, ng] = measure_deep(avp, sv, N);
        H = [H0; zeros(1,N)];
        H(end,3) = -1;
        R = diag([ones(1,ng)*R_rou, ones(1,ng)*R_drou, (0.01/180*pi)^2]);
        [psi,~,~] = dcm2angle(quat2dcm(avp(1:4)'));
        dpsi = angle_pmpi(psi-0);
        Z = [Z0; dpsi];
        [X, P, E, Xc] = kalman_filter(Phi, X, P, Q, H, Z, R);
        filter_Xc(k,:) = Xc';

        %---------adjust----------%
        if norm(X(1:3))>0
            phi = norm(X(1:3));
            qc = [cos(phi/2), X(1:3)'/phi*sin(phi/2)];
            avp(1:4) = quatmultiply(qc, avp(1:4)')';
        end
        avp(5:10) = avp(5:10) - X(4:9);
        X(1:9) = zeros(9,1);
        dgyro = dgyro + X(12:14)/pi*180; %deg
        if N==15
            dacc(3) = dacc(3) + X(15)/g; %g
        elseif N==17
            dacc = dacc + X(15:17)/g; %g
        end
        X(12:end) = zeros(N-11,1);
    end
    
    bias_esti(k,:) = [X(10), X(11), dgyro', dacc']; %m, m/s, deg/s, g
    filter_P(k,:) = sqrt(diag(P))';

%-----------------------------------------------------------------------------------------------%
    
    %--store--%
    nav(k,1:2) = avp(8:9)' /pi*180; %deg
    nav(k,3) = avp(10); %m
    nav(k,4:6) = avp(5:7)'; %m/s
    [r1,r2,r3] = quat2angle(avp(1:4)');
    nav(k,7:9) = [r1,r2,r3] /pi*180; %deg
end
nav = [pva0; nav];

%% Function
function Phi = state_deep(avp, fb, dt, n)
    global a f
    lat = avp(8);
    h = avp(10);
    Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5;
    Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5;
    Cbn = quat2dcm(avp(1:4)')';
    fn = antisym(Cbn*fb);
    A = zeros(n);
    A(4:6,1:3) = fn;
    if n ==14 %3-axis gyroscope bias
        A(1:3,12:14) = -Cbn;
    elseif n == 15 %3-axis gyroscope bias + z-axis accelerometer bias
        A(1:3,12:14) = -Cbn;
        A(4:6,15) = Cbn(:,3);
    elseif n ==17 %3-axis gyroscope bias + 3-axis accelerometer bias
        A(1:3,12:14) = -Cbn;
        A(4:6,15:17) = Cbn;
    end
    A(7:9,4:6) = diag([1/(Rm+h), sec(lat)/(Rn+h), -1]);
    A(10,11) = 1;
    %----------------------------------------------------%
%     global w
%     v = avp(5:7);
%     wien = [w*cos(lat); 0; -w*sin(lat)];
%     wenn = [v(2)/(Rn+h); -v(1)/(Rm+h); -v(2)/(Rn+h)*tan(lat)];
%     A(1:3,1:3) = -antisym(wien+wenn);
%     A(1:3,4:6) = [0,1/(Rn+h),0; -1/(Rm+h),0,0; 0,-tan(lat)/(Rn+h),0];
%     A(4:6,4:6) = -antisym(2*wien+wenn);
    %----------------------------------------------------%
    Phi = eye(n)+A*dt+(A*dt)^2/2;
end

function [H, Z, ng] = measure_deep(avp, sv, n)
% sv = [x,y,z, rou, vx,vy,vz, drou]
    global a f
    lat = avp(8);
    lon = avp(9);
    h = avp(10);
    Cen = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
                    -sin(lon),           cos(lon),         0;
           -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];
    ng = size(sv,1); %ng is the number of selected stars
    rpe = lla2ecef([lat/pi*180, lon/pi*180, h]); %row vector
    rpe = repmat(rpe, ng, 1); %ng row vectors
    vpe = avp(5:7)'*Cen; %row vector
    vpe = repmat(vpe, ng, 1); %ng row vectors
    rse = sv(:,1:3); %ng row vectors
    vse = sv(:,5:7); %ng row vectors
    rps = rse - rpe;
    rou = sum(rps.*rps,2).^0.5; %rou, column vector
    rpsu = rps./(rou*[1,1,1]);
    drou = sum(rpsu.*(vse-vpe),2); %drou, column vector
    He = -rpsu;
    F = [-(a+h)*sin(lat)*cos(lon), -(a+h)*cos(lat)*sin(lon), cos(lat)*cos(lon);
         -(a+h)*sin(lat)*sin(lon),  (a+h)*cos(lat)*cos(lon), cos(lat)*sin(lon);
           (a*(1-f)^2+h)*cos(lat),             0,            sin(lat)        ];
    Ha = He*F;
    Hb = He*Cen';
    H = zeros(2*ng,n);
    H(1:ng,7:9) = Ha;
    H(1:ng,10) = -ones(ng,1);
    H(ng+1:2*ng,4:6) = Hb;
    H(ng+1:2*ng,11) = -ones(ng,1);
    Z = [rou-sv(:,4); drou-sv(:,8)];
end