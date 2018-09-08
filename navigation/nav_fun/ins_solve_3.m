function avp = ins_solve_3(avp, dt, imu0, imu1, imu2)

global a f w
q = avp(1:4);
v = avp(5:7);
v0 = v;
lat = avp(8);
lon = avp(9);
h = avp(10);
Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5;
Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5;
Cnb = quat2dcm(q');
wien = [w*cos(lat); 0; -w*sin(lat)];
wenn = [v(2)/(Rn+h); -v(1)/(Rm+h); -v(2)/(Rn+h)*tan(lat)];
winb = Cnb*(wien+wenn);  
wnbb0 = imu0(1:3) - winb;
wnbb1 = imu1(1:3) - winb;
wnbb2 = imu2(1:3) - winb;
acc0 = imu0(4:6);
acc1 = imu1(4:6);
acc2 = imu2(4:6);
dtheta1 = (wnbb0+wnbb1)*dt/4;
dtheta2 = (wnbb1+wnbb2)*dt/4;
dv1 = (acc0+acc1)*dt/4;
dv2 = (acc1+acc2)*dt/4;
q = RK4(@fun_dq, q, dt, wnbb0,wnbb1,wnbb2);
q = quatnormalize(q')';
dvc = 0.5*cross(dtheta1,dv1) + 7/6*cross(dtheta1,dv2) - 1/6*cross(dtheta2,dv1) + 0.5*cross(dtheta2,dv2);
v = v + Cnb'*(dv1+dv2+dvc) - dt*cross((2*wien+wenn),v) + dt*[0;0;gravity(lat,h)];
lat = lat + dt*(v0(1)+v(1))/2/(Rm+h);
lon = lon + dt*(v0(2)+v(2))/2/(Rn+h)*sec(lat);
h = h - dt*(v0(3)+v(3))/2;
avp = [q; v ; lat; lon; h];

end