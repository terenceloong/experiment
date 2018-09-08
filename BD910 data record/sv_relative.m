function svr = sv_relative(sv, p, v)
% sv = [x,y,z, vx,vy,vz]
% svr = [range, elevation-angle, azimuth, range_rate], m,deg,deg,m/s
% p = [latitude, longitude, altitude], deg,deg,m

svr = zeros(1,4);

Cen = dcmecef2ned(p(1), p(2));
rp = lla2ecef(p)'; %(ecef)
rs = sv(1:3)';
rps = rs-rp;
r = norm(rps);
svr(1) = r; %range, m
rpsu = rps/r; %unit vector (ecef)
rpsu_n = Cen*rpsu; %unit vector (ned)
svr(2) = asind(-rpsu_n(3)); %elevation-angle, deg
svr(3) = mod(atan2d(rpsu_n(2),rpsu_n(1)),360); %azimuth, deg

vs = sv(4:6)';
vp = Cen'*v';
svr(4) = (vs-vp)'*rpsu; %range_rate, m/s

end