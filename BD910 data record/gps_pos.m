function pos = gps_pos(sv)
% sv = [x,y,z, rou, cx,vy,vz, drou]
% pos = [lat, lon, h, vn, ve, vd, clock_offset, doppler], deg,...,ms,Hz

if size(sv,1)<4
    pos = ones(1,8)*NaN;
    return;
end

pos = zeros(1,8);
n = size(sv,1);

G = ones(n,4)*-1;
R = sv(:,4); %rou, m
V = sv(:,8); %drou, m/s
x0 = [0;0;0;0];
S = zeros(n,1);

cnt = 0;
while 1
    for k=1:n
        G(k,1:3) = sv(k,1:3)-x0(1:3)';
        G(k,1:3) = G(k,1:3)/norm(G(k,1:3));
        S(k) = G(k,1:3)*sv(k,1:3)';
    end
    x = (G'*G)\G'*(S-R);
    cnt = cnt+1;
    if cnt == 10
        error('Position iteration exceeds the threshold!');
    end
    if norm(x-x0)<1e-6
        break;
    end
    x0 = x;
end
pos(1:3) = ecef2lla(x(1:3)'); %[lat, lon, h], deg
pos(7) = x(4)/299792458*1000; %clock_offset, ms

for k=1:n
    G(k,1:3) = sv(k,1:3)-x(1:3)';
    G(k,1:3) = G(k,1:3)/norm(G(k,1:3));
    S(k) = G(k,1:3)*sv(k,5:7)';
end
v = (G'*G)\G'*(S-V);
Cen = dcmecef2ned(pos(1), pos(2));
pos(4:6) = (Cen*v(1:3))'; %[vn, ve, vd], m/s
pos(8) = v(4)/299792458*1575.42e6; %doppler, Hz

end