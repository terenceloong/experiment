function tiono = Klobuchar_iono(ion, E, A, lat_u, lon_u, t)
% E,A,lat_u,lon_u: deg

alpha = ion(1:4);
beta = ion(5:8);
%semi-circles
E = E/180;
A = A/180;
lat_u = lat_u/180;
lon_u = lon_u/180;
%---------------------------------------------------
psi = 0.0137/(E+0.11) - 0.022;
lat_i = lat_u + psi*cos(psi*pi);
if lat_i>0.416
    lat_i = 0.416;
elseif lat_i<-0.416
    lat_i = -0.416;
end
lon_i = lon_u + psi*sin(A*pi)/cos(lat_i*pi);
lat_m = lat_i + 0.064*cos((lon_i-1.617)*pi);
%---------------------------------------------------
AMP = sum(alpha.*[1,lat_m,lat_m^2,lat_m^3]);
if AMP<0
    AMP = 0;
end
PER = sum(beta.*[1,lat_m,lat_m^2,lat_m^3]);
if PER<72000
    PER = 72000;
end
t = 43200*lon_i + t; %local time
if t>=86400
    t = t-86400;
elseif t<0
    t = t+86400;
end
x = 2*pi*(t-50400)/PER;
%---------------------------------------------------
F = 1 + 16*(0.53-E)^3;
if abs(x)<pi/2
    tiono = F*(5e-9 + AMP*cos(x));
else
    tiono = F*5e-9;
end

end