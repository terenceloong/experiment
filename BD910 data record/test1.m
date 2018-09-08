% clc;

n = size(pos,1);

pos_ecef = zeros(n,3);
pos_ecef(1,:) = lla2ecef(pos(1,2:4));
for k=2:n
    pos_ecef(k,:) = lla2ecef(pos(k,2:4))-pos_ecef(1,:);
end
pos_ecef(1,:) = [0,0,0];

pos_geog = zeros(n,3);
Ceg = dcmecef2ned(pos(1,2), pos(1,3));
for k=1:n
    pos_geog(k,:) = (Ceg*pos_ecef(k,:)')';
end

figure
plot3(pos_geog(:,1),pos_geog(:,2),pos_geog(:,3))
hold on
plot3(0,0,0, 'o')
grid on
axis equal

figure
plot(pos_geog(:,1),pos_geog(:,2))
hold on
plot(0,0, 'o')
grid on
axis equal