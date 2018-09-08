function imu = fun_parse_imu(filename)

fileID = fopen(filename, 'r');
stream = fread(fileID, 'uint8=>uint8');
fclose(fileID);

k = 1;
while 1
    if stream(k)==85 && stream(k+1)==52 && stream(k+2)==0
        stream(1:k-1) = [];
        break;
    end
    k = k+1;
end

k = 0;
while 1
    if stream(end-k-2)==85 && stream(end-k-1)==52 && stream(end-k)==0
        stream(end-k-2:end) = [];
        break;
    end
    k = k+1;
end

n = length(stream);
pn = n/52;
if mod(pn,1)~=0
    error('Lost data!');
end

for k=1:pn
    sum = uint8(0);
    for ki=1:51
        sum = bitxor(sum, stream(ki+(k-1)*52));
    end
    if sum ~= stream(k*52)
        error('Check failure!');
    end
end

imu = zeros(pn,12);
for k=1:pn
    imu(k,1) = double(typecast(stream((4:7)+(k-1)*52),'uint32'));
    imu(k,2) = double(typecast(stream((8:11)+(k-1)*52),'single'));
    imu(k,3) = double(typecast(stream((12:15)+(k-1)*52),'single'));
    imu(k,4) = double(typecast(stream((16:19)+(k-1)*52),'single'));
    imu(k,5) = double(typecast(stream((20:23)+(k-1)*52),'single'));
    imu(k,6) = double(typecast(stream((24:27)+(k-1)*52),'single'));
    imu(k,7) = double(typecast(stream((28:31)+(k-1)*52),'single'));
    imu(k,8) = double(typecast(stream((32:35)+(k-1)*52),'single'));
    imu(k,9) = double(typecast(stream((36:39)+(k-1)*52),'single'));
    imu(k,10) = double(typecast(stream((40:43)+(k-1)*52),'single'));
    imu(k,11) = double(typecast(stream((44:47)+(k-1)*52),'single'));
    imu(k,12) = double(typecast(stream((48:51)+(k-1)*52),'single'));
end

end