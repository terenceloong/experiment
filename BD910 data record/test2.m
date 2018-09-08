% Initialize ephemeris and ionosphere parameter
ion_c = ion(1,2:end);
ephemeris_sv = []; %the SVs that have ephemeris
for k=1:32
    id = num2str(k);
    if exist(['ephemeris_',id], 'var')
        eval(['ephemeris_c_',id,' = ephemeris_',id,'(1,2:end)'';']);
        eval('ephemeris_sv = [ephemeris_sv, k];');
    end
end

%-------------------------------------------------------------------------%
n = min(size(rou,1), size(pos,1));

% measured value 
range = ones(n,32)*NaN;
range_rate = -(doppler(1:n,:)+pos(1:n,6)*ones(1,32))/1575.42e6*299792458;

% theoretical value
ele0 = ones(n,32)*NaN;
azi0 = ones(n,32)*NaN;
range0 = ones(n,32)*NaN;
range_rate0 = ones(n,32)*NaN;

%-------------------------------------------------------------------------%

for k=1:n
    t = info(k,2); %time
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
    visible_sv = find(rou(k,:)>0);
    for ki=1:length(visible_sv)
        id = num2str(visible_sv(ki));
        if ~exist(['ephemeris_c_',id], 'var')
            visible_sv(ki) = 0;
        end
    end
    visible_sv(visible_sv==0) = [];
    
    %-----------------------------------------------------------------------------------------------%
    for ki=1:length(visible_sv)
        idn = visible_sv(ki); %id number
        id = num2str(idn); %id string
        
        eval(['sv = sv_ecef_ephemeris(ephemeris_c_',id,', t/1000, rou(k,',id,'), info(k,3)/1000);']);
        svr = sv_relative([sv(1:3),sv(5:7)], pos(k,2:4), [0,0,0]);
        range0(k,idn) = svr(1); %range
        ele0(k,idn) = svr(2); %elevation-angle
        azi0(k,idn) = svr(3); %azimuth
        range_rate0(k,idn) = svr(4); %range_rate
%         tiono = Klobuchar_iono(ion_c, svr(2), svr(3), pos(k,2), pos(k,3), t/1000);
        eval(['tiono = Klobuchar_iono(ion_c, ele(k,',id,'), azi(k,',id,'), pos(k,2), pos(k,3), t/1000);']);
        range(k,idn) = sv(4) - tiono*299792458;
    end
    %-----------------------------------------------------------------------------------------------%
end

drange = range-range0;
drange_rate = range_rate - range_rate0;