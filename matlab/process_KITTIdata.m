dir = 'dir_to_the_gps_way_point_file';
M   = csvread([dir,'gps_way_point.csv']);
nM  = size(M,1);
totalLength = 0;
gps_data_interp = [];

for i = 1:nM-1
    latitude1 = M(i,1);
    longitude1= M(i,2);
    latitude2 = M(i+1,1);
    longitude2= M(i+1,2);
    dist = CoordinatesToMeters(latitude1,longitude1,latitude2,longitude2);
    if dist > 1
        np = floor(dist);
        x  = [latitude1,latitude2];
        dm = latitude2 - latitude1;
        y  = [longitude1,longitude2];
        xq = latitude1:dm/np:latitude2;
        yq = interp1(x,y,xq);
        cur =[xq',yq'];
        gps_data_interp = [gps_data_interp;cur];
    end
    totalLength = totalLength+ dist;
end
plot(gps_data_interp(:,1),gps_data_interp(:,2),'ro')
dlmwrite([dir,'gps_prior_map.csv'],gps_data_interp,'precision', 9);
