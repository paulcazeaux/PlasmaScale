clear; clc; clf;

[d1, v1, p1] = DataImport('ionwave_FullPIC.dmp');

%%
[d2, v2, p2] = DataImport('ionwave_Shay_2.dmp');

[d3, v3, p3] = DataImport('ionwave_Shay_3.dmp');
[d4, v4, p4] = DataImport('ionwave_Shay_4.dmp');
% plot density

[ds1, ds3] = synchronize(d1, d3, 'uniform', 'interval', 0.01);
[ds2, ds4] = synchronize(d2, d4, 'uniform', 'interval', 0.01);


figure(1); clf;
surf(ds1.time*ones(1,size(ds1.data, 2)) , ones(length(ds1.time),1) * (1:size(ds1.data, 2))/size(ds1.data, 2), ds1.data);
hold all;
surf(ds2.time*ones(1,size(ds2.data, 2)) , ones(length(ds2.time),1) * (1:size(ds2.data, 2))/size(ds2.data, 2), ds2.data);
surf(ds3.time*ones(1,size(ds3.data, 2)) , ones(length(ds3.time),1) * (1:size(ds3.data, 2))/size(ds3.data, 2), ds3.data);
surf(ds4.time*ones(1,size(ds4.data, 2)) , ones(length(ds4.time),1) * (1:size(ds3.data, 2))/size(ds4.data, 2), ds4.data);


figure(2); clf;

m = min(min(min(d1.data)), min(min(d2.data)));
M = max(max(max(d1.data)), max(max(d2.data)));

while(true)
    for i = 1:length(ds1.time)
        clf;
        n1 = size(ds1.data, 2);
        n2 = size(ds2.data, 2);
        n3 = size(ds3.data, 2);
        n4 = size(ds4.data, 2);

        plot(1./n1*(1:n1), ds1.data(i,:), 'b');
        ylim([m, M]);
        hold all;
        plot(1./n2*(1:n2), ds2.data(i,:), 'g');
        plot(1./n3*(1:n3), ds3.data(i,:), 'r');
        plot(1./n4*(1:n4), ds4.data(i,:), 'k');
        title(['Time: ',num2str(ds1.time(i))]);
        pause;
    end
end

%% plot density error

figure(3); clf;
[ds1, ds2] = synchronize(d1, d2, 'uniform', 'interval', 0.01);
m = min(min(min(d1.data)), min(min(d2.data)));
M = max(max(max(d1.data)), max(max(d2.data)));

while(true)
    for i = 1:length(ds1.time)
        clf;
        n = size(ds1.data, 2);

        plot(1./n*(1:n), abs(ds2.data(i,:) - ds1.data(i,:)), 'r');
        title(['Time: ',num2str(ds1.time(i))]);
        pause;
    end
end

%% plot velocity

figure(2); clf;
[vs1, vs2] = synchronize(v1, v2, 'uniform', 'interval', 0.01);
m = min(min(min(v1.data)), min(min(v2.data)));
M = max(max(max(v1.data)), max(max(v2.data)));

while(true)
    for i = 1:length(vs1.time)
        clf;
        n1 = size(vs1.data, 2);
        n2 = size(vs2.data, 2);

        plot(1./n1*(1:n1), vs1.data(i,:));
        ylim([m, M]);
        hold all;
        plot(1./n2*(1:n2), vs2.data(i,:));
        title(['Time: ',num2str(ds1.time(i))]);
        pause;
    end
end

%% plot pressure

figure(2); clf;
[ps1, ps2] = synchronize(p1, p2, 'uniform', 'interval', 0.01);
m = min(min(min(p1.data)), min(min(p2.data)));
M = max(max(max(p1.data)), max(max(p2.data)));

while(true)
    for i = 1:length(ds1.time)
        clf;
        n1 = size(ps1.data, 2);
        n2 = size(ps2.data, 2);

        plot(1./n1*(1:n1), ps1.data(i,:));
        ylim([m, M]);
        hold all;
        plot(1./n2*(1:n2), ps2.data(i,:));
        title(['Time: ',num2str(ps1.time(i))]);
        pause;
    end
end

%% 

[ds1, ds2] = synchronize(d1, d2, 'uniform', 'interval', 0.01);
ds1 = getsamples(ds1, (1:21));
ds2 = getsamples(ds2, (1:21));
clf
plot(0.25*ds1.data(:,1) + .5*ds1.data(:,2) + .25*ds1.data(:,3))
hold all
plot(ds2.data(:,1))
