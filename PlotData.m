clear; clc; clf; close all;

[dref, vref, pref] = DataImport('TestIonWave/FullPIC.dmp');
range = 1:45;
for n=range
    [d(n), v(n), p(n)] = DataImport(sprintf('TestIonWave/EFPI%d.dmp', n));
end
%%
%[d(46), v(46), p(46)] = DataImport('TestIonWave/EFPItest.dmp');
%range = 1:46;
 %%
[drefs, ~] = synchronize(dref, d(range(1)), 'uniform', 'interval', 0.01);
[vrefs, ~] = synchronize(vref, v(range(1)), 'uniform', 'interval', 0.01);
[prefs, ~] = synchronize(pref, p(range(1)), 'uniform', 'interval', 0.01);

for n=range
    [~, ds(n)] = synchronize(dref, d(n), 'uniform', 'interval', 0.01);
    [~, vs(n)] = synchronize(vref, v(n), 'uniform', 'interval', 0.01);
    [~, ps(n)] = synchronize(pref, p(n), 'uniform', 'interval', 0.01);
end
% 
% figure(1); clf;
% surf(ds1.time*ones(1,size(ds1.data, 2)) , ones(length(ds1.time),1) * (1:size(ds1.data, 2))/size(ds1.data, 2), ds1.data);
% hold all;
% surf(ds2.time*ones(1,size(ds2.data, 2)) , ones(length(ds2.time),1) * (1:size(ds2.data, 2))/size(ds2.data, 2), ds2.data);
% surf(ds3.time*ones(1,size(ds3.data, 2)) , ones(length(ds3.time),1) * (1:size(ds3.data, 2))/size(ds3.data, 2), ds3.data);
% surf(ds4.time*ones(1,size(ds4.data, 2)) , ones(length(ds4.time),1) * (1:size(ds3.data, 2))/size(ds4.data, 2), ds4.data);

% plot density
figure(2); clf;

m = min(min(dref.data));
M = max(max(dref.data));

Nt = length(drefs.time);
for n = range
    Nt = min(Nt, length(ds(n).time));
end

dref_m = drefs.data;
while size(dref_m,2) > size(d(range(1)).data,2)
    dref_m = 0.25*(dref_m(:, 1:2:end) + [dref_m(:, 3:2:end) dref_m(:, 1)]) + 0.5*dref_m(:, 2:2:end);
end

for loop=1:5
    for i = 1:Nt
        clf;
        nx = size(dref_m, 2);

        plot(1./nx*(0:nx-1), dref_m(i,:), '-');
        ylim([m, M]);
        hold all;
        for n = [17 30 42]
            nx = size(ds(n).data, 2);
            plot(1./nx*(0:nx-1), ds(n).data(i,:), '--');
        end
        title(['Time: ',num2str(drefs.time(i))]);
        pause;
    end
end

%% Error
Nt = length(dref.time);
for n = range
    Nt = min(Nt, length(d(n).time));
end

dref_m = dref.data;
while size(dref_m,2) > size(d(range(1)).data,2)
    dref_m = 0.25*(dref_m(:, 1:2:end) + [dref_m(:, 3:2:end) dref_m(:, 1)]) + 0.5*dref_m(:, 2:2:end);
end

for n = range
    err(n) = norm(dref_m(1:Nt-1, :) - d(n).data(2:Nt,:), 2);
end

figure(5); clf;
plot(err);
%% plot density error

figure(3); clf;
[drefs, ds7] = synchronize(dref, d7, 'uniform', 'interval', 0.01);
m = min(min(dref.data));
M = max(max(dref.data));

while(true)
    for i = 1:length(ds1.time)
        clf;
        n = size(ds7.data, 2);

        plot(1./n*(0:n-1), abs(ds7.data(i,:) - drefs.data(i,2:4:512)), 'r');
        title(['Time: ',num2str(ds1.time(i))]);
        pause;
    end
end

%% plot velocity
figure(3); clf;

[vs1, vs2] = synchronize(v1, v2, 'uniform', 'interval', 0.01);
[~, vs3] = synchronize(v1, v3, 'uniform', 'interval', 0.01);
[~, vs4] = synchronize(v1, v4, 'uniform', 'interval', 0.01);
m = min(min(v1.data));
M = max(max(v1.data));

Nt = min([length(vs1.time) length(vs2.time) length(vs3.time) length(vs4.time)]  );

while(true)
    for i = 1:Nt
        clf;
        n1 = size(vs1.data, 2);
        n2 = size(vs2.data, 2);
        n3 = size(vs3.data, 2);
        n4 = size(vs4.data, 2);

        plot(1./n1*(0:n1-1), vs1.data(i,:), 'b');
        ylim([m, M]);
        hold all;
        plot(1./n2*(0:n2-1), vs2.data(i,:), 'g');
        plot(1./n3*(0:n3-1), vs3.data(i,:), 'r');
        plot(1./n4*(0:n4-1), vs4.data(i,:), 'k');
        title(['Time: ',num2str(ds1.time(i))]);
        pause;
    end
end

%% plot pressure

figure(4); clf;

[ps1, ps2] = synchronize(p1, p2, 'uniform', 'interval', 0.01);
[~, ps3] = synchronize(p1, p3, 'uniform', 'interval', 0.01);
[~, ps4] = synchronize(p1, p4, 'uniform', 'interval', 0.01);
m = min(min(p1.data));
M = max(max(p1.data));

Nt = min([length(ps1.time) length(ps2.time) length(ps3.time) length(ps4.time)]  );

while(true)
    for i = 1:Nt
        clf;
        n1 = size(ps1.data, 2);
        n2 = size(ps2.data, 2);
        n3 = size(ps3.data, 2);
        n4 = size(ps4.data, 2);

        plot(1./n1*(0:n1-1), ps1.data(i,:), 'b');
        ylim([m, M]);
        hold all;
        plot(1./n2*(0:n2-1), ps2.data(i,:), 'g');
        plot(1./n3*(0:n3-1), ps3.data(i,:), 'r');
        plot(1./n4*(0:n4-1), ps4.data(i,:), 'k');
        title(['Time: ',num2str(ds1.time(i))]);
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
