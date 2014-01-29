function make_testmodels


%% Very simple model (just constant gradient from center to surface)
radius = linspace(0,3000,301);

vp = linspace(10,5,301)*1000;
vs = vp / sqrt(3);

rho = linspace(8,4,301)*1000;

fid = fopen('test_simple.bm', 'w');
fprintf(fid, ' F F\n 301\n');
for ilayer = 1:301
    fprintf(fid, '%f %f %f %f \n', radius(ilayer), rho(ilayer), vp(ilayer), vs(ilayer));
end
fclose(fid);
figure; hold on
plot(radius, vp, 'r')
plot(radius, vs, 'g')
plot(radius, rho, 'k')
legend('vp', 'vs', 'rho')
title('test simple.bm')
hold off


%% simple model (just constant gradient from center to surface, with one jump)
radius(1:101) = linspace(0,1000000,101);
radius(102:301) = linspace(1000000,3000000,200);

vp = linspace(10,5,301)*1000;
vs = vp / sqrt(3);

vp(1:101) = vp(1:101) + 1000;
vs(1:101) = vs(1:101) + 1000;

rho = linspace(8,4,301)*1000;

fid = fopen('test_one1st.bm', 'w');
fprintf(fid, ' F F\n 301\n');
for ilayer = 1:301
    fprintf(fid, '%f %f %f %f \n', radius(ilayer), rho(ilayer), vp(ilayer), vs(ilayer));
end
fclose(fid);
figure; hold on
plot(radius, vp, 'r')
plot(radius, vs, 'g')
plot(radius, rho, 'k')
legend('vp', 'vs', 'rho')
title('test one1st.bm')
hold off

%% less simple model (three layers, one fluid in the middle)
radius(1:101) = linspace(0,1000000,101);
radius(102:201) = linspace(1000000,2000000,100);
radius(202:301) = linspace(2000000,3000000,100);

vp = linspace(10,5,301)*1000;
vs = vp / sqrt(3);

vp(1:101) = vp(1:101) + 1000;
vs(1:101) = vs(1:101) + 1000;
vp(102:201) = vp(102:201) ;
vs(102:201) = 0; %vs(102:201) + 1000;
vp(202:301) = vp(202:301) + 2000;
vs(202:301) = vs(202:301) + 2000;


rho = linspace(8,4,301)*1000;

fid = fopen('test_two1st_fluid.bm', 'w');
fprintf(fid, ' F F\n 301\n');
for ilayer = 1:301
    fprintf(fid, '%f %f %f %f \n', radius(ilayer), rho(ilayer), vp(ilayer), vs(ilayer));
end
fclose(fid);
figure; hold on
plot(radius, vp, 'r')
plot(radius, vs, 'g')
plot(radius, rho, 'k')
legend('vp', 'vs', 'rho')
title('test two 1st fluid.bm')
hold off

%% not simple model (three layers, one fluid in the middle, 2nd order discs everywhere)
radius(1:101) = linspace(0,1000000,101);
radius(102:201) = linspace(1000000,2000000,100);
radius(202:301) = linspace(2000000,3000000,100);

vp = linspace(10,5,301)*1000;
vs = vp / sqrt(3);

vp(1:100) = vp(1:100) + 3000;
vs(1:100) = vs(1:100) + 3000;
vp(102:201) = vp(102:201) ;
vs(102:201) = 0; %vs(102:201) + 1000;
vp(203:301) = vp(203:301) + 3000;
vs(203:301) = vs(203:301) + 3000;


rho = linspace(8,4,301)*1000;

fid = fopen('test_two1st_two2nd.bm', 'w');
fprintf(fid, ' F F\n 301\n');
for ilayer = 301:-1:1
    fprintf(fid, '%f %f %f %f \n', radius(ilayer), rho(ilayer), vp(ilayer), vs(ilayer));
end
fclose(fid);
figure; hold on
plot(radius, vp, 'r')
plot(radius, vs, 'g')
plot(radius, rho, 'k')
legend('vp', 'vs', 'rho')
title('test two 1st two 2nd fluid.bm')
hold off

end 