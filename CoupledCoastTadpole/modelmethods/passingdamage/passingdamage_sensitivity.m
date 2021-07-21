% measuring the sensitivity to resolution and time step

% take the total strength loss and divide it by the number of shoreline
% cells to get a noramlized strength loss

% input shoreline that is not gridded

load('fetch_input8con','fetch_sl_cells'); % x and y are ordered clockwise if i increases downward, first point != last point
x = fetch_sl_cells{1}(:,1);
y = fetch_sl_cells{1}(:,2);

p.Kuniform = 0.0000123;
p.Kwave = 0.00001;
p.con8 = 1;
p.So = 1;
p.dxo = 100;
p.doStreamPower = 0;
p.doLandslides = 0;
p.doUniformErosion = 1;
p.doWaveErosion = 0;

% grid shoreline and calculate strength loss over 100 years with time steps
% varying from 5 to 100 years
% dx = [5 18 45 90];
dx = 18;
dt = 100;%[10 20 50 100];
tend = 150000;
ncells = 50;%[20 50 100];
sidelength = 100;
totalSloss_frac = zeros(length(ncells),length(dt));
for i = 1:length(ncells)
        [lake,~,~] = gridlake(flipud(x)',flipud(y)',dx(i),dx(i),400);
%     [init,~] = test_circle(p,dx(i),200);lake = ~init;
%     [init,~] = test_square(p,dx(i),900);lake = ~init;
%     [lake,dx(i)] = test_square_repelem(ncells(i),sidelength);
    for t = 1:length(dt)
        [strengthloss{i,t},totalSloss_frac(i,t),test{i,t}] = calcstrengthloss(lake,dt(t),tend,p,dx(i),0);
    end
end
%%
figure()
imagesc(dt,ncells,totalSloss_frac*100)
colormap gray
colorbar
xlabel('dt')
ylabel('dx')
set(gca,'FontSize',12)

figure()
plot(ncells,totalSloss_frac(:,1)*100)
xlabel('ncells')
ylabel('total Strength loss/ domain initial strength')

%%

for i = 1:length(strengthloss)
size_init(i) = size(strengthloss{i},1)^2;
end
Sinit = p.So*dx.^2./p.dxo^2;
% Ssum = Sinit.*size_init;

for i = 1:length(strengthloss)
total_strengthloss(i) = sum(strengthloss{i},'all');
end


for i = 1:length(strengthloss)
shorelinelength(i) = length(find(strengthloss{i}>0));
end
figure()
for i = 1:length(strengthloss)
    test{i,1}(test{i,1}==0) = nan;
    plot(test{i,1})
    hold on
end
ylabel('sum(strength_o_u_t-strength_i_n)')
xlabel('time')

figure()
imagesc(strengthloss{1})

%total strength loss
% circle mean pass = 89.4615, max = 89.4805, no pass = 86.5286
