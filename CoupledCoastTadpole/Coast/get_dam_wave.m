function [dam_matrix,wave_weight_matrix,fetch_matrix,indshoreline_ordered,cells2trash,p] = get_dam_wave(lake,p)

% initialize (in matrix units)
if ~isfield(p,'Nx')
    p.Nx = size(lake,1);
    p.Ny = size(lake,2);
end
x = (1:p.Nx)*p.dx;
y = (1:p.Ny)*p.dy;
[X,Y] = meshgrid(x,y);
wave_weight_matrix = nan(size(lake));
fetch_matrix = nan(size(lake));
dam_matrix = zeros(size(lake));
cells2trash = [];

% find the shoreline
[shoreline] = addidshoreline(lake,~lake,p);

% find number of first order lakes
[F_lake_all,~,~,~] = find_first_order_lakes(lake,p);

for ff = 1:length(F_lake_all)
    F_lake = F_lake_all{ff};
    if length(find(F_lake))<2
        continue
    end
    clearvars fetch_sl_cells indshoreline_ordered WaveArea_cell
    
    %order the shoreline and islands
    % new order the shoreline code
    [indshoreline_ocw,~,cells2trash_ff,p] = order_shoreline_bwbound(F_lake,p);
    if isfield(p,'boundary') 
        dam_matrix= [];
        wave_weight_matrix= [];
        fetch_matrix = [];
        indshoreline_ordered= [];
        cells2trash= [];
        return
    end
    cells2trash = [cells2trash;cells2trash_ff];
    
    if isempty(indshoreline_ocw) % stops eroding when shoreline hits the boundary
        break
    end
    
    for l = 1: length(indshoreline_ocw)
        indshoreline_ordered{l,1} = sub2ind(size(X),indshoreline_ocw{l}(:,1),indshoreline_ocw{l}(:,2));
        fetch_sl_cells{l,1}(:,1) = X(indshoreline_ordered{l,1});
        fetch_sl_cells{l,1}(:,2) = Y(indshoreline_ordered{l,1});
    end
    % calculate wave weighted (sqrt(F)*cos(theta-phi))
%     [WaveArea_cell,FetchArea_cell] = fetch_vis_approx(fetch_sl_cells);% first is wave, second is fetch!!
[WaveArea_cell,FetchArea_cell,~,~,~,~,~,~] = GetFetchArea(fetch_sl_cells,F_lake,p.nrays,p.delta,p.dx);    
clearvars erodedind
    %         end
    
    % Damage the shoreline
    indshoreline_ordered = cell2mat(indshoreline_ordered);
    %         [shoreline] = addidshoreline_cardonly(lake,land); % edges only

    wave_weighting = cell2mat(WaveArea_cell);
    fetch_all = cell2mat(FetchArea_cell);
    
    % when nonunique (on promontories), only will use first value unless
    % add them up. This part of the code finds duplicate indices and adds
    % up the wave & fetch areas to get waves from both sides of
    % promontories
    A = [indshoreline_ordered wave_weighting];
    [a,~,c] = unique(A(:,1),'rows');
    out_wave = [a(:,1),accumarray(c,A(:,2))];
    out_fetch = [a(:,1),accumarray(c,fetch_all)];
    
    unique_ind = out_wave(:,1);
    unique_sum_wave = out_wave(:,2);
    unique_sum_fetch = out_fetch(:,2);
    
    wave_weight_matrix(unique_ind) = unique_sum_wave;
    fetch_matrix(unique_ind) = unique_sum_fetch;
    dam = p.dt*p.Kwave*shoreline(unique_ind).*unique_sum_wave*p.So./p.dxo*p.dx./p.Ao;
    dam_matrix(unique_ind) = dam_matrix(unique_ind)+dam; % have to do this because if a shoreline cell is an island on the big one and has a second order lake in it, we don't want the damage from the second order lake to be the only damage it receives. need to add up damage from both sides.
    
end
    
end