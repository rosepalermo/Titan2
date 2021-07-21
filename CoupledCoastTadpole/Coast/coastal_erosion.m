function [lake,strength,p,dam_matrix,wave_matrix] = coastal_erosion(lake,fetch_on,strength,p,dam_matrix,wave_matrix,cells2trash)

% Titan analogue damage model for coastal erosion of a lake
% Rose Palermo 6-2019
% input is a matrix called "lake" -- this is a logical array of a landscape
% where 1 is lake and 0 is land. takes "lake" and either preforms uniform
% erosion on the shoreline or finds first order lakes

% if the shoreline is on the boundary, we aren't going to do coastal
% erosion because it doesn't work for waves (p.boundary is defined in
% ordershoreline_bw when the shoreline is on the boundary
if isfield(p,'boundary')
    return
end

% initialize lake, land
land = ~lake;

if fetch_on
    % calculate damage matrix -- wave
    if sum(dam_matrix,'all')==0 % if the wave matrix doesn't already exist
        [dam_matrix,wave_matrix,~,ind_sl_old,cells2trash,p] = get_dam_wave(lake,p);
    else
        [shoreline] = addidshoreline(lake,land,p); % corners and edges
        ind_sl_old = find(shoreline);
        %         ind_sl_old = get_ordered_sl(lake,p); % don't need the ordered
        %         shoreline anymore if we already have the damage
    end
else
    % calculate damage matrix -- uniform
    [dam_matrix,ind_sl_old,p] = get_dam_uniform(lake,p);
end
if any(dam_matrix<0,'all')
    pause
end
% if it hit a boundary, quit
if isfield(p,'boundary')
    erodedind_save = [];
    return
end

% avoid grid bias with adaptive timestep if needed
if max(dam_matrix,[],'all')>p.sensitivity*p.strength
    p.doAdaptiveCoastalTimeStep = 1;
end
% damage the shoreline
strength = strength - dam_matrix;

lake(find(strength<=0)) = 1;
land = ~lake;
erodedind = ind_sl_old(find(strength(ind_sl_old)<=0));

% find the points that failed during this time step and erode the
% neighbors for the amount of time since the coastline points failed
[sum_dam_excess,ind_excess,time_excess,dam_matrix] = calc_excess_dam(strength,dam_matrix,p);

while sum_dam_excess > 0
    if any(any(isinf(time_excess)))
        disp('test')
    end
    %update the lake & shoreline & find new sl indices
    % change land to lake at eroded pts
    lake(erodedind) = true;
    land = ~lake;
    [shoreline_new] = addidshoreline(lake,~lake,p); % corners and edges
    ind_sl_new = find(shoreline_new);
    ind_new = ind_sl_new(~ismember(ind_sl_new,ind_sl_old));
    
    %     lake(ind_excess) = 1;
    
    % make strength 0 for excess that I'm about to distribute
    strength(find(strength<0)) = 0;
    % damage the neighbors for the amount of time after the cell failed
    
    % find max time excess of 8 connected neighbors
    [time_dam_excess_,ind_rec_excess_,p] = pass_excess_time(time_excess,ind_new,ind_sl_new,ind_excess,lake,p);    % if it hit a boundary, quit
    if isfield(p,'boundary')
        return
    end
    
    if ~isempty(time_dam_excess_) % will be empty if a whole island is the only thing eroded because no connected shoreline
        % take the mamximum time of duplicates
        A = [ind_rec_excess_' time_dam_excess_'];
        [a,~,c] = unique(A(:,1),'rows');
        time_dam_excess = accumarray(c,A(:,2),[],@mean);
        ind_rec_excess = a(:,1);
        
        if ~isempty(ind_rec_excess) % if there are any points recieving the excess damage..
            if fetch_on
                % interpolate wave_weighting for the new shoreline cells
                wave_weighting_interp = interp_fetch_for_ind(lake,ind_new,wave_matrix); % this already includes
                wave_matrix(ind_new) = wave_weighting_interp;
                % damage for max excess time of 8 con neighbor
                dam_matrix(ind_rec_excess) = time_dam_excess.*p.Kwave.*shoreline_new(ind_rec_excess)*p.So./p.dxo*p.dx./p.Ao.*wave_matrix(ind_rec_excess);
            else
                % damage for max excess time of 8 con neighbor
                dam_matrix(ind_rec_excess) = time_dam_excess.*p.Kuniform.*shoreline_new(ind_rec_excess)*p.So./p.dxo*p.dx;
            end
            
            % damage shoreline
            strength(ind_rec_excess) = strength(ind_rec_excess) - dam_matrix(ind_rec_excess);
        end
        
        % UPDATE LAKE
        lake(find(strength<=0)) = 1;
        land = ~lake;
        ind_sl_old = ind_sl_new; % update the old shoreline
        erodedind = [erodedind;ind_rec_excess(find(strength(ind_rec_excess)<=0))];
    end
    % recalculate sum of excess damage (if the neighboring cells erode,
    % then the time will need to pass on to their neighbors in the next loop)
    % end
    [sum_dam_excess,ind_excess,time_excess,dam_matrix] = calc_excess_dam(strength,dam_matrix,p);
end

if exist('cells2trash') && ~isempty(cells2trash)
    cells2trash = sub2ind(size(lake),cells2trash(:,1),cells2trash(:,2));
    erodedind =[erodedind;cells2trash];
    strength(cells2trash) = 0;
end



if p.doStreamPower % if doing streampower, return strength to IC to represent strength of underlying rock.
    strength(find(strength<=0)) = p.strength;
end

% change land to lake at eroded pts
lake(erodedind) = true;


end