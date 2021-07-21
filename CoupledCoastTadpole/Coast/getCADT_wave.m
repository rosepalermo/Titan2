function [adt,dam_matrix,wave_matrix,cells2trash,p] = getCADT_wave(p,lake);

%This function finds the damage associated with the maximum fetch in the
%landscape. Based on max damage, it calculates an adaptive timestep
if isfield(p,'boundary')
    adt = nan;
    dam_matrix = nan;
    return
end

[dam_matrix,wave_matrix,fetch_matrix,~,cells2trash,p] = get_dam_wave(lake,p);

% because damage is normalized to 1, we need as many time steps as if we
% round up the damage. When <1, only one time step, so not really adaptive.
% When >1, we divide by ceil(damage) so that it is never greater than 0.1.
if max(dam_matrix,[],'all')>(1/p.sensitivity)*p.strength
adt = ceil(max(dam_matrix,[],'all')./p.strength)*(p.sensitivity/100);

dam_matrix = dam_matrix./adt;
else
    adt = 1;
end

end

