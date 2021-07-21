function [p,g] = uniformerosion(p,g)

g.uniform_input = g.U<g.sealevel;
p.dt_save = p.dt;

g.nLakeCells = length(find(g.uniform_input));
if p.doAdaptiveCoastalTimeStep
    i = 0;
    [adt,dam_matrix,p] = getCADT_uniform(p,g.uniform_input);
    p.dt = p.dt_save./adt; %adaptive time step
    p.dt_test = 0;
    
    while i<adt
        i = i+1;
        p.dt_test = p.dt*adt;
        [g.uniform_output,g.Strength,p,dam_matrix,~] = coastal_erosion(g.uniform_input,0,g.Strength,p,[],[],[]);
        % make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
        % should be.
        erodedind = find(g.uniform_output - g.uniform_input);
        if ~isempty(erodedind)
%             disp('test')
        end
        g.uniform_input = g.uniform_output; % update input for the loop!
        g.U(erodedind) = g.sealevel-0.5;
    end
    p.dt = p.dt_save; %return the time step to what it was.
else
    % call uniform erosion function, output updates U and strength
    [g.uniform_output,g.Strength,p,dam_matrix,~] = coastal_erosion(g.uniform_input,0,g.Strength,p,[],[],[]);
    erodedind = find(g.uniform_output - g.uniform_input);
    
    % make new subaqueous points elevation 0.5 below sea level? Talk to Andrew about what this depth
    % should be.
    % g.U(erodedind) = g.sealevel(p.n)-p.deptheroded;
    g.U(erodedind) = g.sealevel-p.deptheroded;
end

g.nLakeCells = length(find(g.uniform_output));
g.dam_uniform = dam_matrix;
g.adt = adt;

end