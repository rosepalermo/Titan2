function [p,g] = TadpoleRun(p,g)

% TadpoleRun.m
%
% Performs main iteration loop of Tadpole

n=0;
p.lastsave = 0;
kill_switch = 0;
keepgoing = 1;
g.nLakeCells_init = length(find(g.U<=p.sealevel_init));
if isfield(p,'rand_gen')
    g.meanrelief = p.relief_init;
end
while keepgoing
    n = n + 1;
    p.n = n;
    %%%%%%%%%%%%%%%%%%%%%% UPDATE ELEVATIONS %%%%%%%%%%%%%%%%%%%%%%%
    
    [p,g] = Update(p,g);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%% INCREMENT TIME %%%%%%%%%%%%%%%%%%%%%%%%%
    
    p.t = p.t + p.dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if it is time to redraw the plot, do so.
    if p.doDrawPlot
        if ~rem(n,p.plotint)
            DrawPlot(n,p,g)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if p.doSaveOutput
        %         if ~rem(n,p.saveint)
        if ~p.noSLR% if sea level rise simulaiton (future runs)
            if p.n == 1 || ((sign(g.SL_slope(p.n)) == -1) && (sign(g.SL_slope(p.n-1)) == 1)) % if first ts iteration or if highstand
                %             p.lastsave = n/p.saveint + 1;
                p.lastsave = p.lastsave+1;
                g.output(:,:,p.lastsave) = g.U;
                g.t(p.lastsave) = p.t;
                g.sealevelsave(p.lastsave) = g.sealevel;
                save(p.runname, '-v7.3', 'p', 'g');
            end
        else % if no sea level rise
            if p.n == 1 || ~rem(n,p.saveint) || kill_switch
                %             p.lastsave = n/p.saveint + 1;
                p.lastsave = p.lastsave+1;
                g.output(:,:,p.lastsave) = g.U;
                if p.doWaveErosion | p.doUniformErosion
                g.adt_save(p.lastsave) = g.adt;
                end
                if p.doWaveErosion
                    if isempty(g.wave_matrix)
                        continue
                    else
                        %                     g.dam_wave_save(:,:,p.lastsave) = g.dam_wave;
                        g.wave_matrix_save(:,:,p.lastsave) = g.wave_matrix;
                    end
                end
                if p.doUniformErosion && ~p.doWaveErosion
                    if isempty(g.dam_uniform)
                        continue
                    else
                        g.dam_uniform_save(:,:,p.lastsave) = g.dam_uniform;
                    end
                end
                g.t(p.lastsave) = p.t;
                g.sealevelsave(p.lastsave) = g.sealevel;
%                 save(p.runname, '-v7.3', 'p', 'g');
                g.nLakeCells = length(find(g.U<=g.sealevel));
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kill_switch = isfield(p,'boundary');
    g.nLakeCells = length(find(g.U<=g.sealevel));
    
    if isfield(p,'rand_gen')
    temp = g.output(:,:,1)-g.U;
    g.meanrelief = mean(g.U,'all')-min(g.U,[],'all');
    keepgoing = g.meanrelief > p.relief_init*p.relief_final; %keepgoing if the relief isn't eroded past final
    elseif isfield(p,'aerrortest')
        if p.t<p.tf
            keepgoing = 1;
        else
            keepgoing = 0;
        end
    
    else
    keepgoing = (~kill_switch && g.nLakeCells < p.size_final*g.nLakeCells_init);%g.nLakeCells < p.Ao_cells*p.size_final) % if either the lake hits the bournday or the max lake size
    end
end

p.iterations = n;