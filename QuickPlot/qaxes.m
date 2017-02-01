function [ axes, axlbl ] = qaxes( project, sim )
    
    % GET INPUT PARAMETERS
    rp = rpinputParser(project, sim);
    
    % AXES AND OFFSETS
    % x and y measured in [um] from the axis
    % z measured from the first bunch (forward is more negative)
    axes.X = linspace(-rp.sim.dim.x/2, rp.sim.dim.x/2, 2^rp.sim.ind.x);
    axes.Y = linspace(-rp.sim.dim.y/2, rp.sim.dim.y/2, 2^rp.sim.ind.y);
    axes.Z = linspace(-rp.beam{1}.offset.z, rp.sim.dim.z-rp.beam{1}.offset.z, 2^rp.sim.ind.z);

    % LABELS (full 6D phase space)
    [axlbl.X, axlbl.Y, axlbl.Z] = deal('x [\mum]', 'y [\mum]', 'z [\mum]');
    [axlbl.XP, axlbl.YP, axlbl.E] = deal('x'' [\murad]', 'y'' [\murad]', 'E [GeV]');
    
end

