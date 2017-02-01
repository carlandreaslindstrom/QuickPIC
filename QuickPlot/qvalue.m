function [ val, label, rp ] = qvalue( project, sim, step, data, projection, x, y, z)
    
    % QuickValue: give the value of a field/force/density at a certain
    %              location, at a certain step.
    % locations are from the axis (in x/y) and from the first beam (in z)
    
    % CONSTANTS
    SI_c = 299792458; % [m/s] speed of light
    SI_e = 1.60217662e-19; % [C] electron charge
    SI_me = 9.10938356e-31; % [kg] electron mass
    SI_eps0 = 8.85418782e-12; % [m^-3 kg^-1 s^4 A^2] permittivity of free space
    
    type = data(1);
    dim = data(2);
    projection = sort(projection);
    hor_proj = projection(2);
    ver_proj = projection(1);
    
    % GET INPUT PARAMETERS
    rp = rpinputParser(project, sim);
    
    % plasma parameters
    omega_p = sqrt(rp.n0 * SI_e^2/(SI_me*SI_eps0));
    
    % axes and offsets
    axes = qaxes( project, sim );
    
    % extraction
    qp = qextract(project, sim, step);
    
    % projection image
    img_proj = qp.(type).(dim).(projection);
    
    % allow list in one dimension
    assert((numel(x)-1)*(numel(y)-1)==0);
    assert((numel(y)-1)*(numel(z)-1)==0);
    assert((numel(z)-1)*(numel(x)-1)==0);
    val_raw = zeros(1,max([numel(x), numel(y), numel(z)]));
    for i = 1:numel(x)
        for j = 1:numel(y)
            for k = 1:numel(z)
                
                % extract the (interpolated) value!
                [loc.X, loc.Y, loc.Z] = deal(x(i), y(j), z(k));
                val_raw(max([i,j,k])) = interp2(axes.(hor_proj),axes.(ver_proj),img_proj,loc.(hor_proj),loc.(ver_proj));
                
            end
        end
    end
    
    % check that slice coincides with coordinate (TODO)
    
    % add label and scaling
    switch type
        case 'E'
            label = 'E-field [GV/m]';
            scale = SI_me*SI_c*omega_p/(SI_e*1e9); % norm to GV/m
        case 'B'
            label = 'B-field [T]';
            scale = SI_me*omega_p/SI_e; % norm to T
        case 'F'
            label = 'Force [GeV/m]';
            scale = SI_me*SI_c*omega_p/(SI_e*1e9); % norm to GeV/m
        case 'Q'
            label = 'Slice density [cm^{-3}]';
            scale = rp.n0*1e-6;
    end

    % scale the value
    val = val_raw*scale;
    
    % plot if desired
    visu = false;
    if visu
        figure(5);
        set(gcf,'color','w');
    
        imagesc(axes.(hor_proj), axes.(ver_proj), img_proj);
        set(gca, 'ydir','normal');
        axis tight equal;
        
        cmax = max(max(max(img_proj)), -min(min(img_proj)));
        colormap(cmap(-cmax, cmax));
        colorbar;
        caxis([-cmax cmax]);
        
        hold on;
        plot(axes.(hor_proj)(ind.(hor_proj)), axes.(ver_proj)(ind.(ver_proj)),'k+');
        hold off;
        
        legend([label ': ' num2str(val)]);
    end
    
end

