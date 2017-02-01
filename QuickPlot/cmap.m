function map = cmap(mn,mx)

    % determine relative position of the zero
    if nargin < 1
        zero = 0.5;
    else
        zero = -mn/(mx-mn);
    end
    
    % COLORS
    if zero == 0 % no plasma
        map = interp1([0, 0.5, 1],[255 255 255; 70 172 255; 0 51 139]/256,linspace(0, 1, 500));
    elseif zero == 1 % no beam
        map = interp1([0, 0.5, 1],[255,148,0; 255 197 124; 255 255 255]/256,linspace(0, 1, 500));
    else
        map = interp1([0 mean([0,zero]) zero mean([zero,1]) 1],[255,148,0; 255 197 124; 255 255 255; 70 172 255; 0 51 139]/256,linspace(0, 1, 500));
    end
    
    % INVERSE COLORS
    %{
    if zero == 0 % no plasma
        map = interp1([0, 0.5, 1],[0 51 139; 70 172 255; 255 255 255]/256,linspace(0, 1, 500));
    elseif zero == 1 % no beam
        map = interp1([0, 0.5, 1],[255 255 255; 255 197 124; 255,148,0]/256,linspace(0, 1, 500));
    else
        map = interp1([0 mean([0,zero]) zero mean([zero,1]) 1],[0 51 139; 70 172 255; 255 255 255; 255 197 124; 255,148,0]/256,linspace(0, 1, 500));
    end
    %}
end

