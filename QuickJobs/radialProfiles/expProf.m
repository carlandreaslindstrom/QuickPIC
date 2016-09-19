function [ output ] = expProf(rmax, r1, sigr, npoints, isarb)
    
    % set up plasma profile
    rs = linspace(0, rmax, npoints+1);
    ns = min(1, exp(-(r1-rs)/sigr));
    
    % apply threshold
    thresh = 0.01;
    ns(ns<thresh) = 0;
    
    % plot it
    visu = true;
    if visu
        set(gcf, 'color','w');
        bar(rs, ns,1);
        axis tight;
        xlabel('r [\mum]'); ylabel('n / n_0');
        title('Radial plasma profile');
    end
    
    % convert to 'arb' profile
    arb = {{'arb', [rs' , ns']}};
    
    % convert to multi-species
    species = cell(1,npoints);
    for i = 1:npoints
        if i == 1 % innermost is a cylinder, not a hollow channel
            species{i} = {'cyl', [rs(i+1), ns(i)]};
        else
            species{i} = {'hc', [rs(i), rs(i+1), ns(i)]};
        end
    end
    
    if exist('isarb','var') && isarb
        output = arb;
    else
        output = species;
    end
    
end

