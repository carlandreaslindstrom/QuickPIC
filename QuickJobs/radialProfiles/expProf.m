function [ output ] = expProf(rInner, rOuter, rBox, sigr, npoints)
   
    % 'arb' rules:
    % 1. nothing less some threshold (between 4e-5 and 1e-9)?
    % 2. end at the box radius
    % 3. equal spacing
    
    %thresh = 1e-4;
    rs = linspace(0, rBox, npoints);
    ns = min(1, exp(-(rInner-rs)/sigr));
    %ns(rs>=rOuter) = 0;
    
    % fix low densities
    %ns(ns<=thresh) = thresh;
    
    % set up plasma profile (points better distributed)
    %{
    rs = linspace(0, rInner, npoints-3);
    ns = min(1, exp(-(rInner-rs)/sigr));
    rs = [rs, rOuter-(rs(2)-rs(1))/2, rOuter+(rs(2)-rs(1))/2, rBox];
    ns = [ns, 1, 1e-9, 1e-9];
    %}
    
    % plot it
    visu = true;
    if visu
        figure(100);
        set(gcf, 'color','w');
        plot(rs, ns,'-o');
        axis tight;
        xlabel('r [\mum]'); ylabel('n / n_0');
        title('Radial plasma profile');
    end
    
    % give it back to the people
    output = {{'arb', [rs' , ns']}};
    
end

