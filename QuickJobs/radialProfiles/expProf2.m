function [ output ] = expProf2(rInner, rOuter, rBox, sigr, npoints)
   
    % MAKES HOLLOW CHANNEL WITH EXPONENTIAL FALL-OFF WALLS
    % Version 2
    
    % general 'arb' r-profile rules:
    % 1. nothing less some threshold (about 1e-4)
    % 2. needs to end at the box radius
    % 3. does not need equal spacing 
    % 4. does not need to start at r=0
    
    % define profile
    thresh = 5e-4;
    rMin = max(rInner + sigr*log(thresh),0);
    rs = linspace(rMin, rBox, npoints);
    ns = max(thresh,min(min(1, exp(-(rInner-rs)/sigr)),exp((rOuter-rs)/sigr)));
    
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

