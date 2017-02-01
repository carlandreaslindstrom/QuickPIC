function [ output ] = trapProf(rInner, rOuter, rBox, rthick)
    
    % TRAPEZOIDAL HOLLOW CHANNEL RADIAL PROFILE
    % r < (rInner-rthick) : n = 0
    % r > (rInner-rthick) => r < rInner : n = 0->1 (linear upramp)
    % r > rInner => r < rOuter : n = 1
    % r > rOuter => r < (rOuter+rthick) : n = 1->0 (linear downramp)
    % r > (rOuter+rthick) : n = 0
    % USAGE: rpinputMaker(..., {trapProf(...)}, ...}
    
    % define profile
    rs = [max(0,rInner-rthick), rInner, rOuter,  min(rBox,rOuter+rthick), rBox];
    ns = [0,1,1,0,0];
    
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
    output = {'arb', [rs' , ns']};
    
end

