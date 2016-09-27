function [output] = besselProf(r1, rBox, nu, npoints)
    
    function [ nfrac ] = ADKfrac( Z, xi, Elas, Tlas )
        n = neff(Z,xi); % effective principle quantum number
        w = (1.52E15)*((4^n)*xi/(n*gamma(2*n))).*((20.5*(xi^1.5)./Elas).^(2*n-1)).*exp(-6.83*(xi^1.5)./Elas); % 1/s, ionization rate
        nfrac = w.*(Tlas*1e-15); % fraction of atoms ionized
    end

    function [ Efield ] = EfromI( Ilaser )
        eps0 = 8.854E-12; % F/m
        c = 3E8; % m/s
        Efield = sqrt((2*Ilaser*1E-14)/(eps0*c)); % GV/m
        Efield=double(Efield);
    end

    function [ neff ] = neff( Z, xi )
        neff = 3.69*Z/sqrt(xi);
    end

    % input parameters
    Tlas = 100; % [fs] laser pulse duration
    xi = 5.39; % [eV] Lithium ionization threshold
    Z = 1; % first ionization only
    %r1 = 250; % [?m] radius of first peak intensity
    I0 = 6*1e12; % [W/cm^2] peak intensity
    %nu = 8; % bessel order
    
    % determine the location and height of the first bessel peak
    zs = linspace(0, 2*nu, 10000);
    [bj_maxval, bj_maxind] = max(besselj(nu, zs).^2);
    bj_maxloc = zs(bj_maxind);
    
    % bessel intensity profile
    rs = linspace(0,rBox,npoints); 
    Is = I0 * besselj(nu,rs/(r1/bj_maxloc)).^2 / bj_maxval; % [W/cm^2]
    
    % calculate density profile
    nfrac = zeros(size(rs));
    for i = 1:numel(rs)
        Elas = EfromI(Is(i));
        nfrac(i) = ADKfrac( Z, xi, Elas, Tlas );
    end
    nfrac(nfrac>1) = 1; % cap at full ionization
    
    thresh = 1e-4;
    nfrac(nfrac<thresh) = thresh;
    output = {{'arb', [rs' , nfrac']}};
    
    % plot intensity and density profiles
    yyaxis left; 
    plot(rs, Is,'r');
    xlabel('r [\mum]');
    ylabel('I [W/cm^-2]');
    
    yyaxis right;
    plot(rs, nfrac,'b');
    xlabel('r [\mum]');
    ylabel('n/n_0');
    
end