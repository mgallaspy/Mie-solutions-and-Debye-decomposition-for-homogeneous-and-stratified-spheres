function [an, bn] = Debye_StratifiedSphere_PWC(size_prms, ns_part, n_med, orderP)
%   Size parameters and indices of refraction should be specified in order
%   starting with the core and ending with the last layer

    maxsz = max(size_prms);
    ind_max = ceil(2+maxsz+4.3*maxsz^(1/3));  

    throw_away = 0;

    ns = [ns_part, n_med];

    num_layers = length(size_prms);
    if num_layers ~= length(ns_part)
        error('Length of arguments "size_prms" and "ns_part" do not match.');
    end
    
    % Compute basic amplitudes for each interface, TE case
    Nbase(1:ind_max,1:num_layers) = 0;
    Dbase(1:ind_max,1:num_layers) = 0;
    Pbase(1:ind_max,1:num_layers) = 0;
    Qbase(1:ind_max,1:num_layers) = 0;
    for L=1:num_layers
        x = size_prms(L)*ns(L+1);
        y = size_prms(L)*ns(L);
        n1 = ns(L);
        n2 = ns(L+1);
        [psi_x,chi_x,throw_away,psi_x_p,chi_x_p,throw_away] = riccBessels(ind_max,x);
        [psi_y,chi_y,throw_away,psi_y_p,chi_y_p,throw_away] = riccBessels(ind_max,y);
        Nbase(:,L) = (n1*psi_x.*psi_y_p - n2*psi_x_p.*psi_y);
        Dbase(:,L) = (n1*chi_x.*psi_y_p - n2*chi_x_p.*psi_y);
        Pbase(:,L) = (n1*psi_x.*chi_y_p - n2*psi_x_p.*chi_y);
        Qbase(:,L) = (n1*chi_x.*chi_y_p - n2*chi_x_p.*chi_y);
    end
    
    % Then calculate the composite amplitudes
    [Ncomp, Dcomp, Pcomp, Qcomp] = compositeAmplitudes(Nbase,Dbase,Pbase,Qbase,num_layers);
    
    % Ans use the composite amplitudes to calculate reflection and
    % transmission amplitudes.
    T_M1 = 2i*prod(ns(2:num_layers+1))./(Ncomp + Qcomp + 1i*(Dcomp - Pcomp));
    T_1M = 2i*prod(ns(1:num_layers))./(Ncomp + Qcomp + 1i*(Dcomp - Pcomp));
    R_M1M = (-Ncomp + Qcomp + 1i*(Dcomp + Pcomp))./(Ncomp + Qcomp + 1i*(Dcomp - Pcomp));
    R_1M1 = (-Ncomp + Qcomp - 1i*(Dcomp + Pcomp))./(Ncomp + Qcomp + 1i*(Dcomp - Pcomp));
    
    % Then get the TE scattering coefficient.
    if orderP < 0
        bn = 0.5*(1 - R_M1M - T_M1.*T_1M./(1-R_1M1));
    elseif orderP == 0
        bn = 0.5*(1 - R_M1M);
    elseif orderP == 1
        bn = T_M1.*T_1M;
    else
        bn = T_M1.*T_1M.*(R_1M1).^(orderP-1);
    end
    
    % ------------------------- TM Case ----------------------
    % Compute basic amplitudes for each interface, TM case
    Nbase(1:ind_max,1:num_layers) = 0;
    Dbase(1:ind_max,1:num_layers) = 0;
    Pbase(1:ind_max,1:num_layers) = 0;
    Qbase(1:ind_max,1:num_layers) = 0;
    for L=1:num_layers
        x = size_prms(L)*ns(L+1);
        y = size_prms(L)*ns(L);
        n1 = ns(L);
        n2 = ns(L+1);
        [psi_x,chi_x,throw_away,psi_x_p,chi_x_p,throw_away] = riccBessels(ind_max,x);
        [psi_y,chi_y,throw_away,psi_y_p,chi_y_p,throw_away] = riccBessels(ind_max,y);
        Nbase(:,L) = (n2*psi_x.*psi_y_p - n1*psi_x_p.*psi_y);
        Dbase(:,L) = (n2*chi_x.*psi_y_p - n1*chi_x_p.*psi_y);
        Pbase(:,L) = (n2*psi_x.*chi_y_p - n1*psi_x_p.*chi_y);
        Qbase(:,L) = (n2*chi_x.*chi_y_p - n1*chi_x_p.*chi_y);
    end
    
    % Then calculate the composite amplitudes
    [Ncomp, Dcomp, Pcomp, Qcomp] = compositeAmplitudes(Nbase,Dbase,Pbase,Qbase,num_layers);
    
    % Ans use the composite amplitudes to calculate reflection and
    % transmission amplitudes.
    T_M1 = 2i*prod(ns(2:num_layers+1))./(Ncomp + Qcomp + 1i*(Dcomp - Pcomp));
    T_1M = 2i*prod(ns(1:num_layers))./(Ncomp + Qcomp + 1i*(Dcomp - Pcomp));
    R_M1M = (-Ncomp + Qcomp + 1i*(Dcomp + Pcomp))./(Ncomp + Qcomp + 1i*(Dcomp - Pcomp));
    R_1M1 = (-Ncomp + Qcomp - 1i*(Dcomp + Pcomp))./(Ncomp + Qcomp + 1i*(Dcomp - Pcomp));
    
    % Then get the TE scattering coefficient.
    if orderP < 0
        an = 0.5*(1 - R_M1M - T_M1.*T_1M./(1-R_1M1));
    elseif orderP == 0
        an = 0.5*(1 - R_M1M);
    elseif orderP == 1
        an = T_M1.*T_1M;
    else
        an = T_M1.*T_1M.*(R_1M1).^(orderP-1);
    end
end

function [Ncomp, Dcomp, Pcomp, Qcomp] = compositeAmplitudes(Nbase,Dbase,Pbase,Qbase,num_layers)

    % Initialize composite values for the iterative computation
    Ncomp = Nbase(:,1);
    Dcomp = Dbase(:,1);
    Pcomp = Pbase(:,1);
    Qcomp = Qbase(:,1);
    
    for L=1:num_layers-1
        Ntemp = Ncomp;
        Dtemp = Dcomp;
        Ptemp = Pcomp;
        Qtemp = Qcomp;
        Ncomp = Dtemp.*Nbase(:,L+1) - Ntemp.*Pbase(:,L+1);
        Dcomp = Dtemp.*Dbase(:,L+1) - Ntemp.*Qbase(:,L+1);
        Pcomp = Qtemp.*Nbase(:,L+1) - Ptemp.*Pbase(:,L+1);
        Qcomp = Qtemp.*Dbase(:,L+1) - Ptemp.*Qbase(:,L+1);
    end    
end

function [psi, chi, xi, psi_p, chi_p, xi_p] = riccBessels(ind_max, args)
    jlength = length(args);
    psi(1:ind_max+1,1:jlength) = 0;
    chi(1:ind_max+1,1:jlength) = 0;
    xi(1:ind_max+1,1:jlength) = 0;
    psi_p(1:ind_max,1:jlength) = 0;
    chi_p(1:ind_max,1:jlength) = 0;
    xi_p(1:ind_max,1:jlength) = 0;
    
    orders = [0:ind_max]';
    
    for j=1:jlength
        % Calculate 0th through ind_max-th order Riccati-Bessels
        psi(:,j) = sqrt(args(j)*pi/2)*besselj(orders+0.5,args(j));
        chi(:,j) = sqrt(args(j)*pi/2)*bessely(orders+0.5,args(j));
        xi(:,j) = psi(:,j) - 1i*chi(:,j);
        
        % Calculate 1st through ind_max-th order Riccati-Bessels
        psi_p(:,j) = psi(1:ind_max,j) - (2:ind_max+1)'.*psi(2:ind_max+1,j)/args(j);
        chi_p(:,j) = chi(1:ind_max,j) - (2:ind_max+1)'.*chi(2:ind_max+1,j)/args(j);
        xi_p(:,j) = xi(1:ind_max,j) - (2:ind_max+1)'.*xi(2:ind_max+1,j)/args(j);
    end
    
    % Remove 0th order Riccati-Bessels from returned arrays
    psi(1,:) = [];
    chi(1,:) = [];
    xi(1,:) = [];
end

