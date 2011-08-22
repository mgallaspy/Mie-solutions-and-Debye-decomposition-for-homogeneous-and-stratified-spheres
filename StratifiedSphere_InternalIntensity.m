function StratifiedSphere_InternalIntensity (n_med,ns_part,size_prms)

% ---------- Initialization -----------
% Set up the spherical coordinate system. 
% These are the points at which to evaluate the electric field.

% Number of points to sample for each coordinate.
phi_num = 1;
theta_num = 180;
rho_num = 60;

% Relative refractive index of the sphere.
% m_relative = n_particle/n_medium;

% Phi is the polar angle. 
% tan(phi) = y/x.
phi = linspace(0,pi,phi_num);

% Theta is the azimuthal angle. 
% cos(theta) = z/r.
theta = linspace(eps,2*pi,theta_num);

% Rho is the unitless radial coordinate, scaled by the radius of the sphere.
% Rho = r/a, where a is the radius of the sphere.
% r = sqrt(x^2 + y^2 + z^2).
rho = linspace(eps,1,rho_num);

% Calculate partiral wave coefficients, and then intensity at each point.
%[an,bn,cn,dn]=HomogeneousSphere_PWC(m_relative,size_prm,ind_max,orderP);
[~, ~, cnj, dnj, enj, fnj] = StratifiedSphere_PWC(size_prms, ns_part, n_med);

%[E_r,E_th,E_phi] = HomogeneousSphere_InternalField(rho,theta,phi,ind_max,n_medium,n_particle,size_prm,cn,dn,gnm_TE,gnm_TM);
[E_r, E_th, E_phi] = StratifiedSphere_InternalField(rho, theta, phi, n_med, ns_part, size_prms, cnj, dnj, enj, fnj);

Intensity=abs(E_r).^2 + abs(E_th).^2 + abs(E_phi).^2;

% ===========================  Plotting========================
% Cartesian coordinates for plotting.
x(1:length(rho),1:length(theta),1:length(phi)) = 0;
y(1:length(rho),1:length(theta),1:length(phi)) = 0;
z(1:length(rho),1:length(theta),1:length(phi)) = 0;
for rho_ind=1:length(rho)
    for th_ind=1:length(theta)
        for phi_ind=1:length(phi)
            x(rho_ind,th_ind,phi_ind) = rho(rho_ind)*sin(theta(th_ind))*cos(phi(phi_ind));
            y(rho_ind,th_ind,phi_ind) = rho(rho_ind)*sin(theta(th_ind))*sin(phi(phi_ind));
            z(rho_ind,th_ind,phi_ind) = rho(rho_ind)*cos(theta(th_ind));
        end
    end
end

% Choosing a scattering plane:
phi_plane_ind=1;

figure('position',[200,200,650,500]);
%axes1 = axes('Fontsize',14,'box','on'); %Taking reference to SU Jinming "Practicle guide to MATLAB",P206 and 207     

titlestr = sprintf('Stratified Sphere (Mie)');
title(titlestr);

xlabel('z/a','Fontsize',14,'Fontname','Times New Roman');
ylabel('x/a','Fontsize',14,'Fontname','Times New Roman');

hold on;                             % hold current figure
pcolor(z(:,:,phi_plane_ind),x(:,:,phi_plane_ind),Intensity(:,:,phi_plane_ind));
shading interp;                      % interpolate between gridpoints
axis auto;
colorbar;                            % add contour legend
colormap jet;                       % set colormap (red=high;blue=low)
hold off;                            % remove hold on figure

figure; 
h = contour(z(:,:,phi_plane_ind),x(:,:,phi_plane_ind),Intensity(:,:,phi_plane_ind),60);
%---------- End HomogeneousSphere_InternalIntensity function --------
