% ===============================================
% Function name: HomogeneousSphere_InternalIntensity
% Authors: Feng Xu, Michael Gallaspy
% Last modified: 5/18/2011
% Inputs: 
%   orderP: The desired Debye order terms used to calculate the electric
%       field. Should be an integer. Negative integers will use Mie
%       coefficients. Note there is also a maximum order calculated beyond
%       which Mie coefficients will be used. This maximum value is
%       changeable in the Sphere_Partial_Wave_Coeff function.
%   n_particle, n_medium: Complex refractive indices of the surrounding
%       medium and the homogeneous sphere, respectively.
%   size_prm: Size parameter of the sphere. That is 2*pi*a/lambda, where a
%       is the radius of the sphere, and lambda is the wavelength of incident
%       light.
%
% Outputs:
%   Produces a 3-dimensional plot of the internal electric field intensity 
%		(proportional to electric field magnitude squared) in the xz-plane.
%
% Description:
% 	Calculates internal electric field intensity per unit incoming intensity of
%		a homogeneous sphere.
% 	Convention is that negative imaginary index corresponds to absorption.
% 	Revised on the subroutine Force_Sphere_Interior for internal 
%   	electric field calculation for Michael Gallaspy (DRI).
%   Feng Xu - Jet Propulsion Lab NASA/Caltech
%   April 13 2011
%
%	Depends on the HomogeneousSphere_PWC function,
%	and HomogeneousSphere_InternalField function.
% ===============================================
function [an, bn, cn, dn] = HomogeneousSphere_InternalIntensity (orderP,n_particle,n_medium,size_prm)

% This function should be easily modified to return other values of
% interest (i.e., field components, partial wave coefficients, etc.)

% ---------- Initialization -----------
% Set up the spherical coordinate system. 
% These are the points at which to evaluate the electric field.

% Number of points to sample for each coordinate.
phi_num = 1;
theta_num = 180;
rho_num = 60;

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
[an,bn,cn,dn]=HomogeneousSphere_PWC(orderP,n_particle,n_medium,size_prm);

[E_r,E_th,E_phi] = HomogeneousSphere_InternalField(rho,theta,phi,n_medium,n_particle,size_prm,cn,dn);

% %%% Take this part out later...
% rho2 = linspace(1+eps(1),2,rho_num);
% [E_r2,E_th2,E_phi2] = HomogeneousSphere_ScatteredField(rho2,theta,phi,n_medium,n_particle,size_prm,an,bn);
% E_r = cat(1,E_r,E_r2);
% E_th = cat(1,E_th,E_th2);
% E_phi = cat(1,E_phi,E_phi2);
% rho = cat(2,rho,rho2);
% %%%%

Intensity=abs(E_r).^2 + abs(E_th).^2 + abs(E_phi).^2;


% ===========================  Plotting========================
% Cartesian coordinates for plotting.
% Initialize the matrices which will hold Cartesian coordinates.
x(1:length(rho),1:length(theta),1:length(phi)) = 0;
y(1:length(rho),1:length(theta),1:length(phi)) = 0;
z(1:length(rho),1:length(theta),1:length(phi)) = 0;
for rho_ind=1:length(rho)
    for th_ind=1:length(theta)
        for phi_ind=1:length(phi)
            % This is a transformation from spherical coordinates to
            % Cartesian coordinates.
            x(rho_ind,th_ind,phi_ind) = rho(rho_ind)*sin(theta(th_ind))*cos(phi(phi_ind));
            y(rho_ind,th_ind,phi_ind) = rho(rho_ind)*sin(theta(th_ind))*sin(phi(phi_ind));
            z(rho_ind,th_ind,phi_ind) = rho(rho_ind)*cos(theta(th_ind));
        end
    end
end

% Choosing a scattering plane:
phi_plane_ind=1;

figure('position',[200,200,650,500]);   

titlestr = sprintf('p=%i',orderP);
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

%---------- End HomogeneousSphere_InternalIntensity function --------
