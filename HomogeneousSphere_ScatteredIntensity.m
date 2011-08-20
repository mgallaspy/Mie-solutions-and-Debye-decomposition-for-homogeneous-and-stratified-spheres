% ===============================================
% Function name: HomogeneousSphere_ScatteredIntensity
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
%   Produces a 3-dimensional plot of the near-field intensity (proportional to
%   	electric field magnitude squared) in the xz-plane.
%
% Description:
% 	Calculates scattered electric field intensity per unit incoming intensity.
% 	Convention is that negative imaginary index corresponds to absorption.
% 	Revised on the subroutine Force_Sphere_Interior for internal 
%   	electric field calculation for Michael Gallaspy (DRI).
%   Feng Xu - Jet Propulsion Lab NASA/Caltech
%   April 13 2011
%
% 	Further revised to perform near-field calculation by Michael Gallaspy.
% 	May 17, 2011.
%	Depends on the HomogeneousSphere_PWC function.
%	and HomogeneousSphere_ScatteredField function.
% ===============================================
function HomogeneousSphere_ScatteredIntensity (orderP,n_particle,n_medium,size_prm)

% ---------- Initialization -----------
% Set up the spherical coordinate system. 
% These are the points at which to evaluate the electric field.

% Number of points to sample for each coordinate.
phi_num = 1;
theta_num = 180;
rho_num = 40;

% Relative refractive index of the sphere.
m_relative = n_particle/n_medium;

% Phi is the polar angle. 
% tan(phi) = y/x.
phi = linspace(0,pi,phi_num);

% Theta is the azimuthal angle. 
% cos(theta) = z/r.
theta = linspace(eps,2*pi,theta_num);

% Rho is the unitless radial coordinate, scaled by the radius of the sphere.
% Rho = r/a, where a is the radius of the sphere.
% r = sqrt(x^2 + y^2 + z^2).
rho = linspace(1,2,rho_num);

% Calculate partiral wave coefficients, and then intensity at each point.
[an,bn,cn,dn]=HomogeneousSphere_PWC(orderP,n_particle,n_medium,size_prm);
% Intensity = ScatteredIntensity(rho,theta,phi,ind_max,m_relative,size_prm,an,bn);
[E_r, E_th, E_phi] = HomogeneousSphere_ScatteredField(rho,theta,phi,n_medium,n_particle,size_prm,an,bn);

Intensity = abs(E_r).^2 + abs(E_th).^2 + abs(E_phi).^2;

% ===========================  Plotting========================
% Cartesian coordinates for plotting.
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
axes1 = axes('Fontsize',14,'box','on');     

titlestr = sprintf('p=%i, x=%d, n=%d - i*%d',orderP,size_prm,real(n_particle),imag(n_particle));
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

%---------- End Homogeneous_Sphere_Near function --------

% ===============================================
% Function name: BeamShapeCoeff
% Authors: Feng Xu, Michael Gallaspy
% Last modified: 5/17/2011
% Inputs: 
%   ind_max: maximum order term to consider in the electric field series.
%
% Outputs:
%	gnm_TE, gnm_TM: Each a vector of length ind_max. TE and TM beam shape coefficiencts for a plane wave,
%		as per the Bromwich formulation of the Mie scattering problem.
%
% Description:
%	A subroutine of the Homogeneous_Sphere_Near function.
%	See "Light scattering from a sphere arbitrarily located in a Gaussian beam, using a Bromwich formulation",
%		By Gouesbet, et al, JOSA A, Vol. 5, Issue 9, pp. 1427-1443 (1988) 
%		for more info on the Bromwich formulation.
% ===============================================
function [gnm_TE,gnm_TM]=BeamShapeCoeff(ind_max)
    %------------ Beam shape coefficients for plane wave ------------%
    % These have been calculated analytically already for the plane wave, 
    % we merely need to assign the proper values.
    % For n, m subscriped matrices, the first index equals n, and the
    % second index maps to a particluar m.
    % Specifically, gnm_TM(p,q) means n = p and m = q - (p+1).
    gnm_TM(1:ind_max,1:ind_max*2+1)=0.0;
    gnm_TE(1:ind_max,1:ind_max*2+1)=0.0;
    for n=1:ind_max
        gnm_TM(n,n+1+1)=0.5;
        gnm_TM(n,n+1-1)=0.5;    
        gnm_TE(n,n+1+1)=-0.5*i;
        gnm_TE(n,n+1-1)=0.5*i;
    end
return