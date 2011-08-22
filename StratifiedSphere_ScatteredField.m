function [E_r, E_theta, E_phi] = StratifiedSphere_ScatteredField(rho, theta, phi, n_med, ...
									size_prms, an, bn)

    % Highest term calculated for the electric field series. Empirically
    % determined by Bohren and Huffman, among others.
    maxsz = max(size_prms);
    ind_max = ceil(2+maxsz+4.3*maxsz^(1/3));
                                
	%------- Calculate the theta-dependent portion of the electric field -------
    PI_Mat(1:ind_max, 1:(2*ind_max)+1, 1:length(theta)) = 0;
    TAU_Mat(1:ind_max, 1:(2*ind_max)+1, 1:length(theta)) = 0;
    for th_ind=1:length(theta)
            [PI_Mat(:,:,th_ind),TAU_Mat(:,:,th_ind)]=PI_TAU_Calc(theta(th_ind),ind_max);
    end

    % Calculate the BSCs
    [gnm_TE,gnm_TM]=BeamShapeCoeff(ind_max);
    
	% Calculate Riccati-Hankel functions for each layer of the sphere
	rho_len = length(rho);
	Xi_n(1:rho_len,1:ind_max) = 0;
	Xi_n_p(1:rho_len,1:ind_max) = 0;
    mxrho(1:rho_len) = n_med*maxsz*rho;
	for rho_ind=1:rho_len
        [Xi_n(rho_ind,:), Xi_n_p(rho_ind,:)] = HankelCalc(ind_max, mxrho(rho_ind));
	end

    % Calculate the plane wave coefficients needed in the Bromwich formulation
    n_arr=[1:1:ind_max];
    Cn_pw=[1i.^(n_arr-1)].*[(-1).^n_arr].*[(2*n_arr+1)./n_arr./(n_arr+1)];
    
    % Calculate the internal electric field
    E_r(1:length(rho),1:length(theta),1:length(phi))=0;
    E_theta(1:length(rho),1:length(theta),1:length(phi))=0;
    E_phi(1:length(rho),1:length(theta),1:length(phi))=0;
    
    for rho_ind=1:length(rho)
    	for th_ind=1:length(theta)
            for phi_ind=1:length(phi)  
                for n=1:ind_max
                    % m_arr corresponds to an array of actual m values we are interested in summing over.
                    m_arr = [-n:1:n];
                    % m_ind is a function of m_arr, which represents an array
                    % of indices which correspond to the desired m_arr values
                    % for all m-subscripted matrices (gnm_TM, PI_Mat, etc.). 
                    m_ind = m_arr + (n+1);
                    
                    % Appears in several expressions.
                    exp_imphi = exp(1i*m_arr.*phi(phi_ind));
                    
                    % Note: the PWCs in the arguments are the "plane wave"
                    % coefficients -- that is, you can always factor out
                    % the BSCs, and thus remove the m-dependence
                    
                    Er_tm = (-Cn_pw(n)*gnm_TM(n,m_ind)) .* (sin(theta(th_ind))*PI_Mat(n,m_ind,th_ind)) .* exp_imphi ...
                        .* ((an(n)*n*(n+1)*Xi_n(rho_ind,n))/(mxrho(rho_ind)^2));
                    
                    E_r(rho_ind,th_ind,phi_ind) = E_r(rho_ind,th_ind,phi_ind) + sum(Er_tm);
                        
                    
                    Eth_tm = (-Cn_pw(n)*gnm_TM(n,m_ind)) .* TAU_Mat(n,m_ind,th_ind) .* exp_imphi ...
                        .* ((an(n)*Xi_n_p(rho_ind,n))/mxrho(rho_ind));
                        
                    Eth_te = (1i*Cn_pw(n)*gnm_TE(n,m_ind)) .* PI_Mat(n,m_ind,th_ind) .* (1i*m_arr.*exp_imphi) ...
                        .* ((bn(n)*Xi_n(rho_ind,n))/mxrho(rho_ind));
                        
                    E_theta(rho_ind,th_ind,phi_ind) = E_theta(rho_ind,th_ind,phi_ind) + sum(Eth_tm) + sum(Eth_te);
                    
                    
                    Ephi_tm = (-Cn_pw(n)*gnm_TM(n,m_ind)) .* PI_Mat(n,m_ind,th_ind) .* (1i*m_arr.*exp_imphi) ...
                        .* ((an(n)*Xi_n_p(rho_ind,n))/mxrho(rho_ind));
                    
                    Ephi_te = (-1i*Cn_pw(n)*gnm_TE(n,m_ind)) .* TAU_Mat(n,m_ind,th_ind) .* (exp_imphi) ...
                        .* ((bn(n)*Xi_n(rho_ind,n))/mxrho(rho_ind));       
                    
                    E_phi(rho_ind,th_ind,phi_ind) = E_phi(rho_ind,th_ind,phi_ind) + sum(Ephi_tm) + sum(Ephi_te);
                end
            end
    	end
    end
    
    % Debugging output
    % Electric field intensity
    Intensity=abs(E_r).^2 + abs(E_theta).^2 + abs(E_phi).^2;
    max(max(max(Intensity)))
return
%End function StratifiedSphere_ScatteredField

function [R_Hankel,R_Hankel_p]=HankelCalc(ind_max,X)
    % Calculate outgoing spherical waves from orders 1 to ind_max.
	% By the Bromwich formulation, this represents an outgoing spherical wave.
    orders=[1:ind_max]; 
	R_Hankel = sqrt(X.*pi./2).*( besselj(orders+0.5,X) + i.*((-1).^(orders)).*besselj(-orders-1+0.5,X) );
	% Then calculate first derivatives from orders 1 to ind_max.
    R_Hankel_p(1:ind_max) = 0;
    for n=1:ind_max
        R_Hankel_p(n) = sqrt(X*pi/2)*(besselj(n+0.5-1,X) - (n/X)*besselj(n+0.5,X));
        R_Hankel_p(n) = R_Hankel_p(n) + i*((-1)^n)*sqrt(X*pi/2)*( besselj(-n-2+0.5,X) + ((n+1)/X)*besselj(-n-1+0.5,X));
    end
return

% ===============================================
% Function name: PI_TAU_Calc
% Authors: Feng Xu, Michael Gallaspy
% Last modified: 5/17/2011
% Inputs:
%	theta: The azimuthal angle at which the Bromwich angular functions are calculated. A scalar.
%   ind_max: maximum order term to consider in the electric field series.
%
% Outputs:
%	gnm_TE, gnm_TM: Each a vector of length ind_max. TE and TM beam shape coefficiencts for a plane wave,
%		as per the Bromwich formulation of the Mie scattering problem, appropriately truncated.
%
% Description:
%	A subroutine of the HomogeneousSphere_InternalField function.
%	Note that the subscript convention used can be a little hard to follow.
%	'n' and 'm' in this function are reserved for use as subscript values,
%		whereas 'n_ind' and 'm_ind' are reserved for use as indices which correspond to
%		the appropriate 'n' and 'm' values.
%	See paper by James A. Lock, JOSA A, Vol. 10, Issue 4, 1993.
%
%	Also note that the m=0 functions diverge for theta=0,pi.
%	In principle, the m=0 beam shape coefficients for a plane wave go to zero, and the limit of
%		the product of the m=0 BSC and the angular functions presumably is finite (and zero) as theta -> 0.
%		Since these quantities always appear as a product, there is no problem with undefined or infinite intensities.
%	In practice, MATLAB gives a divide-by-zero warning for theta=0, but does not complain for theta=pi.
% ===============================================
function [PI_nm,TAU_nm]=PI_TAU_Calc(theta,ind_max)
	% Calculate PI_nm & TAU_nm, where n ranges from 0 to ind_max, and m ranges from -ind_max to ind_max.
	% At the end, we'll drop the n=0 values and return PI_nm & TAU_nm.
	% VERY IMPORTANT: PI_nm(n_ind,m_ind) corresponds to n=n_ind-1, m=m_ind-n_ind.
	% This convention is chosen so that, for a given n-value, the m-values starting at m_ind=1 
	% correspond to -n,-n+1,-n+2,...,n-1,n. The index of a given m-value depends on the n-value!
	PI_nm(1:ind_max+2,1:2*ind_max+2) = 0;
	TAU_nm(1:ind_max+1,1:2*ind_max+1) = 0;
	
	% Initialize some known values, for m>=1.
	% See James A. Lock, JOSA A, Vol. 10, Issue 4, 1993, pp. 697.
	for m=1:ind_max
		% These indices correspond to m=n.
		n_ind = m+1;
		m_ind = 2*m+1;
		odds = 2*[1:m] - 1;
		% Note: the n-value is decreased by the decreased index, but the m-value is constant with the decreased index!
		PI_nm(n_ind-1,m_ind-1) = 0;
		PI_nm(n_ind,m_ind) = prod(odds)*( sin(theta) )^(m-1);
	end
		
	% Now use upwards recursion over n, for each positive m value.
	for m=1:ind_max
		m_start = 2*m+1;
		n_start = m+1;
		for n=m:ind_max
			offset = n-m;
			n_ind = n_start + offset;
			m_ind = m_start + offset;
			% The m index of each PI_nm entry is chosen to maintain a constant m-value.
			PI_nm(n_ind+1,m_ind+1) = ((2*n+1)*cos(theta)*PI_nm(n_ind,m_ind) - (n+m)*PI_nm(n_ind-1,m_ind-1))/(n+1-m);
			% Now we can also get TAU_nm values, using a recurrence on PI_nm.
			TAU_nm(n_ind,m_ind) = n*cos(theta)*PI_nm(n_ind,m_ind) - (n+m)*PI_nm(n_ind-1,m_ind-1);
			% Now also set negative m-values, using symmetry.
			neg_m_ind = (2*n+1) - m_ind + 1;
			PI_nm(n_ind,neg_m_ind) = PI_nm(n_ind,m_ind);
			TAU_nm(n_ind,neg_m_ind) = PI_nm(n_ind,m_ind);
		end
	end
	
	% Finally, we just need to calculate values for m=0, using associated Legendre polynomials.
	% Again, P_0n(n_ind) correponds to n=n_ind-1;
	P_0n(1:ind_max+1) = 0;
	% For the cases n=0 and n=1:
	P_0n(1) = 1; P_0n(2) = cos(theta);
	% For the case n=m=0.
	PI_nm(1,1) = P_0n(1)/sin(theta);
	for n=1:ind_max
		n_ind = n+1;
		P_0n(n_ind+1) = ((2*n+1)*cos(theta)*P_0n(n_ind) - n*P_0n(n_ind-1))/(n+1);
		% m_ind corresponds to m=0
		m_ind = n+1;
		PI_nm(n_ind,m_ind) = P_0n(n_ind)/sin(theta);
		TAU_nm(n_ind,m_ind) = n*(cos(theta)*PI_nm(n_ind,m_ind) - PI_nm(n_ind-1,m_ind-1));
	end

	% Then drop the n=0 parts of each returned matrix.
	PI_nm(1,:) = [];
	TAU_nm(1,:) = [];
    
    % And drop an extra dimension added to the end of PI_NM
    PI_nm(ind_max+1,:) = [];
    PI_nm(:,2*ind_max+1) = [];    
return

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
%	A subroutine of the Homogeneous_Sphere_Interior function.
%	See "Light scattering from a sphere arbitrarily located in a Gaussian beam, using a Bromwich formulation",
%		By Gouesbet, et al, JOSA A, Vol. 5, Issue 9, pp. 1427-1443 (1988) 
%		for more info on the Bromwich formulation.
% ===============================================
function [gnm_TE,gnm_TM]=BeamShapeCoeff(ind_max)
    %------------ Beam shape coefficients for plane wave ------------%
    % These have been calculated already for the plane wave, 
    % we merely need to assign the proper values.
    % For n, m subscriped matrices, the first index equals n, and the
    % second index maps to a particluar m.
    % E.G. gnm_TM(p,q) means n = p and m = q - (n+1).
    gnm_TM(1:ind_max,1:ind_max*2+1)=0.0;
    gnm_TE(1:ind_max,1:ind_max*2+1)=0.0;
    for n=1:ind_max
        gnm_TM(n,n+1+1)=0.5;
        gnm_TM(n,n+1-1)=0.5;    
        gnm_TE(n,n+1+1)=-0.5*1i;
        gnm_TE(n,n+1-1)=0.5*1i;
    end
    %------------ End of beam shape coefficient calculations ------------%
return