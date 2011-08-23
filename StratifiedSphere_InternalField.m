function [E_r, E_theta, E_phi] = StratifiedSphere_InternalField(rho, theta, phi, n_med, ...
									ns_part, size_prms, cnj, dnj, enj, fnj)
                                
    % Highest term calculated for the electric field series. Empirically
    % determined by Bohren and Huffman, among others.
    maxsz = max(size_prms);
    ind_max = ceil(2+maxsz+4.3*maxsz^(1/3));    

    % Calculate the BSCs
    [gnm_TE,gnm_TM]=BeamShapeCoeff(ind_max);
                                    
	%------- Calculate the theta-dependent portion of the electric field -------
    PI_Mat(1:ind_max, 1:(2*ind_max)+1, 1:length(theta)) = 0;
    TAU_Mat(1:ind_max, 1:(2*ind_max)+1, 1:length(theta)) = 0;
    for th_ind=1:length(theta)
            [PI_Mat(:,:,th_ind),TAU_Mat(:,:,th_ind)]=PI_TAU_Calc(theta(th_ind),ind_max);
    end

	% Calculate radially-dependent values for each layer of the sphere
	rho_len = length(rho);
    mxrho(1:rho_len) = 0;
	for rho_ind=1:rho_len
		% Determine which layer we are in...
        layer = getLayer(rho(rho_ind), size_prms);
		size_prm = size_prms(layer);
        n_part = ns_part(layer);
        mxrho(rho_ind) = size_prm*n_part*rho(rho_ind);
	end

    % Calculate the plane wave coefficients needed in the Bromwich formulation
    n_arr=[1:1:ind_max];
    Cn_pw=[1i.^(n_arr-1)].*[(-1).^n_arr].*[(2*n_arr+1)./n_arr./(n_arr+1)];
    
    % Calculate the internal electric field
    E_r(1:length(rho),1:length(theta),1:length(phi))=0;
    E_theta(1:length(rho),1:length(theta),1:length(phi))=0;
    E_phi(1:length(rho),1:length(theta),1:length(phi))=0;
    
    % Relative refractive indices, with the index of the medium tacked onto
    % the end in order to have a consistent indexing convention.
    relRefs = [ns_part, n_med] ./ n_med;
    
    for rho_ind=1:length(rho)
        layer = getLayer(rho(rho_ind), size_prms);
    	for th_ind=1:length(theta)
            for phi_ind=1:length(phi)
                % Need to express Riccati-Bessel functions as numerically
                % stable ratios.
                [D1_mxrho, D2_mxrho, ~, ~, psi_chi_mxrho] = riccBessRatios(ind_max, mxrho(rho_ind));
                [~, ~, ~, ~, psi_chi_xL] = riccBessRatios(ind_max, max(size_prms)*max(ns_part));
                psi1_psi2 = psiRatio(ind_max, mxrho(rho_ind), max(size_prms)*max(ns_part));
                chi_psi_mxrho = 1./psi_chi_mxrho;
                
                chi_ratio = chi_psi_mxrho .* psi1_psi2 .* psi_chi_xL;
                chi_p_ratio = D2_mxrho .* chi_ratio;
                if layer == 1
                   chi_ratio(1:ind_max) = 0;
                   chi_p_ratio(1:ind_max) = 0;
                end
                
                psi_ratio = psi1_psi2 .* psi_chi_xL;
                psi_p_ratio = D1_mxrho .* psi_ratio;
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
                    Er_tm = (relRefs(layer+1)*Cn_pw(n)*gnm_TM(n,m_ind)) .* (sin(theta(th_ind))*PI_Mat(n,m_ind,th_ind)) .* exp_imphi ...
                        .* ((cnj(n,layer)*n*(n+1)*psi_ratio(n) + enj(n,layer)*n*(n+1)*chi_ratio(n))/(mxrho(rho_ind)^2));
                    
                    E_r(rho_ind,th_ind,phi_ind) = E_r(rho_ind,th_ind,phi_ind) + sum(Er_tm);
  
                        
                    Eth_tm = (relRefs(layer+1)*Cn_pw(n)*gnm_TM(n,m_ind)) .* TAU_Mat(n,m_ind,th_ind) .* exp_imphi ...
                        .* ((cnj(n,layer)*psi_p_ratio(n) + enj(n,layer)*chi_p_ratio(n))/mxrho(rho_ind));
                        
                    Eth_te = (-1i*relRefs(layer+1)*Cn_pw(n)*gnm_TE(n,m_ind)) .* PI_Mat(n,m_ind,th_ind) .* (1i*m_arr.*exp_imphi) ...
                        .* ((dnj(n,layer)*psi_ratio(n) + fnj(n,layer)*chi_ratio(n))/mxrho(rho_ind));
                        
                    E_theta(rho_ind,th_ind,phi_ind) = E_theta(rho_ind,th_ind,phi_ind) + sum(Eth_tm) + sum(Eth_te);
                    
                    
                    Ephi_tm = (relRefs(layer+1)*Cn_pw(n)*gnm_TM(n,m_ind)) .* PI_Mat(n,m_ind,th_ind) .* (1i*m_arr.*exp_imphi) ...
                        .* ((cnj(n,layer)*psi_p_ratio(n) +  enj(n,layer)*chi_p_ratio(n))/mxrho(rho_ind));
                    
                    Ephi_te = (1i*relRefs(layer+1)*Cn_pw(n)*gnm_TE(n,m_ind)) .* TAU_Mat(n,m_ind,th_ind) .* (exp_imphi) ...
                        .* ((dnj(n,layer)*psi_ratio(n) +  fnj(n,layer)*chi_ratio(n))/mxrho(rho_ind));   
                    
                    E_phi(rho_ind,th_ind,phi_ind) = E_phi(rho_ind,th_ind,phi_ind) + sum(Ephi_tm) + sum(Ephi_te);
                end
            end
    	end
    end
end
%End function StratifiedSphere_InternalField

function layer = getLayer(rho, size_prms)
    % Get largest size parameter (i.e., the last element)
    max_size = size_prms(length(size_prms));
    layer = 0;
	for layer_iter=1:length(size_prms)
        if rho <= size_prms(layer_iter)/max_size
			layer = layer_iter;
			break;
		end
    end
end

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
end

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
end

function psi1_on_psi2 = psiRatio(ind_max, args1, args2)
    % See Kaiser and Schweiger, 1993, for recurrence formula
    % args1 and args2 should be the same length
    jjlength = length(args1);
    psi1_on_psi2(1:ind_max, 1:jjlength) = 0;
    psi_start(1:jjlength) = 1;
    
    % Find some auxiliary values for initialization of recurrence
    A = exp(imag(2*args1 + args2)) + exp(imag(args2));
    B = exp(imag(2*args1 + args2)) - exp(imag(args2));
    C = exp(imag(2*args2 + args1)) + exp(imag(args1));
    D = exp(imag(2*args2 + args1)) - exp(imag(args1));
    
    psi_start(:) = ( A.*sin(real(args1)) + 1i*B.*cos(real(args1)) ) ...
        ./ ( C.*sin(real(args2)) + 1i*D.*cos(real(args2)) );
    
    % Compute logarithmic derivatives
    [D1x1, ~, ~, ~, ~] = riccBessRatios(ind_max, args1);
    [D1x2, ~, ~, ~, ~] = riccBessRatios(ind_max, args2);
    
    psi1_on_psi2(1,:) = psi_start.*((D1x2(1,:)+1./args2)./(D1x1(1,:)+1./args1));
    for nn=2:ind_max
       psi1_on_psi2(nn,:) = psi1_on_psi2(nn-1,:).*((D1x2(nn,:)+nn./args2)./(D1x1(nn,:)+nn./args1));
    end
end

function [D1, D2, D3, psi_on_xi, psi_on_chi] = riccBessRatios(ind_max, args)   
    
    % Initialize to avoid many dynamic allocations.
    jlength = length(args);
    D1(1:ind_max,1:jlength) = 0;
    D2(1:ind_max,1:jlength) = 0;
    D3(1:ind_max,1:jlength) = 0;
    psi_on_xi(1:ind_max,1:jlength) = 0;
    psi_on_chi(1:ind_max,1:jlength) = 0; 
    
    % For each argument, compute D1 by downward recurrence.
    % Start the downward recurrence a few terms above the desired ind_max.
    % Justification to start 15 terms above desired index comes from Bohren
    % and Huffman.
    D1start(1:jlength)=0;
    nDelta = 60;
    for n=(ind_max+nDelta):-1:ind_max
        D1start = n./args - (D1start + n./args).^(-1); 
    end
    D1(ind_max,:)=D1start;
    for n=ind_max:-1:2
        D1(n-1,:) = n./args - (D1(n,:) + n./args).^(-1);
    end
    
    % Then let's calculate D3 by upwards recurrence.
    D3(1,:) = (1./args - 1i).^(-1) - 1./args;
    for n=2:ind_max
        D3(n,:) = (n./args - D3(n-1,:)).^(-1) - n./args;
    end
    
    % Then we get psi_on_xi by upwards recurrence
    psi_on_xi(1,:) = (sin(args)./(sin(args)-1i*cos(args))).* ...
                        (D3(1,:) + 1./args) ./ (D1(1,:) + 1./args);
    for n=2:ind_max
        psi_on_xi(n,:) = psi_on_xi(n-1,:) .* (D3(n,:) + n./args) ./ (D1(n,:) + n./args);   
    end
    
    % D2 also calculated by upward recurrence...
    D2(1,:) = (1./args + tan(args)).^(-1) - 1./args;
    for n=2:ind_max
        D2(n,:) = (n./args - D2(n-1,:)).^(-1) - n./args;
    end
    
    % Then we can get psi_on_chi by upwards recurrence
    psi_on_chi(1,:) = tan(args) .* (D2(1,:) + 1./args) ./ (D1(1,:) + 1./args);
    for n=2:ind_max
        psi_on_chi(n,:) = psi_on_chi(n-1,:) .* (D2(n,:) + n./args) ./ (D1(n,:) + n./args);
    end
end
