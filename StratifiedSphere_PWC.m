% Calculates the Mie theory partial wave coefficients (PWCs) for a radially
% stratified sphere.
function [an, bn, cnj, dnj, enj, fnj] = StratifiedSphere_PWC(size_prms, ns_part, n_med)

    maxsz = max(size_prms);
    ind_max = ceil(2+maxsz+4.3*maxsz^(1/3)); 

    num_layers = length(size_prms);
    if num_layers ~= length(ns_part)
        error('Length of arguments "size_prms" and "ns_part" do not match.');
    end

    % References are Onofri, et al., 1995
    % and Wu & Yang, 1991
    % Are these arguments incorrect? Check back later...
    relRefs = [ns_part(2:num_layers), n_med] ./ ns_part(1:num_layers);
    xjs = ns_part .* size_prms;
    mj_xj1s = relRefs(2:num_layers) .* xjs(1:num_layers-1);
    mxjs = relRefs .* xjs;
    
    [D1x, D2x, ~, ~, psi_chi_x] = riccBessRatios(ind_max, xjs);
    [D1mx, ~, D3mx, psi_xi_mx, psi_chi_mx] = riccBessRatios(ind_max, mxjs);
    [D1mx1, D2mx1, ~, ~, psi_chi_mx1] = riccBessRatios(ind_max, mj_xj1s);
    L = num_layers;
    
    Q(1:ind_max,1:num_layers) = 0;
    K(1:ind_max,1:num_layers) = 0;
    R(1:ind_max,1:num_layers) = 0;
    H(1:ind_max,1:num_layers) = 0;
    
    K(:,1) = D1x(:,1); Q(:,1) = 0;
    H(:,1) = D1x(:,1); R(:,1) = 0;
    for j=2:num_layers
        Q(:,j) = psi_chi_mx1(:,j-1) ...
                    .* (relRefs(j)*K(:,j-1) - relRefs(j-1)*D1mx1(:,j-1)) ...
                    ./ (relRefs(j)*K(:,j-1) - relRefs(j-1)*D2mx1(:,j-1));
        K(:,j) = psi_chi_mx(:,j).*D1x(:,j)./(psi_chi_x(:,j)-Q(:,j)) - Q(:,j).*D2x(:,j)./(psi_chi_x(:,j)-Q(:,j));
        
        R(:,j) = psi_chi_mx1(:,j-1) ...
                    .* (relRefs(j-1)*H(:,j-1) - relRefs(j)*D1mx1(:,j-1)) ...
                    ./ (relRefs(j-1)*H(:,j-1) - relRefs(j)*D2mx1(:,j-1));
        H(:,j) = psi_chi_mx(:,j).*D1x(:,j)./(psi_chi_x(:,j)-R(:,j)) - R(:,j).*D2x(:,j)./(psi_chi_x(:,j)-R(:,j));
    end  
    
    % Calculate the scattered wave coefficients.
    % Note stability limitations discussed in Wu & Yang.
    an = psi_xi_mx(:,L) .* ( K(:,L) - relRefs(L).*D1mx(:,L) ) ./ ( K(:,L) - relRefs(L).*D3mx(:,L) );
    bn = psi_xi_mx(:,L) .* ( relRefs(L).*H(:,L) - D1mx(:,L) ) ./ ( relRefs(L).*H(:,L) - D3mx(:,L) );
    
    % Let's see if this fixes the internal coefficients...
    % [an,bn] = Debye_StratifiedSphere_PWC(size_prms,ns_part,n_med,-1);
    
    % Calculate internal partial wave coefficients
    cnj(1:ind_max,1:num_layers) = 0;
    dnj(1:ind_max,1:num_layers) = 0;

    % Stable reccurence formulas for ratios of Riccati-Bessel functions of
    % different arguments are given by Kaiser and Schweiger, 1993
    mLxL = mxjs(L);
    psiRat_mLxL_xL = psiRatio(ind_max, mLxL, xjs(L));
    psi_xi_mLxL = psi_xi_mx(:,L);
    
    cnj(:,L) = ((psi_xi_mLxL-an)./(psi_chi_x(:,L)-R(:,L))).*((psi_xi_mLxL).^(-1)).*psi_chi_x(:,L).*psiRat_mLxL_xL;
    dnj(:,L) = ((psi_xi_mLxL-bn)./(psi_chi_x(:,L)-Q(:,L))).*((psi_xi_mLxL).^(-1)).*psi_chi_x(:,L).*psiRat_mLxL_xL;
    
    % These values are multiplied by a factor of chi, so that the internal
    % field components can be expressed as Bessel function ratios
    cnj(:,L) = cnj(:,L) .* ( sqrt(xjs(L)*pi/2)*bessely((1:ind_max)+0.5,xjs(L)) )';
    dnj(:,L) = dnj(:,L) .* ( sqrt(xjs(L)*pi/2)*bessely((1:ind_max)+0.5,xjs(L)) )';
    
    [~, ~, ~, ~, psi_chi_mx1] = riccBessRatios(ind_max, mj_xj1s);
    [~, ~, ~, ~, psi_chi_x1] = riccBessRatios(ind_max, xjs(1:L-1));
    psiRat_mx1_x1 = psiRatio(ind_max, mj_xj1s, xjs(1:L-1));
    for jj=L:-1:2
        cnj(:,jj-1) = cnj(:,jj).*((psi_chi_mx1(:,jj-1)-R(:,jj))./(psi_chi_x1(:,jj-1)-R(:,jj-1))) ...
            .*((psi_chi_mx1(:,jj-1)).^(-1)).*psi_chi_x1(:,jj-1).*psiRat_mx1_x1(:,jj-1);
        dnj(:,jj-1) = cnj(:,jj).*((psi_chi_mx1(:,jj-1)-Q(:,jj))./(psi_chi_x1(:,jj-1)-Q(:,jj-1))) ...
            .*((psi_chi_mx1(:,jj-1)).^(-1)).*psi_chi_x1(:,jj-1).*psiRat_mx1_x1(:,jj-1);
    end
    
    enj = (-R.*cnj);
    fnj = (-Q.*dnj);
end

function psi1_on_psi2 = psiRatio(ind_max, args1, args2)

%     arg_len = length(args1);
%     psi1(1:ind_max,1:arg_len) = 0;
%     psi2(1:ind_max,1:arg_len) = 0;
%     for nn=1:ind_max
%        psi1(nn,:) = sqrt((2/pi)*args1).*besselj(nn+0.5,args1); 
%        psi2(nn,:) = sqrt((2/pi)*args2).*besselj(nn+0.5,args2);
%     end
% 
%     psi1_on_psi2 = psi1./psi2;
    
    % The problem with this approach is that computation of the Bessel
    % function of the first kind becomes unstable when the order >> than
    % the absolute value of the argument.
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
