function [an,bn,cn,dn]=HomogeneousSphere_PWC(order,n_part,n_med,size_prm)
    % Maximum order Debye term to calculate. For example, p_max=200 means that
    % the Debye terms up to p=200 will be calculated for each term in the electric field series.
    p_max=200;
    m_relative = n_part/n_med;
    
    ind_max = ceil(2+size_prm+4.3*size_prm^(1/3));

    % Compute both the Mie coefficients and the reflection and transmission
    % amplitudes with which to calculate the Debye coefficients.
    % Conjugate of the refractive index is taken because the subroutine
    % uses a different sign convention.
    [anMie,bnMie,cnMie,dnMie,R11a,R11b,R22a,R22b,T12a,T21a,T12b,T21b] = mieCoeffAndRTAmpl(m_relative,size_prm,ind_max);


    % Coefficients named with the letters 'a' and 'c' (anMie, cnDebye, etc.) correspond to TM spherical waves.
    % Coefficients named with the letters 'b' and 'd' (bnMie, dnDebye, etc.) correspond to TE spherical waves.
    % Also note that the indices are increased by one for the a and b terms.
    % That is, anDebye(p,:) corresponds to the (p-1)-order Debye terms and similarly for bnDebye, 
    % whereas cnDebye(p,:) corresponds to the p-order Debye terms and similarly for dnDebye.
    anDebye(1,:)=0.5*(1-R22a);
    bnDebye(1,:)=0.5*(1-R22b);
    cnDebye(1,:)=T21a;
    dnDebye(1,:)=T21b;
    if p_max>1
        for j=2:p_max
            anDebye(j,:)=-0.5*T12a.*R11a.^(j-2).*T21a;
            bnDebye(j,:)=-0.5*T12b.*R11b.^(j-2).*T21b;
            cnDebye(j,:)=R11a.^(j-1).*T21a;
            dnDebye(j,:)=R11b.^(j-1).*T21b;
        end
    end

    % Pick out the desired order of Debye terms.
    % Or pick the Mie terms if desired order isn't sensible (<0) or outside
    % of calculated range (>p_max).
    if (order < 0) || (order > p_max)
        an=anMie;
        bn=bnMie;
        cn=cnMie;
        dn=dnMie;
    elseif order == 0
        an=anDebye(1,:);
        bn=bnDebye(1,:);
        % The internal coefficients are by definition 0 for p=0, but they
        % need to be a zero vector of the correct length for computations.
        len = length(anDebye(1,:));
        cn(1:len)=0;
        dn(1:len)=0;
    else
        an=anDebye(order+1,:);
        bn=bnDebye(order+1,:);
        cn=cnDebye(order,:);
        dn=dnDebye(order,:);
    end
return
   
function [anMie,bnMie,cnMie,dnMie,R11a,R11b,R22a,R22b,T12a,T21a,T12b,T21b]= mieCoeffAndRTAmpl(m,x,ind_max)
    % The calling function uses a different sign convention, so take a
    % conjugate in order to correct for that.
    m = conj(m);
    
    % The size parameter x and y = relative refractive index * size parameter
    % will be the arguments for all our various bessel functions.
    y=m*x;
    
    for L=1:ind_max  
        % Compute the Riccati-Bessel and -Neumann functions to be used.
        J1x=sqrt(x*pi/2)*besselj(L+1/2,x); 
        N1x=sqrt(x*pi/2)*bessely(L+1/2,x);
        J1y=sqrt(y*pi/2)*besselj(L+1/2,y);
        N1y=sqrt(y*pi/2)*bessely(L+1/2,y);
        [J1px, N1px] = RiccBesselDerivatives(L,x);
        [J1py, N1py] = RiccBesselDerivatives(L,y);
        
        % Calculations for TE case
        tL1a(L)=J1x*J1py-m*J1px*J1y;
        tL2a(L)=N1x*N1py-m*N1px*N1y;
        tL3a(L)=N1x*J1py-m*N1px*J1y;
        tL4a(L)=J1x*N1py-m*J1px*N1y;
        [R11a(L),R22a(L),T21a(L),T12a(L)] = ReflTranAmpl(m,x,tL1a(L),tL2a(L),tL3a(L),tL4a(L));
        anMie(L)=tL1a(L)/(tL1a(L)+1i*tL3a(L));
        cnMie(L)=-m*1i/(tL1a(L)+1i*tL3a(L));
    
        % Calculations for TM case
        tL1b(L)=m*J1x*J1py-J1px*J1y;
        tL2b(L)=m*N1x*N1py-N1px*N1y;
        tL3b(L)=m*N1x*J1py-N1px*J1y;
        tL4b(L)=m*J1x*N1py-J1px*N1y;
        [R11b(L),R22b(L),T21b(L),T12b(L)] = ReflTranAmpl(m,x,tL1b(L),tL2b(L),tL3b(L),tL4b(L));
        bnMie(L)=tL1b(L)/(tL1b(L)+1i*tL3b(L));
        dnMie(L)=-m*1i/(tL1b(L)+1i*tL3b(L));
    end
    
    % Conjugate returned values to account for different sign convention of
    % the calling function.
    anMie=conj(anMie);bnMie=conj(bnMie);cnMie=conj(cnMie);dnMie=conj(dnMie);
    R11a=conj(R11a);R11b=conj(R11b);R22a=conj(R22a);R22b=conj(R22b);
    T12a=conj(T12a);T21a=conj(T21a);T12b=conj(T12b);T21b=conj(T21b);
return    

function [Jp, Np] = RiccBesselDerivatives(order,arg)
    Jm1 = sqrt(arg*pi/2)*besselj(order-1+1/2,arg); 
    J = sqrt(arg*pi/2)*besselj(order+1/2,arg);
    Jp = Jm1 - (order/arg)*J;
    
    Nm1 = sqrt(arg*pi/2)*bessely(order-1+1/2,arg); 
    N = sqrt(arg*pi/2)*bessely(order+1/2,arg);
    Np = Nm1 - (order/arg)*N;
return

function [R11,R22,T21,T12]=ReflTranAmpl(m,x,tL1,tL2,tL3,tL4)
    R11=(-(tL1-tL2)-1i*(tL3+tL4))/((tL1+tL2)+1i*(tL3-tL4));
    R22=(-(tL1-tL2)+1i*(tL3+tL4))/((tL1+tL2)+1i*(tL3-tL4));
    % Ref. [2] error for the calcualtion of T21_1 and T21_1 (see Manuscript VI), an additional multiplication factor (m*x^2) should be added !
    T21=-(m*x^2)*2*1i/x^2/((tL1+tL2)+1i*(tL3-tL4));
    T12=-(m*x^2)*2*1i/(m*x^2)/((tL1+tL2)+1i*(tL3-tL4));
return