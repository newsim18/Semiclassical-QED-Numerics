function S = SBE(t,y,mq,d,Kmin,N,E0,w0,profileType)
    % System of backreaction equations for A,E with a spin-1/2 field in 1+1 dimensions

    if profileType == 1
        j0 = 2*E0*w0*sech(w0.*t).^2.*tanh(w0.*t); % SP
    else
        j0 = -E0./(1+t).^(2); % AC
    end

    S = zeros(N,1);
    % Integrand for k integral present in quantum current term in backreaction equation
    k_integrand = zeros((N-2)/4,1);

    for n=0:4:N-3
        %We give values to K and W
        K=(n*d)/4 + Kmin; % current mode
        W=(K^2+mq^2).^(1/2); % omega

        % System of equations for the modes
        S(n+1) = (- (K - y(N-1)).*y(n+2) - mq*y(n+4));  % Reh1dot
        S(n+2) = (+ (K - y(N-1)).*y(n+1) + mq*y(n+3));  % Imh1dot
        S(n+3) = (+ (K - y(N-1)).*y(n+4) - mq*y(n+2));  % Reh2dot
        S(n+4) = (- (K - y(N-1)).*y(n+3) + mq*y(n+1));  % Imh2dot
        
        k_integrand(n/4+1,1) = y(n+1).^2 + y(n+2).^2 - y(n+3).^2 - y(n+4).^2 + K./W;
    end

    % Simpson's 1/3 method for k-integral
    
    f1 = 0;
    f2 = 0;

    % Calculating the even k terms
    for j = 2:2:length(k_integrand)-1
        f1 = f1 + k_integrand(j);
    end

    % Calculating the odd k terms    
    for j = 3:2:length(k_integrand)-2
        f2 = f2 + k_integrand(j);
    end

    %Electric Field
    S(N) = j0 - y(N-1)/(pi) + (1/(2*pi))*(d/3)*(k_integrand(1) + 4*f1 + 2*f2 + k_integrand(length(k_integrand)));

    % Adot = -E
    S(N-1) = y(N);      

end
