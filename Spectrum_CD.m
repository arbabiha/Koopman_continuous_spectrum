% Code for the paper "Data-driven spectral analysis of the Koopman operator" by Milan Korda, Mihai Putinar and Igor Mezic
% Preprint: https://arxiv.org/pdf/0806.1528.pdf

% Milan Korda, May 2018

clear all
close all

addpath('./tools')

%% Define system

N = 100; % Number of fourier moments
%data = 'artificial'; % Mix of continuous and discrete spectrum
data = 'catMap'; % Purely continuous spectrum
%data = 'artificial';
%data = 'cavity';
%data = 'lorenz'
contTime = 0;
switch data
    case 'artificial'
        % Artificial data - atomic + AC + singular continuous
        cen = [0, 0.2, 0.6, 0.8 1 ]; % Dirac centers
        weights = [0.05 0.1 0.1 0.1 0.05]; % Weights
        % Compute moments
        y = zeros(N+1,1);
        y_atoms = getDiracMoments_fourier( cen, N, weights );
        y_cantor = getMoments_cantor_fourier( N );
        y_leb = 4*getLebesgueMoments_fourier( N, 0.3, 0.7 );
        y = y + y_atoms;
        y = y + y_leb;
        y = y + y_cantor;
    case 'catMap'
        % Dynamics
        T  = @(x)( mod([2*x(1,:)+x(2,:);x(1,:)+x(2,:)],1) );
        
        % Observable 1
        psi = @(x)(exp(1i*2*pi*(2*x(1,:)+x(2,:))) + 0.5*exp(1i*2*pi*(5*x(1,:)+3*x(2,:))));
        rho_f = @(theta)(5/4 + cos(2*pi*theta));
        
        % Observable 2
        %psi = @(x)(exp(1i*2*pi*(2*x(1,:)+x(2,:))) + 0.5*exp(1i*2*pi*(5*x(1,:)+3*x(2,:))) + 0.25*exp(1i*2*pi*(13*x(1,:)+8*x(2,:))));
        %rho_f = @(theta)(21/16 + (5/4)*cos(2*pi*theta) + 0.5*cos(4*pi*theta));
        
        x0 = [0.132;0.567];
    case 'lorenz'
        contTime = 1;
        Ts = 0.2; % Sampling interval
        sig = 10; beta = 8/3; rho = 28;
        lorenz = @(t,x)([sig*(x(2,:)-x(1,:)) ; x(1,:).*(rho-x(3,:)) - x(2,:) ; x(1,:).*x(2,:) - beta*x(3,:) ] );
        psi = @(x)(x(3,:));
        
    case 'cavity'
        Re = 30;
        %addpath('./Cavity_flow_Hassan/ToolBox&Data/')
        disp(['KMD of cavity flow at Re=',num2str(Re)])
        load('cavity_obs',['obs_raw_Re' num2str(Re) 'k'])
        eval(['obs_raw = obs_raw_Re' num2str(Re) 'k'])
        mv = mean(obs_raw);
        pp = (max(obs_raw) - min(obs_raw))/2;
        Psi = (obs_raw -  mv) / pp ;
        M = size(Psi,2) - N;
end

fprintf('Spectrum approximation for %s data \n',data)

%% Generate data

switch data
    case 'catMap'
        M = 1e5; % Simulation length
        disp('....Generating data for Cat map....')
        X = x0;
        for i = 1:M+N
            X = [X T(X(:,end))];
        end
    case 'lorenz'
        M = 1e6; % Simulation length
        if(M == 1e6) % Preload the data for M = 1e6 (takes some time to simulate)
            disp('....Loading data for Lorenz...')
            load('lorenzDataTs02M1e6','X')
        else
            disp('....Generating data for Lorenz...')
            x0 = [0.0142;0.0422;0.0916];
            [tt,X] = ode45(lorenz,0:Ts:(M+N)*Ts,x0);
            X = X';
        end
end
% Create observable from data
if(~exist('Psi','var'))
    Psi = psi(X);
end



%% Hankel DMD
disp('....Computing DMD operator....')
Data_delayEmbed = hankel(Psi(1:N),Psi(N:end));

Ylift = Data_delayEmbed(:,2:end);
Xlift = Data_delayEmbed(:,1:end-1);

XXt = Xlift * Xlift';
YXt = Ylift * Xlift';
A_DMD = YXt * pinv(XXt);
eigDMD = eig(A_DMD);




%% Compute Fourier moments
if(~strcmp(data,'artificial'))
    disp('....Computing moments....')
    y = zeros(N+1,1);
    for i = 0:N
        y(i+1) = (Psi(1:M)*Psi(1+i:M+i)')/M; % compute moments <psi,U^i psi> using the ergodic averages
    end
end


%% Compute the Christoffel-Darboux kernel
disp('....Computing CD kernel....')

% Normalize
if(strcmp(data,'lorenz') || strcmp(data,'cavity'))
    y = y / y(1); % normalize to a probability measure
end

% Offset
offset_const = 1;
y(1) = y(1) + offset_const; % Prevents M from being singular if the measure is discrete

% Msoment matrix
Mtilde = momentMat_fourier(y,0); % \tilde{M} in the paper

% Evaluate CD kernel on a grid
Minv = inv(Mtilde);
x = 0:0.0001:1; x_cd = x;
c = evalChristoffelPol(Minv,x); % Christoffel-Darboux kernel evaluated at x
if(norm(imag(c)) > 1e-12)
    warning('Norm of the imaginary part of the CD kernel is large!');
end
c = real(c);
spectrum_CD = (N+1)./c - offset_const; % -offset_const to subtract the constant density we added before

% Plot approximation of absolutely continuous part of the spectrum
plot(x,spectrum_CD,'linewidth',2); hold on
switch data
    case 'artificial'
        rho = zeros(1,numel(x));
        rho(x >= 0.3 & x <= 0.7) = 4;
        stairs(x,rho,'--r','linewidth',2)
        legend('spectrum approximation','continuus density','atom locations')
        title('CD kernel - density')
    case 'catMap'
        clf
        plot(x,rho_f(x),'-r','linewidth',2); hold on
        plot(x,spectrum_CD,'linewidth',2); hold on
        LEG = legend('True','Computed');
        LEG.FontSize = 20;
        title('Spectral density of CatMap for observable #2','Fontsize',15)
    case 'lorenz'
        axis([0,1,0,5])
end
xlabel('x','fontsize',20)
ylabel('Spectral density','fontsize',20 )

% Plot approximation of point spectrum
figure
plot(x,spectrum_CD/(N+1),'linewidth',2 ); hold on
title('CD kernel - point spectrum')
if(strcmp(data,'artificial'))
    stem(cen,weights,'--r')
end

%% Zeros of orthogonal polynomials and eigvals of Hankel DMD
yy = y;
yy(1) = yy(1) - offset_const; % Get rid of the offset constant - otherwise the zeros of OPUC and Hankel DMD eigenvalues wouldn't match
Mmu = momentMat_fourier(yy,0); % M without tilde in the paper
try
Linv = inv(chol(Mmu).');
V = evalTrigBasis_comExp( x, N );
P = Linv*V;

coef_N = Linv(end,:); % coefficients of Nth orthonorm poly  (index from lowest to highest degree)
coef_N = coef_N / coef_N(end); % make it monic
R_opuc = roots(coef_N(end:-1:1)); % Zeros of N-th OPUC polynomial

if(exist('eigDMD','var'))
    figure
    plot(real(eigDMD),imag(eigDMD),'o','color','blue'); hold on
    plot(real(R_opuc),imag(R_opuc),'x','color','red')
    plot(cos([0:0.001:2*pi]),sin([0:0.001:2*pi]),'-b')
    title('Roots of OPUC and of the Hankel DMD operator')
end
catch
    disp('Roots of OPUC could not be computed because the moment matrix is not psd due to numerical errors.')
end


%% Cesaro mean approximation to CDF
V = evalTrigBasis_comExp( x, N ); V = V.'; Vc = conj(V);

S = cell(1,N+1);
CM = 0;
for i = 1:N+1
    S{i} = V(:,1:i)*conj(y(1:i)) + Vc(:,2:i)*y(2:i);
    CM = CM + S{i}/(N+1);
end
CM = CM - offset_const;

IFT = V*conj(y) + Vc(:,2:end)*y(2:end) - offset_const;

%% Singularity indicator
figure
dconst = 1/(x(2)-x(1));
F_CS = cumsum(CM)*(x(2)-x(1));
Delta = dconst*(diff(cumsum(CM+dconst)*(x(2)-x(1))).'./diff(cumsum(spectrum_CD+dconst)*(x_cd(2)-x_cd(1))) - 1);
Delta = (1/F_CS(end))*Delta;
plot(x(1:end-1),Delta,'linewidth',2); hold on
LEG = legend('Singularity indicator'); LEG.FontSize = 20;

%% Quadrature (needs yalmip to formulate the Linear program)

if(exist('yalmiptest','file') == 2)
    
    % Moments
    yy = y; yy(1) = yy(1) - offset_const;
    
    % Grid
    nq = 10*N;
    xx = 0:1/nq:1;
    
    %DFT matrix
    M_DFT = exp(1i*2*pi*[0:N]'*xx);
    
    % Compute weights using LS
    w_comp = M_DFT\yy;
    
    
    % Compute weights usin LP
    w_pos = sdpvar(numel(xx),1);
    solvesdp([w_pos>=0,sum(w_pos)==yy(1), w_pos(1) == w_pos(end)],norm(M_DFT*w_pos - yy,'fro'))
    w_pos = double(w_pos);
    
    figure
    stem(xx,double(w_pos));
    title('Atomic approximation - weights')
    
    if(strcmp(data,'artificial'))
        cdf = 0;
        xx1 = 0:0.001:1;
        for i = 2:numel(xx1)
            cdf_inc = 0;
            t = xx1(i);
            for k = 1:numel(cen)
                if(abs(t-cen(k)) < 1e-12)
                    cdf_inc = cdf_inc + weights(k);
                end
            end
            if(t>= 0.4 && t < 0.7);
                cdf_inc = cdf_inc + 4*(abs(xx1(i)-xx1(i+1)));
            end
            cdf = [cdf cdf(end)+cdf_inc];
        end
        figure
        [xcant, cdf_cant] = cantor_function(10);
        plot(xcant,cdf_cant,'-r','linewidth',3); hold on
    end
    plot(xx,cumsum(w_pos),'--b','linewidth',2);hold on
    plot(x,cumsum(CM)*(x(2)-x(1)),'--g','linewidth',2); hold on
    plot(x,cumsum(spectrum_CD)*(x_cd(2)-x_cd(1)),'--m','linewidth',2)
    if(strcmp(data,'artificial'))
        LEG = legend('True','Quadrature (LP)','Cesaro', 'N/K(z,z)');
    else
        LEG = legend('Quadrature (LP)','Cesaro', 'N/K(z,z)');
    end
    LEG.FontSize = 20;
    title('Approximation of the cumulative density function','Fontsize',15)
    
else
    warning('Yalmip not found. If you want to to get the quadrature results, install yamip from https://yalmip.github.io/download/')
    figure
    plot(x,cumsum(CM)*(x(2)-x(1)),'--g','linewidth',2); hold on
    plot(x,cumsum(spectrum_CD)*(x_cd(2)-x_cd(1)),'--m','linewidth',2)
    LEG = legend('Cesaro', 'N/K(z,z)'); LEG.FontSize = 20;
    title('Approximation of the cumulative density function','Fontsize',15)
end




