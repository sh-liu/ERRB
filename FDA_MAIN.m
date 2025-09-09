clc; 
clear ; 
close all ;

%% ------------------------------------------Radar parameter settings

j = sqrt(-1);
c = 3e8;                                                % Speed of light
M = 8;                                                  % The number of tramsmit antennas
f0 = 10e9;                                              % Reference carrier frequency
f = 10e9;
fs = 20e9;                                                
lamda0 = c/f0;                                          % Reference carrier wavelength
dt = lamda0/2;                                          % Transmit element spacing
Dt = [0,1,2,3,4,5,6,7]*dt;                              % Transmit antenna spacing configuration
delta_f = 150e3;                                        % Unit frequency offset
Delta_f = [0,1,3,4,6,7,9,10]*delta_f;                   % Transmit antenna frequency offset configuration
K = 1000;                                               % Snapshot count 
SNR = 20;                                               % SNR

%% -------------------------------------------Parameter settings of targets

target = 2;                                             % Target count
R1 = 200;                                               % Range of target 1
theta1 = (-30)/180*pi;                                  % Angle of target  
R2 = 600;                                               % Range of target 2
theta2 = (30)/180*pi;                                   % Angle of target 2

%% --------------------------------------------Receive signal model

noise = 1/sqrt(2)*(randn(M,K)+j*randn(M,K));
a1 = steer_vector(f0,Delta_f,Dt,theta1,R1);             % Steering vector of Target 1
a2 = steer_vector(f0,Delta_f,Dt,theta2,R2);             % Steering vector of Target 2
sig = [a1 a2]*[sqrt(10^(SNR/10))*exp(1i*2*pi*f*(0:(K-1))/fs);sqrt(10^(SNR/10))*exp(1i*2*pi*f*(0:(K-1))/fs)];
data = sig+noise;                                       % Receive signal

%% --------------------------------------------Reference covariance matrix

Rx = data*data'./K;                                     % Covariance matrix
Rm = toeplitz(Rx(:,1)');
Gu=[1 1 1 1 1 1 1 1;                                    % Mask matrix
    1 1 0 1 0 1 0 1;
    1 0 1 1 1 1 1 1;
    1 1 1 1 0 1 0 1;
    1 0 1 0 1 1 1 1;
    1 1 1 1 1 1 0 1;
    1 0 1 0 1 0 1 1;
    1 1 1 1 1 1 1 1];
R_ref = Rm.*Gu;                                         % Reference covariance matrix

%% --------------------------------------------CVX

% gu=size(Gu,1);
% delta=0.004;                                           % Regularization parameter \delta 
% cvx_begin sdp quiet
%     variable Ru_re(gu,gu) hermitian  complex
%     minimize (1/2*square_pos(norm(Ru_re.*Gu-R_ref,'fro'))+delta*(trace(Ru_re)))
%     subject to 
%         Ru_re == hermitian_semidefinite(gu)
% cvx_end

%% --------------------------------------------ADMM

gu=size(Gu,1);
A = eye(gu,gu);
U = eye(gu,gu);
delta = 0.004;                                            % Regularization parameter \delta 
rho = 0.045;                                              % Regularization parameter \rho
Ru_re = zeros(gu,gu);
for i = 1:1000
    Ru_re = (R_ref.*Gu+rho.*((1/rho).*A+U))./(rho+Gu);
    [EV1,D1] = eig((1/rho).*(-A-delta*eye(8,8))+Ru_re);
    D1 = real(D1);
    for jj = 1:length(D1)
        if D1(jj,jj) < 0
            D1(jj,jj) = 0;
        end
    end
    U = EV1*D1*EV1';
    A = A+rho.*(U-Ru_re);
end

%% --------------------------------------------2D-MUSIC

Rmax = c/(2*delta_f);                                    % Maximum detection range
Df = 1;                                                  % Angle search interval
Dr = 2;                                                  % Range search interval
theta = (-90:Df:90)*pi/180; 
R = (0:Dr:Rmax); 
P = zeros(length(theta),length(R)); 

[EV,D] = eig(Ru_re);                  
EVA = diag(D)';                      
[EVA,I] = sort(EVA);                 
EV = fliplr(EV(:,I));                 
En = EV(:,target+1:M);              

for n = 1 : length(theta)
   for m = 1 : length(R)
        a = steer_vector(f0,Delta_f,Dt,theta(n),R(m)); 
        J = a'*En*En'*a;
        P(n,m) =1/J;
   end
end
P=P';

%% --------------------------------------------Figure

figure(); 
imagesc(theta*180/pi,R,10*log(abs(P)/max(max(abs(P)))));
colormap(jet(1e+6));
xlabel('DoA [deg]'); ylabel('Range [m]'); 
axis tight; axis xy;
hold on
scatter([theta1,theta2]*180/pi,[R1,R2],50,'black','^');



