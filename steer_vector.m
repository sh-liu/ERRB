function a = steer_vector(f0,Delta_f,Dt,theta,R)
% steer_vector function: steering vector for FDA radar
% f0 is the reference carrier frequency, Dt is the transmit antenna spacing configuration,
% Delta_f is the frequency offset setting between adjacent transmit antenna, 
% theta and R are the azimuth angle and range, respectively.

j = sqrt(-1);
c = 3e8;
a_t_r = exp(-j*2*pi*2*Delta_f'*R/c);            
a_t_theta = exp(j*2*pi*f0/c*Dt'*sin(theta));        
a = a_t_r.*a_t_theta;

end
