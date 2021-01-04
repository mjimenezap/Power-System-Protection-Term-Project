
clear all
close all
clc

%% Reading the files
%file_dat_event_01 = readtable('ECE6323_PROJECT_02_CTCHANNEL_EVENT_01_MAIN.dat');

file_dat_event_01 = readtable('ECE6323_PROJECT_02_CTCHANNEL_EVENT_02_MAIN.dat');

% %% Information from cfg file
% % Channel 1, identifier "Burden_Voltage", no info about channel or circuit
% % monitoring, voltage units are V, there is offset. 
% %Channel multipplier
% v_multiplier = 0.000038357;
% %Channel offset
% v_offset = -0.204407170; 
% 
% %The nominal frequency of the system is 60 Hz. 
% frequency = 60;
% %There is only one sampling rate, which is 4799.9846400492Hz. There are 7200 samples.
% sampling_rate = 4799.9846400492;
% n_samples = 4700;
% 
% sampling_period = 1/sampling_rate; %In seconds

%% Information from cfg file
% Channel 1, identifier "Burden_Voltage", no info about channel or circuit
% monitoring, voltage units are V, there is offset. 
%Channel multipplier
v_multiplier = 0.000167598;
%Channel offset
v_offset = -0.731290042; 

%The nominal frequency of the system is 60 Hz. 
frequency = 60;
%There is only one sampling rate, which is 4799.9846400492Hz. There are 7200 samples.
sampling_rate = 4799.9846400492;
n_samples = 4700;

sampling_period = 1/sampling_rate; %In seconds

%% Data conversion
time_vector = file_dat_event_01{:,2}/1000000; %seconds
%Multiply the data by the coefficients
voltage_vector = file_dat_event_01{:,3}*v_multiplier+v_offset ;

%% Plots
figure;
plot(time_vector,voltage_vector);
xlabel('Time (s)');ylabel('Voltage (V)');title('Voltage over time');

%% System constants

h  = 1/sampling_rate ;

n = 400;
gm = 0.001; 
L1 = 26.526e-6;
L2 = 348e-6; 
L3 = 348e-6; 
M23 = 287e-6;
gs1 = 1.9635;
gs2 = 0.1497;
gs3 = 0.1497;
r1 = 0.005;
r2 = 0.4469;
r3 = 0.4469;
Rb = 0.1;
gb = 1/Rb;
h0 = 0.1876;
i0 = 6.09109;
L0 = 2.36;

% Define all X state vector
syms  v1t v2t v3t v4t et ht y1t y2t y3t y4t ipt imt iL1t iL2t iL3t
% Define all X state vector for the previous iteration
syms  v1h v2h v3h v4h eh hh y1h y2h y3h y4h iph imh iL1h iL2h iL3h

%% State estimation

% Rmatrix
R_vec = [0.005 0.005 0.005 0.005 0.0005 0.0005 0.0005 0.005 0.005 0.0005 0.00005 0.00005 0.00005 0.00005 0.0003 0.05 0.05 0.05 0.05 0.05 1];
R2_vec = R_vec.*R_vec;
R = repmat(R2_vec,21,1).*eye(21);
W = inv(R);

x_est_prev =  zeros(15,1);
x_est = zeros(size(voltage_vector,1),15);

% For Chi Square estimation
residuals_vec = zeros(1,21);
residuals = zeros(size(voltage_vector,1),1); 
results_chi2 = zeros(size(voltage_vector,1),1); 
confidence_level=0.1;
time_fault = 0;

residuals_mat = zeros(4700,21) ;

% for each data point of Vout 
for v=1:size(voltage_vector,1)
        disp(v);
        
        % Zmed
        Vout = voltage_vector(v);
        Ib =  -voltage_vector(v)/Rb; 
        
         if v==1
              ibh =  0; 
         else
             ibh = -voltage_vector(v-1)/Rb;
         end
        Zmed = [Vout 0 0 0 0 0 0 0 0 0 0 0 0 0 0 Ib Ib Ib Ib Ib 0]'; 
        
        % Initialize variables for first iteration 
        xk= zeros(15,1);
        
        v1h = x_est_prev(1);  v2h = x_est_prev(2); v3h = x_est_prev(3); v4h = x_est_prev(4); eh = x_est_prev(5); hh =x_est_prev(6); y1h  = x_est_prev(7); y2h = x_est_prev(8); 
        y3h = x_est_prev(9); y4h = x_est_prev(10); iph = x_est_prev(11);  imh = x_est_prev(12);  iL1h = x_est_prev(13);  iL2h = x_est_prev(14);  iL3h = x_est_prev(15); 
          
        Dx = ones(15,1)*10000; % just a very high value to ensure that the algorithm enters in the while loop
        % while loop until Dx is less than my tolerance
        while max(Dx)>0.0001
            
                v1t = xk(1);  v2t = xk(2); v3t = xk(3); v4t = xk(4); et = xk(5); ht = xk(6); y1t  = xk(7); y2t = xk(8); 
                y3t = xk(9); y4t = xk(10); ipt = xk(11);  imt = xk(12);  iL1t = xk(13);  iL2t = xk(14);  iL3t = xk(15); 
        
                % State vector 
                % X = [ v1(t)  v2(t)  v3(t)    v4(t)     e(t)                              h(t)                                                 y1(t)                    y2(t)                       y3(t)                     y4(t)                       ip(t)         im(t)           iL1(t)                iL2(t)                     iL3(t)]

                %Compute H matrix
%                 H =      [  0        0            1        -1          0                                 0                                                      0                            0                              0                           0                              0               0                  0                         0                            0  ;                              %1
%                                 0        0           0         0     -gm*h/2                          0                                                       0                            0                              0                           0                         h/2/n        -h/2     (h/2+gs1*L1)             0                            0  ;                                %2
%                                 0        0           0         0      gm*h/2                           0                                                      0                             0                              0                           0                        -h/2/n        h/2                 0               (-h/2-gs2*L2)      gs2*M23  ;                       %3
%                                 0        0           0         0     -gm*h/2                          0                                                       0                            0                              0                           0                          h/2/n       -h/2                0                -gs3*M23        (h/2+gs3*L3) ;                      %4
%                               -h/2   h/2        0         0   (h/2+r1*gm*h/2)           0                                                      0                             0                              0                           0                     -r1*h/2/n    r1*h/2           L1                      0                            0   ;                                 %5
%                                h/2     0        -h/2      0            0                                 0                                                      0                             0                              0                            0                             0              0                   0       (L2+r2*(h/2+gs2*L2))   (-M23-r2*gs2*M23) ;      %6
%                                0      -h/2         0      h/2          0                                 0                                                       0                            0                              0                            0                             0              0                   0        -M23*(1+r3*gs3)       (L3*(1+r3*gs3)+r3*h/2) ;      %7
%                                0        0     gb*h/2  -gb*h/2   0                                0                                                      0                              0                             0                            0                             0              0                   0               (h/2+gs2*L2)         -gs2*M23    ;                       %8
%                                0        0   -gb*h/2   gb*h/2   0                                 0                                                      0                              0                             0                           0                              0              0                    0                 gs3*M23            (-h/2-gs3*L3);                       %9
%                                0        0            0        0            h/2                             -1                                                     0                             0                              0                            0                             0              0                    0                        0                              0 ;                                %10
%                                0        0            0        0             0           -1/(h0^2)*(hh+2*ht)*h/3                             h/2                          0                              0                            0                              0              0                   0                        0                               0 ;                                %11
%                                0        0            0        0             0                                 0                                           -(y1h+2*y1t)*h/3          h/2                         0                            0                              0             0                    0                         0                              0 ;                          %12
%                                0        0            0        0             0                                 0                                                       0              -(y2h+2*y2t)*h/3           h/2                          0                             0              0                    0                        0                              0  ;                          %13
%                                0        0            0        0             0                                 0                                            -(y3h/2+y3t)*h/3        0              -(y1h/2+y1t)*h/3            h/2                           0              0                    0                        0                              0  ;                          %14
%                                0        0            0        0             0    (-h/2/L0-i0/h0*(y4h+y4t/2)*h/3)                       0                            0                             0        -i0/h0*(hh/2+ht)*h/3           0             h/2                 0                        0                              0 ;                          %15
%                                0        0   -gb*h/2   gb*h/2   0                                  0                                                      0                            0                              0                           0                              0              0                    0                         0                              0 ;                       %16
%                                0        0           0        0             0                                    0                                                     0                            0                             0                           0                               0               0         (h/2+gs1*L1)             0                              0  ;                     %17
%                                0        0           0        0             0                                    0                                                     0                            0                             0                           0                               0               0                   0              (h/2+gs2*L2)            -gs2*M23 ;                     %18
%                                0        0           0        0             0                                    0                                                     0                             0                             0                           0                              0               0                   0                  gs3*M23             (-h/2-gs3*L3);                     %19
%                                0        0           0        0         gm*h/2                           0                                                     0                             0                             0                           0                         -h/2/n        h/2                  0                        0                               0 ;                     %19
%                                0        0           0      h/2          0                                    0                                                     0                             0                             0                           0                               0               0                   0                         0                               0 ;
%                                ];
                           
                H =      [  0        0            1        -1          0                                 0                                                      0                            0                              0                           0                              0               0                  0                         0                            0  ;                              %1
                                0        0           0         0     -gm*h/2                          0                                                       0                            0                              0                           0                         h/2/n        -h/2     (h/2+gs1*L1)             0                            0  ;                                %2
                                0        0           0         0      gm*h/2                           0                                                      0                             0                              0                           0                        -h/2/n        h/2                 0               (-h/2-gs2*L2)      gs2*M23  ;                       %3
                                0        0           0         0     -gm*h/2                          0                                                       0                            0                              0                           0                          h/2/n       -h/2                0                -gs3*M23        (h/2+gs3*L3) ;                      %4
                              -h/2   h/2        0         0   (h/2+r1*gm*h/2)           0                                                      0                             0                              0                           0                     -r1*h/2/n    r1*h/2           L1                      0                            0   ;                                 %5
                               h/2     0        -h/2      0            0                                 0                                                      0                             0                              0                            0                             0              0                   0       (L2+r2*(h/2+gs2*L2))   (-M23-r2*gs2*M23) ;      %6
                               0      -h/2         0      h/2          0                                 0                                                       0                            0                              0                            0                             0              0                   0        -M23*(1+r3*gs3)       (L3*(1+r3*gs3)+r3*h/2) ;      %7
                               0        0     gb*h/2  -gb*h/2   0                                0                                                      0                              0                             0                            0                             0              0                   0               (h/2+gs2*L2)         -gs2*M23    ;                       %8
                               0        0   -gb*h/2   gb*h/2   0                                 0                                                      0                              0                             0                           0                              0              0                    0                 gs3*M23            (-h/2-gs3*L3);                       %9
                               0        0            0        0            h/2                             -1                                                     0                             0                              0                            0                             0              0                    0                        0                              0 ;                                %10
                               0        0            0        0             0                 -1/(h0^2)*(2*ht)                                        1                              0                             0                            0                              0              0                   0                        0                               0 ;                                %11
                               0        0            0        0             0                                 0                                                -(2*y1t)                      1                             0                            0                              0               0                    0                         0                              0 ;                          %12
                               0        0            0        0             0                                 0                                                       0                     -(2*y2t)                       1                             0                             0              0                    0                        0                              0  ;                          %13
                               0        0            0        0             0                                 0                                                   -(y3t)                        0                        -(y1t)                         1                             0              0                    0                        0                              0  ;                          %14
                               0        0            0        0             0                  -(i0/h0)*y4t-1/L0                                       0                            0                             0                 -(i0/h0)*ht                      0              1                    0                        0                              0 ;                          %15
                               0        0         -gb     gb             0                                  0                                                      0                            0                              0                           0                              0              0                    0                         0                              0 ;                       %16
                               0        0           0        0             0                                    0                                                     0                            0                             0                           0                               0               0    (h/2+gs1*L1)*(2/h)     0                              0  ;                     %17
                               0        0           0        0             0                                    0                                                     0                            0                             0                           0                               0               0                   0     (h/2+gs2*L2)*(2/h)   -gs2*M23*(2/h)    ;                     %18
                               0        0           0        0             0                                    0                                                     0                             0                             0                           0                              0               0                   0         gs3*M23*(2/h)        (-h/2-gs3*L3)*(2/h)   ;                     %19
                               0        0           0        0           gm                                 0                                                     0                             0                             0                           0                           -1/n            1                  0                        0                               0 ;                     %20
                               0        0           0        1          0                                    0                                                     0                             0                             0                           0                               0               0                   0                         0                               0 ;
                               ];
                    
                           %Compute f(X)
%                     f = [ v3t-v4t;                                                                                                                                                                          %1
%                             -gm*h/2*(et+eh)-h/2*(imt+imh)+h/2/n*(ipt+iph)+h/2*(iL1t+iL1h)+gs1*L1*(iL1t-iL1h);                                                                   %2
%                             gm*h/2*(et+eh)+h/2*(imt+imh)-h/2/n*(ipt+iph)-h/2*(iL2t+iL2h)-gs2*L2*(iL2t-iL2h)+gs2*M23*(iL3t-iL3h);                            %3  
%                             -gm*h/2*(et+eh)-h/2*(imt+imh)+h/2/n*(ipt+iph)+h/2*(iL3t+iL3h)+gs3*L3*(iL3t-iL3h)-gs3*M23*(iL2t-iL2h);                               %4
%                             -h/2*(v1t+v1h)+h/2*(v2t+v2h)+h/2*(et+eh)+L1*(iL1t-iL1h)+r1*gm*h/2*(et+eh)+r1*h/2*(imt+imh)-r1*h/2/n*(ipt+iph);                  %5  
%                             -h/2*(v3t+v3h)+h/2*(v1t+v1h)+r2*h/2*(iL2t+iL2h)+r2*gs2*L2*(iL2t-iL2h)-r2*gs2*M23*(iL3t-iL3h)+L2*(iL2t-iL2h)-M23*(iL3t-iL3h);        %6
%                             -h/2*(v2t+v2h)+h/2*(v4t+v4h)+r3*h/2*(iL3t+iL3h)+r3*gs3*L3*(iL3t-iL3h)-r3*gs3*M23*(iL2t-iL2h)+L3*(iL3t-iL3h)-M23*(iL2t-iL2h);            %7
%                             h/2*(iL2t+iL2h)+gs2*L2*(iL2t-iL2h)-gs2*M23*(iL3t-iL3h)+gb*h/2*(v3t+v3h)-gb*h/2*(v4t+v4h);                                                                           %8
%                             -h/2*(iL3t+iL3h)-gs3*L3*(iL3t-iL3h)+gs3*M23*(iL2t-iL2h)+gb*h/2*(v4t+v4h)-gb*h/2*(v3t+v3h);                                                              %9
%                             h/2*(et+eh)-(ht-hh);                                                                                                                                                                                                            %10
%                             h/2*(y1t+y1h)-1/(h0^2)*(hh^2*h/3+hh*ht*h/6+ht*hh*h/6+ht^2*h/3);                                     %11
%                             h/2*(y2t+y2h)-(y1h^2*h/3+y1h*y1t*h/6+y1t*y1h*h/6+y1t^2*h/3);                                            %12
%                             h/2*(y3t+y3h)-(y2h^2*h/3+y2h*y2t*h/6+y2t*y2h*h/6+y2t^2*h/3);                                            %13
%                             h/2*(y4t+y4h)-(y3h*y1h*h/3+y3h*y1t*h/6+y3t*y1h*h/6+y3t*y1t*h/3);                                    %14
%                             h/2*(imt+imh)-i0/h0*(hh*y4h*h/3+hh*y4t*h/6+ht*y4h*h/6+ht*y4t*h/3)-h/2/L0*(ht+hh);           %15
%                             -gb*(v3t - v4t);                                                           %16
%                             (h/2*(iL1t+iL1h)+gs1*L1*(iL1t-iL1h))*(2/h)-ibh;                                                                     %17
%                             (h/2*(iL2t+iL2h)+gs2*L2*(iL2t-iL2h)-gs2*M23*(iL3t-iL3h))*(2/h)-ibh;                                 %18
%                             (-h/2*(iL3t+iL3h)-gs3*L3*(iL3t-iL3h)+gs3*M23*(iL2t-iL2h))*(2/h)-ibh;                    %19
%                             gm*et+imt-1/n*ipt;                               %20
%                             h/2*(v4t+v4h)                                                                                           %21
%                     ];
                
                     f = [ v3t-v4t;                                                                                                                                                                          %1
                            -gm*h/2*(et+eh)-h/2*(imt+imh)+h/2/n*(ipt+iph)+h/2*(iL1t+iL1h)+gs1*L1*(iL1t-iL1h);                                                                   %2
                            gm*h/2*(et+eh)+h/2*(imt+imh)-h/2/n*(ipt+iph)-h/2*(iL2t+iL2h)-gs2*L2*(iL2t-iL2h)+gs2*M23*(iL3t-iL3h);                            %3  
                            -gm*h/2*(et+eh)-h/2*(imt+imh)+h/2/n*(ipt+iph)+h/2*(iL3t+iL3h)+gs3*L3*(iL3t-iL3h)-gs3*M23*(iL2t-iL2h);                               %4
                            -h/2*(v1t+v1h)+h/2*(v2t+v2h)+h/2*(et+eh)+L1*(iL1t-iL1h)+r1*gm*h/2*(et+eh)+r1*h/2*(imt+imh)-r1*h/2/n*(ipt+iph);                  %5  
                            -h/2*(v3t+v3h)+h/2*(v1t+v1h)+r2*h/2*(iL2t+iL2h)+r2*gs2*L2*(iL2t-iL2h)-r2*gs2*M23*(iL3t-iL3h)+L2*(iL2t-iL2h)-M23*(iL3t-iL3h);        %6
                            -h/2*(v2t+v2h)+h/2*(v4t+v4h)+r3*h/2*(iL3t+iL3h)+r3*gs3*L3*(iL3t-iL3h)-r3*gs3*M23*(iL2t-iL2h)+L3*(iL3t-iL3h)-M23*(iL2t-iL2h);            %7
                            h/2*(iL2t+iL2h)+gs2*L2*(iL2t-iL2h)-gs2*M23*(iL3t-iL3h)+gb*h/2*(v3t+v3h)-gb*h/2*(v4t+v4h);                                                                           %8
                            -h/2*(iL3t+iL3h)-gs3*L3*(iL3t-iL3h)+gs3*M23*(iL2t-iL2h)+gb*h/2*(v4t+v4h)-gb*h/2*(v3t+v3h);                                                              %9
                            h/2*(et+eh)-(ht-hh);                                                                                                                                                                                                            %10
                            y1t-(1/(h0^2))*(ht^2);                                     %11
                            y2t-(y1t^2);                                            %12
                            y3t-(y2t^2);                                            %13
                            y4t-(y3t*y1t);                                    %14
                            imt-(i0/h0)*ht*y4t-(1/L0)*(ht);           %15
                            -gb*(v3t-v4t);                                                           %16
                            (h/2*(iL1t+iL1h)+gs1*L1*(iL1t-iL1h))*(2/h)-ibh;                                                                     %17
                            (h/2*(iL2t+iL2h)+gs2*L2*(iL2t-iL2h)-gs2*M23*(iL3t-iL3h))*(2/h)-ibh;                                 %18
                            (-h/2*(iL3t+iL3h)-gs3*L3*(iL3t-iL3h)+gs3*M23*(iL2t-iL2h))*(2/h)-ibh;                    %19
                            gm*et+imt-1/n*ipt;                               %20
                            v4t                                                               %21
                    ];

                    % Compute delta X (Dx)
                    Dx  = inv(H'*W*H)*H'*W*(Zmed-f);
                    
                    % Compute estimation
                    xk1 = xk+Dx; 
                    
                    % Update vectors 
                    xk = xk1; 
                    
                    % Chi square test
                    
                    residuals_vec = (Zmed - f).^2; 
                    residuals_vec = residuals_vec';
                    
                    residuals_mat(v,:) = residuals_vec./R2_vec;
                    residuals(v) = sum(residuals_vec./R2_vec); 
                    results_chi2(v) = 1-chi2cdf(residuals(v), 6);
                    results_chi2(v)
                    if results_chi2(v) < (1-confidence_level);
                        time_fault = time_vector(v);
                        disp(time_fault);  
                        disp("HOla"); 
                   end
                    
        end

x_est(v,:) = xk1;
x_est_prev = x_est(v,:);

end

% %% Chi Square estimation
% residuals_vec = zeros(1,21);
% residuals = zeros(size(voltage_vector,1),1); 
% results_chi2 = zeros(size(voltage_vector,1),1); 
% 
% % for each data point of Vout 
% for v=1:size(voltage_vector,1)
%     % Zmed
%     Vout = voltage_vector(v);
%     Ib =  -voltage_vector(v)/Rb; 
%         
%     
%     xk = x_est(v,:); 
%     
%     if v == 1
%         x_est_prev = zeros(1,15);
%         ibh = 0; 
%     else
%         x_est_prev = x_est(v-1,:); 
%         ibh = (x_est(v-1,3)-x_est(v-1,4))/Rb;   %estimated Vout (t - h)/Rb
%     end
%     
%     Zmed = [Vout 0 0 0 0 0 0 0 0 0 0 0 0 0 0 Ib Ib Ib Ib Ib 0]'; 
%     
%     v1h = x_est_prev(1);  v2h = x_est_prev(2); v3h = x_est_prev(3); v4h = x_est_prev(4); eh = x_est_prev(5); hh =x_est_prev(6); y1h  = x_est_prev(7); y2h = x_est_prev(8); 
%     y3h = x_est_prev(9); y4h = x_est_prev(10); iph = x_est_prev(11);  imh = x_est_prev(12);  iL1h = x_est_prev(13);  iL2h = x_est_prev(14);  iL3h = x_est_prev(15); 
% 
%     v1t = xk(1);  v2t = xk(2); v3t = xk(3); v4t = xk(4); et = xk(5); ht = xk(6); y1t  = xk(7); y2t = xk(8); 
%     y3t = xk(9); y4t = xk(10); ipt = xk(11);  imt = xk(12);  iL1t = xk(13);  iL2t = xk(14);  iL3t = xk(15); 
% 
%    f = [ v3t-v4t;                                                                                                                                                                          %1
%                             -gm*h/2*(et+eh)-h/2*(imt+imh)+h/2/n*(ipt+iph)+h/2*(iL1t+iL1h)+gs1*L1*(iL1t-iL1h);                                                                   %2
%                             gm*h/2*(et+eh)+h/2*(imt+imh)-h/2/n*(ipt+iph)-h/2*(iL2t+iL2h)-gs2*L2*(iL2t-iL2h)+gs2*M23*(iL3t-iL3h);                            %3  
%                             -gm*h/2*(et+eh)-h/2*(imt+imh)+h/2/n*(ipt+iph)+h/2*(iL3t+iL3h)+gs3*L3*(iL3t-iL3h)-gs3*M23*(iL2t-iL2h);                               %4
%                             -h/2*(v1t+v1h)+h/2*(v2t+v2h)+h/2*(et+eh)+L1*(iL1t-iL1h)+r1*gm*h/2*(et+eh)+r1*h/2*(imt+imh)-r1*h/2/n*(ipt+iph);                  %5  
%                             -h/2*(v3t+v3h)+h/2*(v1t+v1h)+r2*h/2*(iL2t+iL2h)+r2*gs2*L2*(iL2t-iL2h)-r2*gs2*M23*(iL3t-iL3h)+L2*(iL2t-iL2h)-M23*(iL3t-iL3h);        %6
%                             -h/2*(v2t+v2h)+h/2*(v4t+v4h)+r3*h/2*(iL3t+iL3h)+r3*gs3*L3*(iL3t-iL3h)-r3*gs3*M23*(iL2t-iL2h)+L3*(iL3t-iL3h)-M23*(iL2t-iL2h);            %7
%                             h/2*(iL2t+iL2h)+gs2*L2*(iL2t-iL2h)-gs2*M23*(iL3t-iL3h)+gb*h/2*(v3t+v3h)-gb*h/2*(v4t+v4h);                                                                           %8
%                             -h/2*(iL3t+iL3h)-gs3*L3*(iL3t-iL3h)+gs3*M23*(iL2t-iL2h)+gb*h/2*(v4t+v4h)-gb*h/2*(v3t+v3h);                                                              %9
%                             h/2*(et+eh)-(ht-hh);                                                                                                                                                                                                            %10
%                             y1t-(1/(h0^2))*(ht^2);                                     %11
%                             y2t-(y1t^2);                                            %12
%                             y3t-(y2t^2);                                            %13
%                             y4t-(y3t*y1t);                                    %14
%                             imt-(i0/h0)*ht*y4t-(1/L0)*(ht);           %15
%                             -gb*(v3t-v4t);                                                           %16
%                             (h/2*(iL1t+iL1h)+gs1*L1*(iL1t-iL1h))*(2/h)-ibh;                                                                     %17
%                             (h/2*(iL2t+iL2h)+gs2*L2*(iL2t-iL2h)-gs2*M23*(iL3t-iL3h))*(2/h)-ibh;                                 %18
%                             (-h/2*(iL3t+iL3h)-gs3*L3*(iL3t-iL3h)+gs3*M23*(iL2t-iL2h))*(2/h)-ibh;                    %19
%                             gm*et+imt-1/n*ipt;                               %20
%                             v4t                                                               %21
%                     ];
%     
%     
%     residuals_vec = (Zmed - f).^2; 
%     sprintf('Zmed 17: %f', Zmed(17))
%     sprintf('f 17: %f', f(17))
%     sprintf('Residuals vector')
%     disp(residuals_vec)
%     residuals_vec = residuals_vec';
%     residuals(v) = sum(residuals_vec./R2_vec); 
%     sprintf('Residuals vector/R2 Iteration: %d', v)
%     disp(residuals_vec./R2_vec)
%     sprintf('Sum of Residuals/R2')
%     disp(residuals(v))
%     results_chi2(v) = chi2cdf(residuals(v), 6);
%     
%     
% end 

%% Results

figure; 
plot(time_vector,x_est(:,3)-x_est(:,4));
xlabel('Time (s)');ylabel('Estimated Voltage (V)');title('State estimation of Vout');

figure; 
plot(time_vector,x_est(:,11));
xlabel('Time (s)');ylabel('Estimated Current (A)');title('State estimation of Ip');

figure; 
plot(time_vector,x_est(:,11));
hold on 
plot(time_vector, n*voltage_vector/Rb);
xlabel('Time (s)');ylabel('Current (A)');title('CT*Ib and estimated Ip');
legend('Estimated Ip','CT*Ib');

figure; 
plot(time_vector,x_est(:,11)-n*voltage_vector/Rb);
xlabel('Time (s)');ylabel('Current (A)');title('Difference between CT*Ib and estimated Ip');

figure; 
plot(time_vector,residuals);
xlabel('Time (s)');ylabel('Residuals');title('State estimation residuals');

figure; 
plot(time_vector,results_chi2);
xlabel('Time (s)');ylabel('Residuals');title('Chi Square test');

figure; 
plot(time_vector,x_est(:,6));
xlabel('Time (s)');ylabel('Flux linkage');title('Flux linkage over time');


%%

figure; 
plot(residuals_mat(:,1:21));
legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10',  '11',  '12',  '13' , '14', '15', '16','17','18','19','20','21')

%legend('v1(t)', 'v2(t)', 'v3(t)', 'v4(t)', 'e(t)', 'h(t)', 'y1(t)', 'y2(t)', 'y3(t)', 'y4(t)',  'ip(t)',  'im(t)',  'iL1(t)' , 'iL2(t)', 'iL3(t)')
xlabel('Time (s)');ylabel('Residuals');title('Residuals');


figure; 
plot(residuals_mat(:,4));

figure; 
plot(residuals_mat(:,9));
