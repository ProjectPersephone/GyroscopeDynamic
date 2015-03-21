function [dx E] = HW5_2g2 (t,x)
% Function for the ODE45 and also used for computing total energy.
% Wei Ding
% 11/8/2014

    % get in constants
    constants;
    eta = deta*t + eta_0;
    
       theta = x(1);
       phi   = x(2);
       
    THETA_ba = [x(1) x(2)]';

    % --------------------------Kinematics---------------------------------
    %       1.DCM
    Cqa = [cos(theta) sin(theta)    0;
          -sin(theta) cos(theta)    0
                   0           0    1];

    Cbq = [cos(phi)     0   -sin(phi)
                  0     1           0
           sin(phi)     0    cos(phi)];

    Csb = [cos(eta)     sin(eta)    0;
          -sin(eta)     cos(eta)    0
                  0     0           1];

    Csa = Csb*Cbq*Cqa;
    %       2.Angular Velocity
    w_ba_b = [x(3) x(4) x(5)]';
    w_ba_s = Csb*w_ba_b;
    w_sb = [0 0 deta]';
    w_sa = w_sb + w_ba_s;

    %       S_ba
    S_ba = Csb*[Cbq*[0 0 1]' [0 1 0]'];

    %       3.Positon
    r_cw_s = [0 0 l]';
    r_cw_b = r_cw_s;

    %       4.Velocity
    %       5.Acceleration

    % --------------------------------FBD---------------------------------
    %       1.Forces
    %       2.Moments

    m_Bw_b = [0 sin(phi)*l*m*g 0]';
    m_Bw_s = Csb*m_Bw_b;

    % ---------------------------------E2L-------------------------------
    %       1. Momentum
    %       2. Angular Momentum
    %       3. E2L
    % ------------------->Get the DEs

    % ----------------------------------DE-------------------------------
    %          Inertia
    J_MbW_b   = -m*CROSS(r_cw_b)*CROSS(r_cw_b);
    J_MbW_s   = Csb*J_MbW_b*Csb';
    J_BC_s    = [1/12*m*(3*(D/2)^2+h^2)                        0                  0;
                                    0   1/12*m*(3*(D/2)^2+h^2)                  0;
                                    0                        0      1/2*m*(D/2)^2];


    %           DE of w_ba
    dw_ba_s = inv(J_BC_s+J_MbW_s)*(m_Bw_s - CROSS(w_sa)*J_BC_s*w_sa - CROSS(w_ba_s)*J_MbW_s*w_ba_s - J_MbW_s*CROSS(w_sb)*w_ba_s);
    dw_ba_b = Csb'*(dw_ba_s+CROSS(w_sb)*w_ba_s); % From s frame to b frame
    %           DE of theta_ba
    dTHETA_ba = S_ba'*w_ba_s;
    %           spits out time derivatives
    dx = [dTHETA_ba;dw_ba_b]; 
    
    %---------------------------Total Energy------------------------------
    % Kinetic
    T_Bw = 1/2*w_sa'*J_BC_s*w_sa+1/2*w_ba_b'*J_MbW_b*w_ba_b;
    % Potential
    g_a  = [0 0 -g]';
    U_Bw = -r_cw_s'*Csa*g_a*m;
    % Total
    E = T_Bw+U_Bw;

end

function X_cross = CROSS( X )
% Cross of a Matrix

X_cross = zeros(3,3);
X_cross = [  0    -X(3)   X(2);
           X(3)     0    -X(1);
           -X(2)   X(1)   0
];

end
