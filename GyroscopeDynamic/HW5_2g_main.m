clear;

% IC
format long;
THETA_0 = [0 15*pi/180]';

Csb_0 = eye(3,3);   % Initial DCM
w_ba_s_0  = [-0.002588190451025  0  0.009659258262891]';
w_ba_b_0 = Csb_0'*w_ba_s_0;

X_0 = [THETA_0;w_ba_s_0];

% time span
t = [0 2];

% ODE45
options = odeset('RelTol',1e-6);
[T,Y] = ode45(@HW5_2g_fun, t, X_0,options);

% ------------------------------------plot---------------------------------
figure(1);
subplot(3,2,1)
plot(T,rad2deg(Y(:,1))); 
xlabel('Time [s]');ylabel('\theta [{\circ}]','Interpreter','Tex');
title('\theta versus time','Interpreter','Tex')
subplot(3,2,2)
plot(T,rad2deg(Y(:,2))); 
xlabel('Time [s]');ylabel('\phi [{\circ}]','Interpreter','Tex');
title('\phi versus time','Interpreter','Tex')
subplot(3,2,3)
plot(T,Y(:,3));
xlabel('Time [s]');ylabel('\omega^{ba}_{b1} [{rad/s}]','Interpreter','Tex');
title('\omega^{ba}_{b1} versus time','Interpreter','Tex')
subplot(3,2,4)
plot(T,Y(:,4));
xlabel('Time [s]');ylabel('\omega^{ba}_{b2} [{rad/s}]','Interpreter','Tex');
title('\omega^{ba}_{b2} versus time','Interpreter','Tex')
subplot(3,2,5)
plot(T,Y(:,5));
xlabel('Time [s]');ylabel('\omega^{ba}_{b3} [{rad/s}]','Interpreter','Tex');
title('\omega^{ba}_{b3} versus time','Interpreter','Tex')

% --------------------Compute the total  energy---------------------------

SIZE = size(Y);
 E = zeros(SIZE(1),1);
dE = zeros(SIZE(1),1);

for i=1:(SIZE(1))
    [dx E_tmp] = HW5_2g_fun(T(i),Y(i,:));
    %caculate dE versus time
     E(i) = E_tmp;
     if i ==1
    dE(i) = 0;
     else
    dE(i) = (E(i)-E(i-1));
     end
end
        
% plot energy stuffs
figure(2);
plot(T,E)
xlabel('Time [s]');ylabel('E [J]');title('total energy versus time'); 
ylim([.45833 .45834])
figure(3)
plot(T,dE);
xlabel('Time [s]');ylabel('\Delta E_{Bw/a} [J]','Interpreter','Tex');
title('Change of energy versus time');

save AE_540_HW5.mat;