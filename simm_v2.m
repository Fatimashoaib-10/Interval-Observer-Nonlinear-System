%% simulating the system


function [t,x] = simm_v2(L, H ,A,Jphi_upp)
% x0 = [-0.5;-0.5;-0.5;-0.5;-0.5;-0.5;-0.5;-0.5;-0.5];

% x0 = [10;10;10;11;11;11;9;9;9];

x0 = [10;5;2;11;6;3;9;4;1];

tspan = 0:0.1:5;%linspace(0,0.5,50);

[t,x] = ode45(@dots, tspan, x0); % plotting the solution to linear system

figure(1);hold on;grid on;%box on;
%line_color = ['b' 'g' 'y' 'b-' 'g-' 'y-','b*-','g*-','y*-'];
plot(t,x,'.-');shg
xlabel('time');ylabel('x')
legend('x(1)','x(2)','x(3)',...
    'x(1)_{upp}','x(2)_{upp}','x(3)_{upp}',...
    'x(1)_{low}','x(2)_{low}','x(3)_{low}','Location','Best')

figure(2);hold on;
plot(t,x(:,1),'--','color','blue');
hold on;
plot(t,x(:,4),':','color','black');
hold on;
plot(t,x(:,7),':','color','red');
legend('x(1)','x(1)_{upp}','x(1)_{low}');

figure(3);hold on;
plot(t,x(:,2),'--','color','blue');
hold on;
plot(t,x(:,5),':','color','black');
hold on;
plot(t,x(:,8),':','color','red');
legend('x(2)','x(2)_{upp}','x(2)_{low}');

figure(4);hold on;
plot(t,x(:,3),'--','color','blue');
hold on;
plot(t,x(:,6),':','color','black');
hold on;
plot(t,x(:,9),':','color','red');
legend('x(3)','x(3)_{upp}','x(3)_{low}');


%%
function dxs = dots(t,X)

%load data.mat L H A J_phi Jphi_upp Jphi_low

n = 3;
x    = X(1:3);
xupp = X(4:6);
xlow = X(7:9);

A = [2 0 0; 1 -4 sqrt(3); -1 -sqrt(3) -4];
B = 0.15*[-2*5.3; 0; 3.4];
C = [1 0 0];
A_new = A+H;

y  = C*x(:);
u  = 1+sin(2*t);

dx     = A*x(:) + B*x(1)*x(2)*u;
A_LC_d  = diag(diag(A_new-L*C));
A_LC_nd = (A_new-L*C)-A_LC_d;

A_LC_ndpos = max(A_LC_nd,0);
A_LC_ndneg = A_LC_ndpos - A_LC_nd;

% Lpos = max(L,0);
% Lneg = Lpos - L;

%phi = @(x,u) B*x(1)*x(2)*u;
phi = @(x) B*x(1)*x(2)*u - H*x;

%%
% D = zeros(n);
% D_upp = max(sign(Jphi_low),zeros(n));
% D_low = ones(n) - D_upp;

D = max(sign(Jphi_upp),0);

% phi = @(x,y)(D_upp*x + (eye(n) - D_upp)*y);

% phi_d = @(x,y)(D*x + (eye(n) - D)*y);


% D = sign (Jphi_low) + ones(n);
% D_upp = D;
% D_low = -sign(Jphi_low);
% for i = 1:n
%     for j = 1:n
%         if i == j
%             D_upp(i,j)= 1;
%             D_low(i,j)= 0;
%         end
%     end
% end
% phi_upper = phi(D_upp * xupp + (eye(n) - D_upp)* xlow);
% phi_lower = phi((eye(n) - D_low) * xupp + D_low * xlow);


dx_upp = A_LC_d*xupp(:) + A_LC_ndpos*xupp(:) - A_LC_ndneg*xlow(:) + L*y +...
    phi_d_gen(xupp,xlow,D,phi,n); %phi(xupp,xlow);
dx_low = A_LC_d*xlow(:) + A_LC_ndpos*xlow(:) - A_LC_ndneg*xupp(:) + L*y +...
    phi_d_gen(xlow,xupp,D,phi,n); %phi(xlow,xupp);

dxs = [dx;dx_upp;dx_low];

end

function phi_d = phi_d_gen(x,y,D,phi,n)
    
    phi_d=zeros(n,1);
    for i=1:n
        Di = diag(D(i,:));
        phi_di=phi(Di*x+(eye(n)-Di)*y);
        phi_d(i)=phi_di(i);
    end
end
end
