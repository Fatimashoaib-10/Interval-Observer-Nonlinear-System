%% Nonlinear discrete time system
close all;clc;clear all

n = 2;
flagH = 1;

A = [0 1;0.3 0]; 

r = 0.05*[1;0];

amin = 1;
amax = 1.4;
a = 1; %mean([amin amax]);
d = 0;%.015/2;

C  = [1 0];
p = size(C,1);


%% discretized system
func = @(x) r*(1-a*x(1)^2 + d); % A*x + r*(1-a*x(1)^2 + d);
Outfunc = @(x) C*x;

interval_X = interval([-2;-1],[2;1]);
X_upper = interval_X.sup;
X_lower = interval_X.inf;

%% Jacobian
% Gradient
x = sym('x',[n 1]);
grad = jacobian(func(x),x);
gradfunc = matlabFunction(grad,'Vars',{x});   
Jf= gradfunc(interval_X);
Jf_upp = Jf.sup;
Jf_low = Jf.inf;

%grad_y = jacobian(Outfunc(x),x);
%gradfunc_y = matlabFunction(grad_y,'Vars',{x}); 
%Jy = gradfunc_y(interval_X);
Jy_upp = [1 0];%Jy.sup;
Jy_low = [1 0];%Jy.inf;

switch flagH
    case 1
        H=zeros(n,n);
        for i=1:n
            for j=1:n
                if Jf_low(i,j)>=0 || Jf_upp(i,j)<=0
                    H(i,j)=0;
                elseif abs(Jf_upp(i,j))>abs(Jf_low(i,j))
                    H(i,j)=Jf_low(i,j);
                else
                    H(i,j)=Jf_upp(i,j);
                end
            end
        end
    case 2
        JJ       = cell(1,n);
        fupp     = cell(1,n);
        floww    = cell(1,n);
        HH       = cell(1,n);
        H = [];
        for i = 1:n
            JJ{i} = Jf(i,:);
            ind   = cell(1,n);
            
            [~,~,fupp{i},floww{i},HH{i}] = decomp_signstable_modified(func,...
                interval_X.sup,interval_X.inf,JJ{i}.sup,JJ{i}.inf);
            [~,ind{i}] = min(fupp{i}-floww{i});
            H = [H;HH{i}(ind{i},:)]; % Minimum H case
        end
end

J_phi=Jf-H;
Jphi_upp = J_phi.sup;
Jphi_low = J_phi.inf;

A_old=A;
A=A+H;

%%
F_upp_phi = 2*max(Jphi_upp, zeros(n,n)) - Jphi_low ;
F_upp_ksi = 2*max(Jy_upp-C, zeros(p,n)) - Jy_low + C;

%%
P = sdpvar(n);
X = sdpvar(n);
%X=diag(Xd);
J = sdpvar(1,n);
% alpha = 0.5;

M11 = -P;
M12 = ((abs(A))'+ F_upp_phi')*X + (C' + F_upp_ksi')*J;
M22 = P - X - X';

MAT = [ M11 M12 ; M12' M22];
F = [MAT <= -1e-3*eye(2*n)];

[m,n]=size(X);
constraints = [];
for i = 1:m
    for j =1:n
        if( i ~= j)
            constraints = [constraints,(X(i,j)) <= 0  ];
        end
    end
end

F = [F; J'*C >= zeros(n,n) ; constraints; J >=zeros(1,n) ];
%options = sdpsettings('verbose',0,'solver','MOSEK');
diagnostics = optimize(F)
P = value(P)
X = value(X)
J = value(J)
L = inv((X)')*J'
J'*C