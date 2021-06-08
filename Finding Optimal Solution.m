%Assignment 3
%Ryan Klughart 36627875
%----------------1----------------------
x = sdpvar(1,1);
y = sdpvar(1,1);

ops1 = sdpsettings('solver','quadprog','debug',1);
Objective1 = (x+2*y-7)^2 + (2*x + y-5)^2;
Constraints1 = [x>=1,y>=1];
fval1 = optimize(Constraints1,Objective1,ops1);


x1 = value(x);

y1 = value(y);

%----------------1990P1----------------------



x = sdpvar(5,1);
c = [42 44 45 47 47.5]';
Q = 100*eye(5);

% Define constraints 
Constraints2 = [20*x(1) + 12*x(2) + 11*x(3) + 7*x(4) + 4*x(5) <= 40, (0 <= x) <= 1];

% Define an objective
Objective2 = c'*x - 0.5*x'*Q*x;

% Set some options for YALMIP and solver
options = sdpsettings('solver','quadprog');

% Solve the problem
sol = optimize(Constraints2,Objective2, options);
x2 = value(x);
x2
%---------------------------4---------------------

x = sdpvar(1,1);
y = sdpvar(1,1);

Constraints4  = [x>=0,y>=0];
Objective4 = x.^2 +y.^2 +2*x*y;
fval4 = optimize(Constraints4,Objective4,options);


%---------------------------5---------------------

x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);
Constraints5 = [x>=1,y>=-5];
Objective5 = x.^2 +y;
fval5 = optimize(Constraints5,Objective5,options);


%---------------------------6---------------------

x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);
Constraints6 = [x>=z,y>=x];
Objective6 = x.^2 +(y+z).^2 +5*y;
fval6 = optimize(Constraints6,Objective6,options);



%------------------sphere------------------------
x = sdpvar(1,2);
ObjectiveSp = sum(x.^2,2);
ConstraintsSp = x>=-1;

fvalSp= optimize(ConstraintsSp,ObjectiveSp,options);
fvalSp
x = value(x);
x

%------------------booth------------------------

x = sdpvar(1,1);
y = sdpvar(1,1);
ObjectiveBo = (x+2*y-7).^2 + (2*x + y -5).^2;
ConstraintsBo = x>=10;

fvalBo= optimize(ConstraintsBo,ObjectiveBo,options);

x = value(x);

%---------------1990 p2-------------------------

%1990 Christodoulos A. Floudas Panos M. Pardalos  test problem 2
% Define variables
x = sdpvar(5,1);
y = sdpvar(1,1);
c = [-10.5 -7.5 -3.5 -2.5 -1.5]';
Q = eye(5);
d = -10;

% Define constraints 
ConstraintsP2 = [6*x(1) + 3*x(2) + 3*x(3) + 2*x(4) + x(5) <= 6.5, 10*x(1) + 10*x(3) + y <= 20, (0 <= x) <= 1, y >= 0];

% Define an objective
ObjectiveP2 = c.'*x - 0.5*x.'*Q*x + d*y;
% Solve the problem
fvalP2 = optimize(ConstraintsP2,ObjectiveP2, options);
x = value(x);

%---------------1990 p3-------------------------

%1990 Christodoulos A. Floudas Panos M. Pardalos  test problem 3
% Define variables
x = sdpvar(4,1);
y = sdpvar(9,1);
c = [5 5 5 5]';
Q = 10*eye(4);
d = [-1 -1 -1 -1 -1 -1 -1 -1 -1];

% Define constraints 
ConstraintsP3 = [2*x(1) + 2*x(2) + y(6) + y(7) <= 10; 2*x(1) + 2*x(3) + y(6) + y(8) <= 10; 2*x(2) + 2*x(3) + y(7) + y(8) <= 10; -8*x(1) + y(6) <= 0; -8*x(2) + y(7) <= 0; -8*x(3) + y(8) <= 0; -2*x(4) - y(1) + y(6) <= 0; -2*y(2) - y(3) + y(7) <= 0; -2*y(4) - y(5) + y(8) <= 0;
    (0<=x(1))<=1; (0<=x(2))<=1; (0<=x(3))<=1; (0<=x(4))<=1; (0<=y(1))<=1; (0<=y(2))<=1; (0<=y(3))<=1; (0<=y(4))<=1; (0<=y(5))<=1; (0<=y(9))<=1;0<=y(6);0<=y(7);0<=y(8)];

% Define an objective
ObjectiveP3 =  c.'*x - 0.5*x.'*Q*x + d*y;

% Set some options for YALMIP and solver

% Solve the problem
fvalP3 = optimize(ConstraintsP3,ObjectiveP3, options);
xP3 = value(x);
xP3

%------------------------10--------------------

% sum squares d = 20

x = sdpvar(20,1);
Objective10 = sumsquare(x);
Constraints10 = (-10 <= x) <= 10;

fval10 = optimize(Constraints10,Objective10,options)



%---------------------------9-----------------------


x = sdpvar(2,1);
Objective9 = 0.26*(x(1)^2 + x(2)^2) - 0.48*x(1)*x(2);
Constraints9 = (-10 <= x) <= 10;

fval9 = optimize(Constraints9,Objective9,options)


T = [NaN	0.27  NaN	0.002	0.0064	0.0018	0.002	0.0022	0.0026	0.0044;
    0.238	0.203	0.611	0.126	0.123	0.128	0.13	0.139	0.14	0.417;
    NaN	NaN	NaN	0.075	0.062	0.081	0.077	0.09	0.069	0.063;
    1.068	0.279	0.758	0.123	0.209	0.108	0.096	0.081	0.061	2.463];
T = T.';

logplot = 1;
perf(T,logplot)





%-----------------functions----------------

function perf(T,logplot)
%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More', 
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
% Failures on a given problem are represented by a NaN.
% The optional argument logplot is used to produce a 
% log (base 2) performance plot.
%
% This function is based on the perl script of Liz Dolan.
%
% Jorge J. More', June 2004

if (nargin < 2) logplot = 0; end

%colors  = ['m' 'b' 'r' 'g' 'c' 'k' 'y'];
colors  = [ [0 0 0];[0 0 1];[0 1 0];[0 1 1];[1 0 0];[1 0 1];
            [0 0 0.5];[0 0.5 0];[0 0.5 0.5];[0.5 0 0];[0.5 0 0.5];[0.5 0.5 0];
            [0.3 0.9 0.5];[0.7 0.3 0.5];[0.5 0.3 0.5];[0.3 0.3 0.7];[0.2 0.5 0.9];[0.1 0.4 0.3]
          ];
%lines   = cellstr(char( '-.', '--', ':', '-'));
lines   = [ '-.' '--' ':' '-'];
markers = ['x' '>'  'o' 'v' 's' 'd' '*'];

[np,ns] = size(T);

% Minimal performance per solver

minperf = min(T,[],2);

% Compute ratios and divide by smallest element in each row.

r = zeros(np,ns);
for p = 1: np
  r(p,:) = T(p,:)/minperf(p);
end

if (logplot) r = log2(r); end

max_ratio = max(max(r));

% Replace all NaN's with twice the max_ratio and sort.

r(find(isnan(r))) = 2*max_ratio;
r = sort(r);

% Plot stair graphs with markers.
clf;
for s = 1: ns
 [xs,ys] = stairs(r(:,s),[1:np]/np); 
 temp = floor((s-1)./3)+1;
 %set(gca, 'LineStyle', '-');
 plot(xs,ys,'Color', colors(s,:),'LineWidth' ,1.3, 'Marker', markers(temp), 'LineStyle', '-');
 hold on;
 end
%legend('Strategy1','Strategy2','Strategy3','Strategy4','Strategy5','Strategy6','Strategy7','Strategy8', 'Strategy9','Strategy10','Strategy11','Strategy12','Strategy13','Strategy14', 'Strategy15','Strategy16','Strategy17','Strategy18');
%legend('Strategy10','Strategy11','Strategy12','Strategy13','Strategy14', 'Strategy15','Strategy16','Strategy17','Strategy18');
%legend('QBB', 'QBN',	'QBS',	'QNB',	'QNN',	'QNS',	'QSB',	'QSN',	'QSS',	'NBB',	'NBN',	'NBS',	'NNB',	'NNN',	'NNS',	'NSB',	'NSN',	'NSS');
%legend('NBB',	'NBN',	'NBS',	'NNB',	'NNN',	'NNS',	'NSB',	'NSN',	'NSS');
%legend('QBB', 'QBN',	'QBS',	'QNB',	'QNN',	'QNS',	'QSB',	'QSN',	'QSS');
legend('quadprog',	'kktqp',	'osqp',	'bmibnb');

if (logplot) 
    xlabel(strcat('log_2(\tau)')); 
    ylabel('\rho_s(log_2(\tau))');
else
    xlabel('\tau');
    ylabel('\rho_s(\tau)')
end


% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

axis([ 0 1.1*max_ratio 0 1 ]);

% Legends and title should be added.
end