%this code takes in data from c++ "forward_2D_ML_MCdensity" or related, and
%uses the MC approximation of the time-dependent probability density to
%solve a non-linear least squares problem to generate MC data for the first
%nontrivial forward mode
%
%it also requires an eigenvalue estimate - we will probably have to use
%EDMD
%% Generate training set using LHS
N_train = 64000;
d = 6;

% Latin Hypercube Sampling in [0,1]^d, then scale to [-pi, pi]^d
X_train = lhsdesign(N_train, d);                  % N x 6 matrix in [0,1]^6
X_train = X_train * (2*pi) - pi;                  % scale to [-pi, pi]^6

% Save to .txt file (space-delimited, full precision)
writematrix(X_train, 'Kuramoto_backward_LHS_X_train.txt', ...
    'Delimiter', ' ', ...
    'FileType', 'text');

%% load in the data

%establish parameters as in c++ code

%domain
lowx = -pi;
Spx = 2*pi;
%N = 40;
%hx = Spx / N;

%time
T = 200;
dt = 0.01;
gap = 20;
dt = dt*gap;

%num samples
N_sample = 10000;
X_train = load('Kuramoto_backward_LHS_X_train.txt');
Y_train = load('Kuramoto_backward_LHS_Y_train.txt');

%% Rescaling
Y_train = Y_train/N_sample;

%% plot decaying probability density





obs = Y_train(1:10,:);

%discretize time
tt = 0:dt:T*dt;

%plot the decaying exponential (also the mostly periodic curve)
plot(tt,obs(1:10,:),'LineWidth',1);
set(gca,'fontsize',15)
title('Decaying Oscillation')
xlabel('time t')
ylabel('p(x,t)')
box on
axis square
hold on


%% estimate eigenvalues

%convenience
Tstart = 20;
Tend = 200;
tt = 0:dt:(T-1)*dt;
tt = tt(Tstart:Tend);

N_index = 1000;
N_train = 64000;
index = randsample(N_train, N_index);

freq_est = zeros(1,length(index));
rate_est = zeros(1,length(index));

%here is our data
for i = 1:length(index)
    data = Y_train(index(i),Tstart:Tend);


%find useful statistics about data
    datau = max(data);
    datal = min(data);
    datar = (datau-datal);   
    datam = mean(data); 

%centralize the data
    dataz = data - datau + (datar/2);

%estimating the period
    [~,loc] = findpeaks(dataz,'MinPeakDistance',10);
    my_per = mean(diff(tt(loc)));

%define a function to fit MC data to
    fit  = @(b,x)  exp(b(5)*x).*b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);

%define a least squares cost function
    fcn = @(b) sum((fit(b,tt) - data).^2); 

%solve least squares problem (!!!!)
    s = fminsearch(fcn, [datar;  my_per;  1;  datam; -0.1]);
    freq_est(i) = 2*pi/s(2);
    rate_est(i) = s(5);
    if i == 1
        plot(tt,data,'b',tt,fit(s,tt), 'r','linewidth',3)
        legend({'$\int_{y\in A} (p(y,t) - p_{inv}(y)) dy$', 'direct fitting'},'Interpreter','latex')
        set(gca,'FontSize', 20)
        hold on;
    end
end

index_accept = find(abs(freq_est - 1.0) < 0.05);

%plot the result

%plot(tt,data,'b',tt,fit(s,tt), 'r','linewidth',3)
%legend({'$\int_{y\in A} (p(y,t) - p_{inv}(y)) dy$', 'direct fitting'},'Interpreter','latex')
%set(gca,'FontSize', 20)

mean_freq = mean(freq_est(index_accept))
mean_rate = mean(rate_est(index_accept))


%% solve for eigenfunctions

%we know what the eigenvalues are
%lamR = rate_est;
%lamI = freq_est;

Tstart = 20;
Tend = 200;
tt = 0:dt:(T-1)*dt;
tt = tt(Tstart:Tend);

lamR = -0.1875;
lamI = 1.0;

%organize them
para = [lamR lamI];

%choose some timeslots to look at
array_t = Tstart:1:Tend;
array_t = array_t*dt;

N_train = 64000;
X_train_selected = zeros(N_train,6);
Y_train_R = zeros(N_train,1);
Y_train_I = zeros(N_train,1);

A = zeros(length(array_t),3);
A(:,1) = (exp(lamR*array_t).*cos(lamI*array_t))';
A(:,2) = (exp(lamR*array_t).*sin(lamI*array_t))';
A(:,3) = ones(length(array_t),1);
count = 0;
for i = 1:N_train
    array_y = Y_train(i,Tstart:1:Tend);
    sol = A\array_y';
    SSR = norm(A*sol - array_y')/sqrt(length(array_t));
    if SSR < 0.1
        count = count + 1;
        X_train_select(count,:) = X_train(i,:);
        Y_train_R(count) = sol(1);
        Y_train_I(count) = sol(2);
    end
end
X_train_select = X_train_select(1:count,:);
Y_train_R = Y_train_R(1:count);
Y_train_I = Y_train_I(1:count);
%% Save data

writematrix(X_train_select,'backward_Kuramoto_X_sig_1_5.txt','delimiter','space');
writematrix(Y_train_R,'backward_Kuramoto_REAL_sig_1_5.txt','delimiter','space');
writematrix(Y_train_I,'backward_Kuramoto_IMAG_sig_1_5.txt','delimiter','space');

%% Generate second training set using LHS
N_train = 100000;
d = 6;

% Latin Hypercube Sampling in [0,1]^d, then scale to [-pi, pi]^d
X_train = lhsdesign(N_train, d);                  % N x 6 matrix in [0,1]^6
X_train = X_train * (2*pi) - pi;                  % scale to [-pi, pi]^6

% Save to .txt file (space-delimited, full precision)
writematrix(X_train, 'Kuramoto_backward_LHS_X_train_noY.txt', ...
    'Delimiter', ' ', ...
    'FileType', 'text');

%% Generate testing set 1
N_test = 100;
X_test = zeros(N_test*N_test, 6);
Sp = 2*pi;
ht = Sp/N_test;
for i = 1:N_test
    for j = 1:N_test
        xx = -Sp/2 + (i-1)*ht + ht/2;
        yy = -Sp/2 + (j-1)*ht + ht/2;
        X_test((i-1)*N_test + j, :) = [xx yy 0 0 0 0];
    end
end

writematrix(X_test,'backward_Kuramoto_X_test.txt','delimiter','space');

%% Generate testing set 2
N_test = 200;
X_test = zeros(N_test*N_test, 6);
Sp = 4*pi;
ht = Sp/N_test;
v1 = [1 1 1 1 1 1];
v2 = [1 -1 1 -1 1 -1];
v1 = v1/norm(v1);
v2 = 0.5*v2/norm(v2);
for i = 1:N_test
    for j = 1:N_test
        xx = -Sp/2 + (i-1)*ht + ht/2;
        yy = -Sp/2 + (j-1)*ht + ht/2;
        X_test((i-1)*N_test + j, :) = wrapToPi(xx*v1 + yy*v2);
    end
end

writematrix(X_test,'forward_Kuramoto_X_test2.txt','delimiter','space');

%% Compute the Q-trajectory
eig_R = load('forward_2D_SNIC_REAL.txt');
eig_I = load('forward_2D_SNIC_IMAG.txt');
Traj = load('SNIC_2D_traj.txt');

%grid
x_cor = lowx + hx/2:hx:lowx + Spx - hx/2;
y_cor = lowy + hy/2:hy:lowy + Spy - hy/2;

[XX, YY] = meshgrid(x_cor, y_cor);


R_traj = interp2(XX,YY,eig_R, Traj(:,1), Traj(:,2));
I_traj = interp2(XX,YY,eig_I, Traj(:,1), Traj(:,2));

%% plot Q-trajectory
t_index = 400001:500000;
t_span = t_index*0.001;
subplot(1,3,1)
plot(t_span, R_traj(t_index),'linewidth',1);
title('real part of Q-coordinate')
subplot(1,3,2)
plot(t_span, I_traj(t_index),'linewidth',1);
title('imag part of Q-coordinate')
subplot(1,3,3)
plot(R_traj(t_index),I_traj(t_index),'linewidth',0.5);
xlabel('real')
ylabel('imag')
title('Q-coordinate on the complex plane')

%% Mesh of eigenfunctions
subplot(1,2,1);
mesh(XX,YY,eig_R)
title('real part')
subplot(1,2,2);
mesh(XX,YY,eig_I)
title('imag part')

%% functionals

%nonlinear equation to solve for the forward eigenfunction
function res = three_pt(array_t, array_y, para, var)

%a nonlinear equation of 3 variables
res = zeros(length(array_t),1);

%real and imaginary parts of the eigenfunction
Pr = var(1);
Pim = var(2);

%coefficient from orthogonality
c1 = var(3);

%eigenvalues
lamR = para(1);
lamI = para(2);

%nonlinear equation
for i = 1:length(array_t)
    lamrt = lamR*array_t(i);
    lamit = lamI*array_t(i);
    res(i) = exp(lamrt)*(cos(lamit)*Pr - sin(lamit)*Pim ) ...
        + c1*exp(lamrt)*(cos(lamit)*Pim + sin(lamit)*Pr) - array_y(i);
end

end


