% Version with Integer Programming
% Andrew Ptaszek, Florida Polytechnic University, Spring 2024

% Constants 
t = 24 / (24) ; % 1 minute delta t
N = 24 / t; % Number of steps

% Ambient temperature options
Amb_const = 21;
Amb_periodic = 21 + 4 * sin(pi * ((0:3*N-1)*t - 8) / 12);
Amb_sudden = zeros(3*N, 1);
for n = 1:N
    if n*t <= 12
        Amb_sudden(n, 1) = 21;
    else
        Amb_sudden(n, 1) = 17;
    end
end

option = 3; % Amb temp choice
k = .68; % Expirimental k value
c = 30;

w1 = 200; % Weight against deviation from 37 degrees C
w2 = 1; % Weight against high voltage consumption
w3 = 1; % Weight against freqency of turning the heating element on/off

% linprog init 
cost = [w1*ones(1, N) w2*ones(1,N) 0*ones(1,N) w3*ones(1,N) w3*ones(1,N)]; % (dn, un, xn, vnplus, vnminus)
beq = zeros(2*N, 1); 
bineq = zeros(3*N, 1); 
lb = zeros(5*N, 1);
ub = [Inf(N, 1); ones(N, 1); Inf(N, 1); ones(N, 1); ones(N, 1)]; % un, vnplus, vnminus -> [0, 1]
indicies = 1:(5*N);
integerIndicies = [N+1:2*N, 3*N+1:4*N, 4*N+1:5*N];

% Aeq / Aineq components ====================================
dneq = sparse(2*N, N); 

%uneq setup
rows = [];
for i = 1:2*N
    if mod(i, 2) == 1
        rows = [rows i];
    else
        rows = [rows i i];
    end
end

rows(end) = [];
cols = repelem(2:N, 3);
cols = [1 1 cols];

uneqVals = repmat([-c*t 1 -1], 1, N);
uneqVals(end) = [];
uneq = sparse(rows, cols, uneqVals, 2*N, N);

%xn setup
rows = [];
for i = 1:2*N
    if mod(i, 2) == 1
        rows = [rows i i];
    end
end
rows(1) = [];

cols = repelem(1:N, 2);
cols(end) = [];

xnconst = (k*t) - 1;
xnVals = repmat([1 xnconst], 1, N);
xnVals(end) = [];
xneq = sparse(rows, cols, xnVals, 2*N, N);

%vnplus / vnminus setup
rows = 2:2:(2*N);
cols = 1:N;

vnpluseq = sparse(rows, cols, 1, 2*N, N);
vnminuseq = sparse(rows, cols, -1, 2*N, N);


Aeq = horzcat(dneq, uneq, xneq, vnpluseq, vnminuseq);

% dnineq / unineq / xnineq setup
rows = 1:3*N;
toRemove = mod(rows, 3) == 1;
rows(toRemove) = [];
cols = repelem(1:N, 2);

dnineq = sparse(rows, cols, -1, 3*N, N);
unineq = sparse(3*N, N);

xnVals = (-1).^(0:2*N-1);
xnineq = sparse(rows, cols, xnVals, 3*N, N);

% vnineq

rows = 1:3:(3*N);
cols = 1:N;
vnineq = sparse(rows, cols, 1, 3*N, N); %vnplusineq = vnminusineq

Aineq = horzcat(dnineq, unineq, xnineq, vnineq, vnineq);
% ========================================================

%beq setup for different ambient temps
for n = 1:2*N
    if mod(n, 2) == 1
        if option == 1
            beq(n, 1) = t*k*Amb_const;
        end
        if option == 2
            beq(n, 1) = t*k*Amb_periodic(1, n);
        end
        if option == 3
            beq(n, 1) = t*k*Amb_sudden(n, 1);
        end
    
        if n == 1
            if option == 1
                beq(n, 1) = t*k*Amb_const - 21*((k*t) -1);
            end
            if option == 2
                beq(n, 1) = t*k*Amb_periodic(1, n) - 21*((k*t) -1);
            end
            if option == 3
                beq(n, 1) = t*k*Amb_sudden(n, 1) - 21*((k*t) -1);
            end
        end
    else 
        beq(n, 1) = 0;
    end
end 

%bineq setup
for n = 1:3*N 
    if mod(n, 3) == 1
        bineq(n, 1) = 1;
    end

    if mod(n, 3) == 2
        bineq(n, 1) = 37;
    end

    if mod(n, 3) == 0
        bineq(n, 1) = -37;
    end
end

[x, obj] = intlinprog(cost, integerIndicies, Aineq, bineq, Aeq, beq, lb, ub);
disp(x)

stepsNeeded = N;

% Extracting dn, un, xn, vnplus, and vnminus
dn = x(1:stepsNeeded);
un = x(N+1:N+stepsNeeded);
xn = x(2*N+1:2*N+stepsNeeded);
vnplus = x(3*N+1:3*N+stepsNeeded);
vnminus = x(4*N+1:4*N+stepsNeeded);

% Plotting dn, un, xn, vnplus, and vnminus for the first stepsNeeded steps
figure; 

subplot(5, 1, 1);
plot((1:stepsNeeded)*t, dn, 'b-', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('dn');
title('Plot of dn over Time');

subplot(5, 1, 2);
plot((1:stepsNeeded)*t, un, 'r-', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('un');
title('Plot of un over Time');

subplot(5, 1, 3);
plot((1:stepsNeeded)*t, xn, 'g-', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('xn');
title('Plot of xn over Time');

subplot(5, 1, 4);
plot((1:stepsNeeded)*t, vnplus, 'm-', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('vnplus');
title('Plot of vnplus over Time');

subplot(5, 1, 5);
plot((1:stepsNeeded)*t, vnminus, 'c-', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('vnminus');
title('Plot of vnminus over Time');

sgtitle('Solution Components over Time');