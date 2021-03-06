clear; close all; clc;
%% forward simulation
%% Load models and drive cycle
% load Profile_adapted_2
load SS_LPV
load RC_1st_ARX_smooth_4
load curve_t2_smooth
load variations

%% Scale capacity
C0 = C0/10;
C0_varied = C0*C0_var;

%% make smooth
A.Method = 'spline';
B.Method = 'spline';
D.Method = 'spline';

%% Scenario
Rb = 0.5;
y_low = 2.6;
u_max = 0.5;
s0 = 1;

%% Load maximal results
addpath('C:\Users\s136161\Documents\ACC_2019_Dual_results\trust_region_in_constraint')
load Succesful_result_Edited_LPV_longer_3

opts             = optimset('fminsearch');
opts.TolX        = 1.e-12;

%% fmincon settings
options_fmincon = optimset('fmincon');
options_fmincon.Display = 'off'; 
options_fmincon.FunctionTolerance = 1.e-8;
options_fmincon.ConstraintTolerance = 1.e-8;

options = optimoptions('quadprog');
options.Display = 'off';
options.Diagnostics = 'off';

clear y x

%% Create reference
% start = 477;

name = {'V','SoC','C','CV'};
savefolder = 'MPC_results\';

for method = 1:4
    for hor = 1:3

%% Initialisation for known starting point
% Initial condition
x = [s0;0];
x = repmat(x,N,1);
As = [];
Bs = [];

% System definition
C_avg = mean(C0*C0_var);
y_low_violated = 0;

for k = 1:length(P)
%     disp(k)
    for n = 1:N
        As(2*n-1:2*n,2*n-1:2*n) = [1 0; 0 A(x(2*n-1,k))^AB_var(n)];
        Bs(2*n-1:2*n,1:n+1) = [1/(C0*C0_var(n)) zeros(1,n-1) 1/(C0*C0_var(n));B(x(2*n-1,k))*AB_var(n) zeros(1,n-1) B(x(2*n-1,k))*AB_var(n)];
        R0(n) = D(x(2*n-1,k))*D_var(n);
        Vemf(n) = EMF(x(2*n-1,k));
    end

    a = sum(R0);
    b = sum(Vemf)+repmat([0 1],1,N)*x(:,k);
    c = -P(k);
    w(k) = (-b + sqrt(b^2-4*a*c))/(2*a);

    lb = w(k)-3*bal_lim;
    ub = w(k)+3*bal_lim;
    H_cons_1 = zeros(N+1,N+1);
    H_cons_1(1,1) = 2*a;
    k_cons_1 = zeros(N+1,1);
    k_cons_1(1,1) = b;
    d_cons_1 = -P(k);
    H_cons_2 = zeros(N+1,N+1);
    k_cons_2 = zeros(N+1,1);
    d_cons_2 = 0;

    for n=2:N+1
        lb(n) = -bal_lim;
        ub(n) = +bal_lim;
        H_cons_1(1,n) = R0(n-1);
        H_cons_1(n,1) = R0(n-1);
        if n>1
            H_cons_2(1,n) = R0(n-1);
            H_cons_2(n,1) = R0(n-1);
            H_cons_2(n,n) = 2*(R0(n-1)+Rb);
            k_cons_2(n,1) = Vemf(n-1)+x(2*(n-1),k);
        end
    end
    
%     lb = w(k)-3*bal_lim;
%     ub = w(k)+3*bal_lim;
%     H_cons_1 = zeros(N+1,N+1);
%     H_cons_1(1,1) = 2*a;
%     k_cons_1 = zeros(N+1,1);
%     k_cons_1(1,1) = sum(Vemf)+repmat([0 1],1,N)*x(:,k);
%     d_cons_1 = -P(k);
%     H_cons_2 = zeros(N+1,N+1);
%     k_cons_2 = zeros(N+1,1);
%     d_cons_2 = 0;
% 
%     for n=2:N+1
%         lb(n) = -bal_lim;
%         ub(n) = +bal_lim;
%         H_cons_1(1,n) = R0(n-1);
%         H_cons_1(n,1) = R0(n-1);
%         if n>1
%             H_cons_2(1,n) = R0(n-1);
%             H_cons_2(n,1) = R0(n-1);
%             H_cons_2(n,n) = 2*(R0(n-1)+Rb);
%             k_cons_2(n,1) = Vemf(n-1)+x(2*(n-1),k);
%         end
%     end

    U_bal(:,k) = u_opt(2*k,:)';
    u0 = zeros(N+1,1);
    u0(1,1) = w(k);
    nonlconstr = @(x)quadconstr(x,H_cons_1,k_cons_1,d_cons_1,H_cons_2,k_cons_2,d_cons_2);
    x_opt = fmincon(@(x_opt) MPC_acceptable_bal_currents(x_opt,N,U_bal(:,k)),...
        u0,[],[],[],[],lb,ub,nonlconstr,options_fmincon);
    w(k) = x_opt(1); U_bal(:,k) = x_opt(2:end);
   
    U = [w(k); U_bal(:,k)];

    x(:,k+1) = As*x(:,k) + Bs*U;
    for n=1:N
        y(n,k) = x(2*n,k) + [R0(n) zeros(1,n-1) R0(n) zeros(1,N-n)]*(U) + Vemf(n);
    end
    if max(y(:,k)<2.6)
        ending_search = k;
        y_low_violated = 1;
        break
    else
        ending_search = k;
    end
end

figure;hold on;
for n=1:N
   plot(y(n,1:ending_search)) 
end

figure;hold on;
for n=1:N
   plot(x(2*n-1,1:ending_search)) 
end
% 
% figure;hold on;
% for n=1:N
%    plot(C_left(n,1:ending_search)) 
% end

figure;hold on;
for n=1:N
   plot(u_opt(2:2:end,n),'--') 
   plot(U_bal(n,1:ending_search)) 
end
legend

w_base = w;
x_base = x;
y_base = y;

close all

%% Scenario length
scenario_length = numel(P);
horizon_choices = [10,25,50];
horizon = horizon_choices(hor);
hn = horizon;

%% Tuning
Q_par = [1;5];
s_bottom = 0.1;

%%
% Initial condition
x_pred = x(:,1:hn);
clear x y ubal C_left

x = [s0;0];
x = repmat(x,N,1);
% x = repmat(x,1,20);
As = [];
Bs = [];

M = eye(N*hn);
M = M+repmat(-1/N*eye(hn),N);

% System definition
for k = 1:length(P)   
    %% Compute expected pack current
    if k==1
        for ip=1:hn
            for n=1:N
                u_opt(2*ip,n) = 0;
            end
        end
    else
        for ip=1:hn
            for n=1:N
                u_opt(2*ip,n) = u_opt_MPC(n);
            end
        end
    end
    if k<length(P)
        clear i
        w_pred = zeros(hn,1)+i;
        cut_back = 1;
        while isreal(w_pred)==0 || cut_back==1
            clear w_pred
            cut_back = cut_back-1;
            hn = min(max(1,hn+cut_back),numel(P)-k+1);
            for j=1:hn
                for n = 1:N
                    Af = A(x_pred(2*n-1,j))^AB_var(n);
                    Bf = B(x_pred(2*n-1,j))*AB_var(n);
                    As(2*n-1:2*n,2*n-1:2*n) = [1 0; 0 Af];
                    Bs(2*n-1:2*n,1:n+1) = [1/(C0*C0_var(n)) zeros(1,n-1) 1/(C0*C0_var(n));Bf zeros(1,n-1) Bf];
                    R0(n) = D(x_pred(2*n-1,j))*D_var(n);
                    Vemf(n) = EMF(x_pred(2*n-1,j));
                end
                a = sum(R0);
                b = sum(Vemf)+repmat([0 1],1,N)*x_pred(:,j)+R0*u_opt(j,:)';
                c = -P(k+j-1);
                w_pred(j,1) = (-b + sqrt(b^2-4*a*c))/(2*a);
            end
        end

        if k==1
            for ip=1:hn
                for n=1:N
                    u_opt(2*ip-1,n) = w_pred(ip);
                end
            end
        else
            for ip=1:hn
                for n=1:N
                    u_opt(2*ip-1,n) = w_pred(ip);
                end
            end
        end

        clear x_pred
        x_pred(:,1) = x(:,k);
        for j=1:hn
            if k == 1
                U = [w_pred(j); zeros(10,1)];
            else
                U = [w_pred(j); u_opt_MPC];
            end
            for n = 1:N
                Af = A(x_pred(2*n-1,j))^AB_var(n);
                Bf = B(x_pred(2*n-1,j))*AB_var(n);
                As(2*n-1:2*n,2*n-1:2*n) = [1 0; 0 Af];
                Bs(2*n-1:2*n,1:n+1) = [1/(C0*C0_var(n)) zeros(1,n-1) 1/(C0*C0_var(n));Bf zeros(1,n-1) Bf];
            end
            x_pred(:,j+1) = As*x_pred(:,j) + Bs*U;
        end

        %% Linearisation 
        epsilon = 0.005;
        range = linspace(-epsilon,epsilon,3);
        clear beta0 beta1 alpha0 alpha1
        for n=1:N
            for j=1:hn      
                for i=1:numel(range)
                    x_min = x_pred(2*n-1,j)+range(i);

                    Y(i) = EMF(x_min)+D(x_min)*D_var(n)*(u_opt(2*j-1,n)+u_opt(2*j,n))-D(x_pred(2*n-1,j))*D_var(n)*(u_opt(2*j-1,n)+u_opt(2*j,n));  
                    y_perfect(i) = EMF(x_min)+x_pred(2*n,j)+D(x_min)*D_var(n)*(u_opt(2*j-1,n)+u_opt(2*j,n));

                    Af = A(x_min)^AB_var(n);
                    Bf = B(x_min)*AB_var(n);
                    x_perfect(i) = Af*x_pred(2*n,j) + Bf*(u_opt(2*j-1,n)+u_opt(2*j,n));

                    Af = A(x_pred(2*n-1,j))^AB_var(n);
                    Bf = B(x_pred(2*n-1,j))*AB_var(n);
                    Y_X(i) = x_perfect(i)-(Af*x_pred(2*n,j) + Bf*(u_opt(2*j-1,n)+u_opt(2*j,n)));
                end
                X = x_pred(2*n-1,j)+range;
                gradY = (Y(end)-Y(1))/(sum(diff(range)));
                beta0(n,j) = -gradY*x_pred(2*n-1,j)+Y(ceil(numel(range)/2));
                beta1(n,j) = gradY;

                gradX = (Y_X(end)-Y_X(1))/(sum(diff(range)));
                alpha0(n,j) = -gradX*x_pred(2*n-1,j);
                alpha1(n,j) = gradX;

                xboundup(j,n) = min(1,x_pred(2*n-1,j)+epsilon);
                xbounddown(j,n) = max(0,x_pred(2*n-1,j)-epsilon);

                alpha2 = ((A(x_pred(2*n-1,j)+epsilon)^AB_var(n)-(A(x_pred(2*n-1,j)-epsilon))^AB_var(n))*x_pred(2*n,j))/(2*epsilon) + ((u_opt(2*j-1,n)+u_opt(2*j,n))*(B(x_pred(2*n-1,j)+epsilon)*AB_var(n)-(B(x_pred(2*n-1,j)-epsilon)*AB_var(n))))/(2*epsilon);
                alpha3 = -(alpha2)*x_pred(2*n-1,j);

                beta2 = ((EMF(x_pred(2*n-1,j)+epsilon)-(EMF(x_pred(2*n-1,j)-epsilon)))/(2*epsilon) + ((u_opt(2*j-1,n)+u_opt(2*j,n))*(D(x_pred(2*n-1,j)+epsilon)*D_var(n)-(D(x_pred(2*n-1,j)-epsilon)*D_var(n))))/(2*epsilon));
                beta3 = -(beta2)*x_pred(2*n-1,j);
           end
        end

        %% Construct prediction matrices
        clear Aineq X_predictor_A bineq X_predictor_b s_predictor_A s_predictor_b Q_predictor_A Q_predictor_b
        for n=1:N
            [phi,gamma,Cm,Dm,Xiu] = MPC_predictMatrix_constant_LTV_N_lin_full(beta1(n,:),hn,n,C0_varied,AB_var,D_var,A,B,D,x_pred,alpha1,alpha0);

            x0 = x_pred(2*n-1:2*n,1);

            Aineq(hn*(n-1)+1:hn*n,n) = -[Cm*gamma(:,1:2:2*hn)+Dm(:,1:2:2*hn)]*ones(hn,1);
            bineq(hn*(n-1)+1:hn*n,1) = -y_low + Cm*(phi*x0)+ beta0(n,:)'+ Xiu'+[Cm*gamma(:,1:2:2*hn)+Dm(:,1:2:2*hn)]*w_pred(1:hn);
            
            X_predictor_A(2*hn*(n-1)+1:2*hn*n,n) = [-gamma(1:2*hn,1:2:2*hn)*ones(hn,1)];
            X_predictor_b(2*hn*(n-1)+1:2*hn*n,1) = phi(1:end,:)*x0+gamma(1:end,1:2:end)*w_pred(1:hn);
            
            s_predictor_A(hn*(n-1)+1:hn*n,n) = [-gamma(1:2:2*hn,1:2:2*hn)*ones(hn,1)];
            s_predictor_b(hn*(n-1)+1:hn*n,1) = phi(1:2:end,:)*x0+gamma(1:2:end,1:2:end)*w_pred(1:hn);
            
            Q_predictor_A(hn*(n-1)+1:hn*n,n) = C0_varied(n)*[-gamma(1:2:2*hn,1:2:2*hn)*ones(hn,1)];
            Q_predictor_b(hn*(n-1)+1:hn*n,1) = C0_varied(n)*((phi(1:2:end,:)*x0+gamma(1:2:end,1:2:end)*w_pred(1:hn))-(s_bottom+(s0-s_bottom)*((numel(P)-k-[-1:hn-2])./(numel(P)))')); 
            if method == 4
                Q_predictor_b(hn*(n-1)+1:hn*n,1) = C0_varied(n)*((phi(1:2:end,:)*x0+gamma(1:2:end,1:2:end)*w_pred(1:hn))-s_bottom); 
            end
        end
        
        %% Construct M
        if hn*N < numel(diag(M)) || k ==1
            M = eye(N*hn);
            M = M+repmat(-1/N*eye(hn),N);  
            Q = diag(repmat([Q_par(1)*ones(1,hn-1),Q_par(2)*ones(1,1)],1,N));
        end
        
        if method == 1
            %% Cost function voltage balancing
            R = hn/25*0.01;
            H = (-Aineq)'*M'*Q*M*(-Aineq)+R*eye(N);
            F = 2*(bineq+y_low)'*M'*Q*M*(-Aineq);
            
            H_end = (-Aineq)'*M'*Q*M*(-Aineq)+R*eye(N);
            F_end = 2*(bineq+y_low)'*M'*Q*M*(-Aineq);
        
        elseif method == 2
            %% Cost function SoC balancing
            R = (hn/25)*1e-3;
            H = (-s_predictor_A)'*M'*Q*M*(-s_predictor_A)+R*eye(N);
            F = 2*(s_predictor_b)'*M'*Q*M*(-s_predictor_A);
            
            R = hn/25*0.01;
            H_end = (-Aineq)'*M'*Q*M*(-Aineq)+R*eye(N);
            F_end = 2*(bineq+y_low)'*M'*Q*M*(-Aineq);

        elseif method == 3
            %% Cost function charge balancing reference governor
            R = 0;
            H = (-Q_predictor_A)'*M'*(Q)*M*(-Q_predictor_A)+R*eye(N);
            F = 2*(Q_predictor_b)'*M'*(Q)*M*(-Q_predictor_A);
            
            R = hn/25*0.01;
            H_end = (-Aineq)'*M'*Q*M*(-Aineq)+R*eye(N);
            F_end = 2*(bineq+y_low)'*M'*Q*M*(-Aineq);
            
        elseif method == 4
            %% simple charge balancing
            if hor == 1
                R = 6e4;
            elseif hor == 2
                R = 2.5^1.9*6e4;
            else
                R = 5^1.9*6e4;
            end
            H = (-Q_predictor_A)'*M'*(Q)*M*(-Q_predictor_A)+R*eye(N);
            F = 2*(Q_predictor_b)'*M'*(Q)*M*(-Q_predictor_A);
            
            R = hn/25*0.01;
            H_end = (-Aineq)'*M'*Q*M*(-Aineq)+R*eye(N);
            F_end = 2*(bineq+y_low)'*M'*Q*M*(-Aineq);
        end

%         if hn<horizon
%             R = 0.01;
%             H = (-Aineq)'*M'*Q*M*(-Aineq)+R*eye(N);
%             F = 2*(bineq+y_low)'*M'*Q*M*(-Aineq);
%         end

        %% Additional constraints
        Aeq = [ones(1,N)];
        beq = [0];

        lb = -u_max*ones(N,1);
        ub = u_max*ones(N,1);

        %% Tests
    %         y_matrix = -Aineq(1:hn*N,1:N)*u_opt_MPC+ bineq(1:hn*N,1)+y_low;
%             x_matrix_test_first = -Aineq(3*N*hn+1:3*N*hn+hn*N,1:N)*u_opt+bineq(3*N*hn+1:3*N*hn+hn*N,1);

        %% Optimization
        [u_opt_MPC,fval,exitflag,output] = quadprog((H+H')/2,F,[Aineq],[bineq],[Aeq],[beq],lb,ub,[],options);
        
        
%         [u_opt_MPC,fval,exitflag,output] = quadprog((H+H')/2,F,[],[],[Aeq],[beq],lb,ub,[],options);
        
        if isempty(u_opt_MPC) == 1
            disp('Bounds lost!')
            [u_opt_MPC,fval,exitflag,output] = quadprog((H_end+H_end')/2,F_end,[],[],[Aeq],[beq],lb,ub,[],options);
        end
        
        ubal_pred(:,k) = u_opt_MPC(:,1);
        
        %% compute C_left
        temp = M*(-Q_predictor_A*u_opt_MPC+Q_predictor_b);
        C_left_avg(:,k) = temp(1:hn:end);
        
        %% Predict new x
        clear x_pred_check y_pred
        for n=1:N
            x_pred_check(2*n-1:2*n,:) = reshape(-X_predictor_A(2*hn*(n-1)+1:2*hn*n,n)*u_opt_MPC(n)+X_predictor_b(2*hn*(n-1)+1:2*hn*n,1),[2,hn]);
            y_pred(n,:) = reshape(-Aineq(hn*(n-1)+1:hn*n,n)*u_opt_MPC(n)+ bineq(hn*(n-1)+1:hn*n,1)+y_low,[1,hn]);
        end
    else
        u_opt_MPC = zeros(10,1);
    end

    %% Compute real x
    for n = 1:N
        Af = A(x(2*n-1,k))^AB_var(n);
        Bf = B(x(2*n-1,k))*AB_var(n);
        As(2*n-1:2*n,2*n-1:2*n) = [1 0; 0 Af];
        Bs(2*n-1:2*n,1:n+1) = [1/(C0*C0_var(n)) zeros(1,n-1) 1/(C0*C0_var(n));Bf zeros(1,n-1) Bf];
        R0(n) = D(x(2*n-1,k))*D_var(n);
        Vemf(n) = EMF(x(2*n-1,k));
    end
    
    a = sum(R0);
    b = sum(Vemf)+repmat([0 1],1,N)*x(:,k);
    c = -P(k);
    w(k) = (-b + sqrt(b^2-4*a*c))/(2*a);

    % Compute closest acceptable balancing currents
    lb = w(k)-3*bal_lim;
    ub = w(k)+3*bal_lim;
    H_cons_1 = zeros(N+1,N+1);
    H_cons_1(1,1) = 2*a;
    k_cons_1 = zeros(N+1,1);
    k_cons_1(1,1) = sum(Vemf)+repmat([0 1],1,N)*x(:,k);
    d_cons_1 = -P(k);
    H_cons_2 = zeros(N+1,N+1);
    k_cons_2 = zeros(N+1,1);
    d_cons_2 = 0;

    for n=2:N+1
        lb(n) = -bal_lim;
        ub(n) = +bal_lim;
        H_cons_1(1,n) = R0(n-1);
        H_cons_1(n,1) = R0(n-1);
        if n>1
            H_cons_2(1,n) = R0(n-1);
            H_cons_2(n,1) = R0(n-1);
            H_cons_2(n,n) = 2*(R0(n-1)+Rb);
            k_cons_2(n,1) = Vemf(n-1)+x(2*(n-1),k);
        end
    end

    u0 = zeros(N+1,1);
    u0(1,1) = w(k);
    nonlconstr = @(x)quadconstr(x,H_cons_1,k_cons_1,d_cons_1,H_cons_2,k_cons_2,d_cons_2);
    x_opt = fmincon(@(x_opt) MPC_acceptable_bal_currents(x_opt,N,u_opt_MPC),...
        u0,[],[],[],[],[lb],[ub],nonlconstr,options_fmincon);
    w(k) = x_opt(1); ubal(:,k) = x_opt(2:end);

    U = [w(k); ubal(:,k)];

    for n=1:N
        C_left(n,k) = (x(2*n-1,k)-0.11)*C0*C0_var(n);
    end

    x(:,k+1) = As*x(:,k) + Bs*U;
    for n=1:N
        y(n,k) = x(2*n,k) + [R0(n) zeros(1,n-1) R0(n) zeros(1,N-n)]*(U) + Vemf(n);
    end
    if max(y(:,k)<2.6)
        ending = k-1;
        break
    else
        ending = k;
    end
end
tosave = strcat(savefolder,name(method),'_',num2str(horizon));
save(string(tosave),'y','ubal','x','w')
    end
    
end

fprintf('Voltage spread is %4.3f V.\n', max(y(:,end)-min(y(:,end))))
fprintf('Voltage spread of base result is %4.3f V.\n', max(y_base(:,end)-min(y_base(:,end))))

figure;hold on;grid on;
for i=1:10
    plot(y(i,:))
    plot(y_base(i,:),'--')
end

figure;hold on;grid on;
plot(sum(y(:,:)).*w)
plot(sum(y_base(:,:)).*w_base)
plot(P)


% figure;hold on;grid on;
% for i=1:10
%     plot(y_pred(i,:))
% end
% legend

figure;hold on;
for i=1:10
    plot(x(2*i-1,1:end-1))
    plot(x_base(2*i-1,1:end-1),'--')
end

figure;hold on;grid on;
plot(w(1,:))
plot(w_base(1,:),'--')

figure;hold on;
for n=1:N
   plot(C_left(n,:)) 
end

figure;hold on;
for n=1:N
   plot(C_left_avg(n,:)) 
end

figure;hold on;
for n=1:N
   plot(ubal(n,:)) 
   plot(U_bal(n,:),'--')
%    plot(ubal_pred(n,:),'--')
end
legend


