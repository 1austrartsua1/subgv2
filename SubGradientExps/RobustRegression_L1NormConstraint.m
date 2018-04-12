%% Robust Regression Subject to an ell_1 norm constraint
%% working on AvA
fprintf('Hi There');
% min \sum_{i=1}^m |e_i^{\top} x- b_i| : \|x\|_1\leq \tau
%Good experiments
%-------------Exp1----------------
%m=100, dim = 50, tau = 0.5, 
%c_theta = 2e1;
%epsilon=1e-8;
%beta=4
%theiralpha=2
%---------------------------------
%-------------Exp2----------------
%m=100, dim = 50, tau = 2, 
%c_theta = 1e1;
%epsilon=1e-8;
%beta=4
%theiralpha=2
%---------------------------------
%---------------------------------
%-----------Exp3------------------
%m=100, dim = 50, tau = 0.5, 
%c_theta = 2e1;
%epsilon=1e-12;
%beta=4
%theiralpha=2
%---------------------------------

clear all;close all;
%format long
addpath(genpath(''));

m=100;
dim=50;
E=randn(m,dim);
b=randn(m,1);
Esvd=svd(E);
Lmax = Esvd(1);

x_init = randn(dim,1);
tau = 1;
x_init= L1BallProjection(x_init,tau);

v_in=randn(dim,1);



%% cvx

cvx_begin
   variable x(dim);
   maximize  -norm(E*x-b,1);
   subject to
      norm(x,1)<=tau;
cvx_end

f_opt = -cvx_optval;
x_opt = x;

%% regular subgradient
%profile on
Run_decay=1;

tic
T=3e4;

fprintf('\nRunning Subgradient with Decaying Stepsize %d iterations...\n',T);

A=0.1;d=1;
[x,f_sgDecay1,d_sq_decayd1] = decaySG_RobustRegr_L1Constraint(T,x_init,A,d,E,b,tau,x_opt);

A=0.01;d=0.5;
[x,f_sgDecay05,d_sq_decayd05,d_sq_hat] = decaySG_RobustRegr_L1Constraint(T,x_init,A,d,E,b,tau,x_opt);
toc

semilogy(1:T,f_sgDecay1-f_opt,1:T,f_sgDecay05-f_opt)
%T=1000;A=1;d=0.75;
%[x,f_sgDecay3] = decaySG_RobustRegr_L1Constraint(T,x_init,A,d,E,b,tau);
%% Save data for the nonsummable plot




%%

%% DS-SG NONERG Descending Staircase
tic
beta=4;
%fudge=0.5;
B = Esvd(1)*sqrt(m);

min_norm = norm(E(1,:));
for i=2:m
   norm_of_row = norm(E(i,:));
   if(min_norm >norm_of_row)
       min_norm=norm_of_row;
   end
end
c_true = min_norm;
%B=0;
%for i=1:m
%   B=B+norm(E(i,:),2); 
%end
c_theta = 100;
epsilon=1e-8;
p_1 = 2*tau;%p_1>=norm(x_init-x_opt,2)^2;
M= ceil(log2(p_1/epsilon)/log2(beta));M=1200;
K=ceil(sqrt(2*beta)*log(2*beta)*(B/c_theta)^2);
%K=2800;
%M=3;
alpha_1 = sqrt(2)*c_theta*sqrt(p_1)/(sqrt(beta)*B^2);
total_sg_iters = M*K;
Ergodic=0;

fprintf('\nRunning DS-SG, M= %d and K = %d. Total of %d iterations...\n'...
    ,M,K,total_sg_iters);
Early_stop=0;
[x,f_DSSG,d_sq_DS,f_erg,sparsity,sparsity_erg,d_erg_sq] = ...
DS_SG_RobustRegre_L1Constraint(beta,M,alpha_1,K,x_init,E,b,tau,x_opt,Ergodic);
                              
toc
semilogy(d_sq_DS)
%% performance dsq no ergodic
perform_dsqNoErg=1;
if(perform_dsqNoErg)
Tplot=min([total_sg_iters T]);

%Tplot=min([T,total_sg_iters]);

end
%%
%semilogy(d_sq_1)
%%
%plot(flag

%%%%%%%%%%%%%%%%%%%%%%%% RSG Ergodic Descending Staircase

runErg=1;
if(runErg)
tic
their_alpha=2;
beta=their_alpha^2;
%epsilon=1e-12;

B = Esvd(1)*sqrt(m);
c_thetaRSG =64;
epsilon_2=sqrt(epsilon)/2;
epsilon_0=B*sqrt(p_1);
alpha_1 = epsilon_0/(their_alpha*B^2);
M= ceil(log2(epsilon_0/epsilon_2));M=2e3;
%M=4e2;
K=ceil(((B*their_alpha)/c_thetaRSG)^2);
total_sg_iters2 = M*K;
Ergodic=1;

fprintf('\nRunning RSG with M=%d and K=%d, total of %d iterations...\n',M,K,total_sg_iters2);
[x,f_rsg,d_sq,f_ergs,sparsity0,sparsity_erg,d_sq_rsg] = ...
DS_SG_RobustRegre_L1Constraint(beta,M,alpha_1,K,x_init,E,b,tau,x_opt,Ergodic);
                              
toc   
end
semilogy(d_sq_rsg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% R^2SG

K_init = ceil(K);
increase_factor=1.05;
ThetEq = 1-0.5*log2(1.35);
num_rounds=25;
c_thetaRSG =120;
K=ceil(((B*their_alpha)/c_thetaRSG)^2);
M= ceil(log2(epsilon_0/epsilon_2));
R2SG_num_iters = round(M*(1-increase_factor^num_rounds)/(1-increase_factor)*K_init);
fprintf('\nRunning R2SG with M=%d, K_init=%d, num_rounds=%d, for a total of approx %d iterations...\n',M,K_init,num_rounds,R2SG_num_iters);
[dsq_r2sg,f_r2sg] = R2SG(beta,M,alpha_1,K_init,x_init,E,b,tau,x_opt,num_rounds,increase_factor);

T_r2sg = length(dsq_r2sg);
%semilogy(dsq_r2sg);
%semilogy(f_r2sg-f_opt);
%semilogy(mycummin(f_r2sg)-f_opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DS2-SG

run_Adapt=1;
if(run_Adapt)
tic
%epsilon=1e-12;
M = ceil(log2(p_1/epsilon));
c_1 = B; 
Omega = p_1;
T_Ad = 5;
min_c = c_1*2^(-(T_Ad-1));
K = sqrt(beta/2)*log(2*beta)*(B^2)/(c_1^2);
t_sgs_3 = (1/3)*K*M*(4^(T_Ad)-1);
fprintf('\nRunning D2SG-SG, %d total SG iters...\n',ceil(t_sgs_3));
[xtilde,f_ds2sg,d_sq_D2DSG,sparsity_Ad] =...
Ad_DS_SG_RobustRegre_L1Constraint(beta,B,M,c_1,Omega,x_init,T_Ad,...
E,b,tau,x_opt,K);
toc
Total_ad = length(f_ds2sg);
end
semilogy(f_ds2sg-f_opt)
%% Shor's Method
c_thetaShor=64;
Tshor = 50000;
fprintf('\nRunning Shors Method, %d total SG iters...\n',Tshor);
[d_sq_shor,f_shor,hs] = Shor(c_thetaShor,B,p_1,E,b,tau,x_init,Tshor,x_opt);

%semilogy(f_shor - f_opt)
semilogy(sqrt(d_sq_shor))
%Tshort=2500;
%semilogy(1:Tshort,B*hs(1:Tshort)/c_theta,1:Tshort,sqrt(d_sq_shor(1:Tshort)))





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% STOP RUNNING ALGORITHMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% START PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Shor vs DS-SG vs RSG vs decay
plotDistNoAdapt=0;
if(plotDistNoAdapt)
%if(loadData_NoAdapt==0)
%   Tplot=min([total_sg_iters Tshor total_sg_iters2 T]);
%end
%Tplot=5e4;
%Tplot = min([length(d_sq_decayd1) length(d_sq_D2DSG)]);
Grid=1;
start1=1;
start2=1;
start3=1;
start4=1;
start5=1;
h=semilogy(...
start1:Grid:Tplot,sqrt(d_sq_DS(start1:Grid:Tplot)),'-',...
start2:Grid:Tplot,sqrt(d_sq_shor(start2:Grid:Tplot)),'-',...
start3:Grid:Tplot,sqrt(d_sq_decayd1(start3:Grid:Tplot)),'-',...
start4:Grid:Tplot,sqrt(d_sq_decayd05(start4:Grid:Tplot)),'-',...
start5:Grid:Tplot,sqrt(d_sq_rsg(start5:Grid:Tplot)),'-',...
...1:Tplot,sqrt(d_sq_DS(1:Tplot)),'-',...
...%1:Tplot,d_sq_shor(1:Tplot),'k-',...
...1:Tplot,sqrt(d_sq_decayd1(1:Tplot)),'-',...
...1:Tplot,sqrt(d_sq_D2DSG(1:Tplot)),'-'...
...%1:Tplot,d_sq_decayd05(1:Tplot),'k-'...
...1:Tplot,d_sq_rsg(1:Tplot),'k-'...
'LineWidth',2.5);
set(h,'MarkerSize',17);
legend('DS-SG','Shor','\alpha_k=0.1 k^{-1}','\alpha_k=0.01 k^{-0.5}','RSG')
%defense legend v1
%legend('DS-SG (our method)','Shor (1979)','(\alpha_1,p)=(0.1,1)','(\alpha_1,p)=(0.01,0.5)','Yang (2015)')
%defense legend
%legend('\alpha_k=0.1 k ^{-1}','DS2-SG')
xlabel('num iterations k');
ylabel('$d(x_k,\mathcal{X})$',...
'Interpreter', 'Latex', 'FontSize', 15, 'Color', 'k');
axis([0 Tplot 2e-7 1e2])
grid on
end

%%
%Grid=[1 10 100 1000 1e4];
%T=max([length(d_sq_decayd1) length(d_sq_decayd05)]); 
%h=loglog(...
%Grid,d_sq_decayd1(Grid),'ko',...
%Grid,d_sq_decayd05(Grid),'kv',...
%1:T,d_sq_decayd1,'k-',...
%1:T,d_sq_decayd05,'k-'...
%);
%set(h,'MarkerSize',20)
%legend('(A,d)=(0.1,0.99)','(A,d)=(0.01,0.5)')
%grid on
%xlabel('Subgradient Evals');
%ylabel('$d(x_k,\mathcal{X})^2$',...
%'Interpreter', 'Latex', 'FontSize', 15, 'Color', 'k');
%axis([1 1e5 1e-10 10])
%% Save data for Plotting DistNoAdapt
%saveData_NoAdapt=0;
%if(saveData_NoAdapt)
%SaveDataStruct.Tplot=Tplot;
%SaveDataStruct.d_sq_DS=d_sq_DS;
%SaveDataStruct.d_sq_shor=d_sq_shor;
%SaveDataStruct.d_sq_decayd1=d_sq_decayd1;
%SaveDataStruct.d_sq_decayd05=d_sq_decayd05;
%SaveDataStruct.d_sq_rsg=d_sq_rsg;
%SaveDataStruct.c_theta=c_theta;
%SaveDataStruct.c_thetaShor=c_thetaShor;
%SaveDataStruct.c_thetaRSG=c_thetaRSG;
%save C:/Users/prjohns2.UOFI/Dropbox/MATLAB_F14on/Work_By_Semester/Fall16/SubGradientExps/NewDataHold/DataDistNoAdapt.mat -struct SaveDataStruct
%end

%% Load Data for FuncAdapt plot
%loadData_funcAdapt=1;
%if(loadData_funcAdapt)
%    load C:/Users/prjohns2.UOFI/Dropbox/MATLAB_F14on/Work_By_Semester/Fall16/SubGradientExps/SaveDataHold/DataFunc_Adapt.mat
%end
%%
%semilogy(mycummin(f_ds2sg(1:Tplot))-f_opt)




%% PLot function values Adaptive
plot_functionVals=1;
 Tplot = min([length(f_DSSG) length(f_shor) length(f_rsg)...
     length(f_r2sg) length(f_ds2sg) length(f_sgDecay1) length(f_sgDecay05)]);
if(plot_functionVals)
         fprintf('plotting function values');
	 %Tplot=min([total_sg_iters total_sg_iters2 Total_ad Tshor T_r2sg]);
	 %   Tplot = 7.5e3;
	 %Tplot=7000;
	 f_opt = min(mycummin(f_ds2sg));
	 % addpath(genpath('C:\Users\prjohns2.UOFI\Dropbox\MATLAB_F14on\Work_By_Semester\Fall16\SubGradientExps'));   
	 Grid=1;
	 f_ds2sg_cummin = mycummin(f_ds2sg(1:Tplot));
	 f_DSSG_cummin = mycummin(f_DSSG(1:Tplot));
	 f_shor_cummin = mycummin(f_shor(1:Tplot));
	 f_rsg_cummin = mycummin(f_rsg(1:Tplot));
	 f_r2sg_cummin = mycummin(f_r2sg(1:Tplot));
	 f_sgDecay1_cummin = mycummin(f_sgDecay1(1:Tplot));
	 f_sgDecay05_cummin = mycummin(f_sgDecay05(1:Tplot));
	 start1=1;
	 start2=1;
	 start4=1;

	h=semilogy(...
	start1:Grid:Tplot,f_DSSG_cummin(start1:Grid:Tplot) - f_opt,'-',...
	start2:Grid:Tplot,f_shor(start2:Grid:Tplot)-f_opt,'-',...
	1:Grid:Tplot,f_rsg_cummin(1:Grid:Tplot)-f_opt,'-',...
	start4:Grid:Tplot,f_r2sg_cummin(start4:Grid:Tplot) - f_opt,'-',...
	1:Grid:Tplot,f_ds2sg_cummin(1:Grid:Tplot)-f_opt,'-',...
	1:Grid:Tplot,f_sgDecay1_cummin(1:Grid:Tplot)-f_opt,'-',...
	1:Grid:Tplot,f_sgDecay05_cummin(1:Grid:Tplot)-f_opt,'-',...
	'LineWidth',5);
	%
	%defense legend
	%legend('DS-SG','Yang (2015)','DS2-SG','\alpha_k=0.1 k ^{-1}')
	legend('DS-SG','Shor','RSG','R2SG','DS2-SG','\alpha_k=0.1 k ^{-1}','\alpha_k=0.01 k ^{-0.5}');
	grid on
	set(h,'linewidth',2,'MarkerSize',12)
	xlabel('Number of iterations k');
	ylabel('$h(x_k) - h^*$',...
	'Interpreter', 'Latex', 'FontSize', 15, 'Color', 'k');
	axis([0 Tplot 1e-10 1e2])
end
%% Save Data function values adaptive
SaveData_FuncAdapt=0;
if(SaveData_FuncAdapt)
SaveDataStruct.Tplot=Tplot;
SaveDataStruct.f_DSSG=f_DSSG;
SaveDataStruct.f_shor=f_shor;
SaveDataStruct.f_rsg=f_rsg;
SaveDataStruct.f_r2sg=f_r2sg;
SaveDataStruct.f_ds2sg=f_ds2sg;
SaveDataStruct.c_theta=c_theta;
SaveDataStruct.c_thetaRSG=c_thetaRSG;
SaveDataStruct.c_1=c_1;
SaveDataStruct.c_thetaShor=c_thetaShor;
SaveDataStruct.f_sgDecay1=f_sgDecay1;
SaveDataStruct.f_sgDecay05=f_sgDecay05;
save C:/Users/prjohns2.UOFI/Dropbox/MATLAB_F14on/Work_By_Semester/Fall16/SubGradientExps/NewDataHold/DataFunc_Adapt.mat -struct SaveDataStruct
end

%%
%semilogy(1:Tplot,d_sq(1:Tplot))

%%
%v_in=randn(dim,1);
%[v_out] = L1BallProjection(v_in,tau);
%%

% tau=1; 
% clear t
% start=50;
% change=50;
% endDim=5000;
% dim=start:change:endDim;
% for i=1:length(dim)
%     v_in=randn(dim(i),1);
%    tic
%    [v_out] = L1BallProjection(v_in,tau);    
%    t(i)=toc;
% end
% 
% plot(dim,t)
% 
% %plot(1:dim,v_in,1:dim,v_out)


%% Outdated Code
%% performance dsq
run_outdated=0;
if(run_outdated)
perform_dsqAll=1;
if(perform_dsqAll)
Tplot=min([total_sg_iters total_sg_iters2 T]);

%Tplot=min([T,total_sg_iters]);

semilogy(1:Tplot,d_sq_DS(1:Tplot),1:Tplot,d_sq_rsg(1:Tplot),...1:Tplot,d_sq_3(1:Tplot),...
    1:Tplot,d_sq_decayd1(1:Tplot),1:Tplot,d_sq_decayd05(1:Tplot),...
    1:Tplot,d_sq_shor(1:Tplot))
grid on
axis([0 Tplot 1e-15 1e0])
xlabel('Subgradient Evals');
ylabel('$d(x_k,\mathcal{X})^2$',...
'Interpreter', 'Latex', 'FontSize', 15, 'Color', 'k');
legend('DS-SG','RSG','(A,d)=(0.1,0.99)','(A,d)=(0.01,0.5)','Shors Method')
end

end
%% Load Data for Adapt Plot
loadData_Adapt=0;
if(loadData_Adapt)
   load C:/Users/prjohns2.UOFI/Dropbox/MATLAB_F14on/Work_By_Semester/Fall16/SubGradientExps/SaveDataHold/DataDist_Adapt.mat
end

%% Plot Adapt
%figure;
plotAdapt=0;
if(plotAdapt)
if(loadData_Adapt==0)
   Tplot=min([total_sg_iters total_sg_iters2 Total_ad Tshor T_r2sg]);
   Tplot = 6e3;
end
h=semilogy(...1:Tplot,d_sq_DS(1:Tplot),1:Tplot,d_sq_rsg(1:Tplot)...
        1:Tplot,d_sq_D2DSG(1:Tplot),1:Tplot,d_sq_shor(1:Tplot),...
          1:Tplot,dsq_r2sg(1:Tplot))
set(h,'linewidth',2)
legend('DS-SG', 'RSG','D2S-SG','Shor','R2SG')
xlabel('Subgradient Evals');
ylabel('$d(x_k,\mathcal{X})^2$',...
'Interpreter', 'Latex', 'FontSize', 15, 'Color', 'k');
grid on
end
%% Save Data for Adapt Plot
saveData_Adapt=0;
if(saveData_Adapt)
    SaveDataStruct.Tplot=Tplot;
    SaveDataStruct.d_sq_DS=d_sq_DS;
    SaveDataStruct.d_sq_rsg=d_sq_rsg;
    SaveDataStruct.d_sq_D2DSG=d_sq_D2DSG;
    SaveDataStruct.d_sq_shor=d_sq_shor;
    save C:/Users/prjohns2.UOFI/Dropbox/MATLAB_F14on/Work_By_Semester/Fall16/SubGradientExps/NewDataHold/DataDist_Adapt.mat -struct SaveDataStruct
end
%% sparsity
perform_sparsity=0;
if(perform_sparsity)
Tplot=min([total_sg_iters total_sg_iters2,Total_ad]);
Tplot=5e3;
%figure;
h=plot(1:Tplot,sparsity(1:Tplot),1:Tplot,sparsity_erg(1:Tplot));%,1:Tplot,sparsity_Ad(1:Tplot));
grid on
set(h,'linewidth',3)
ylabel('Number of Nonzeros');
xlabel('Number of subgradient evals');
legend('DS-SG','RSG');
end

%%
beta=2:0.1:10;
C = 10000;
f = sqrt(beta).*log(3*beta).*(C./log(beta)+1);

%plot(beta,f)
