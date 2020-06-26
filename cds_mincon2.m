%%
%
% DIFFERENTIAL SEARCH ALGORITHM (DSA) (in MATLAB) for Constrained
% optimization  % STANDARD VERSION of DSA (16.July.2013)
%
% usage : > ds(method,fnc,mydata,popsize,dim,low,up,maxcycle)
% ds(@Prg1f,@Prg1c,4,10,13,zeros(1,13),(1,1,бн,1,100,100,100,1),500)
% method
%--------------------------------------------------------------
%--------------------------------------------------------------
% Please cite this article as;
% P.Civicioglu, "Transforming geocentric cartesian coordinates to geodetic coordinates by using differential search algorithm",  Computers & Geosciences, 46 (2012), 229-247.
% P.Civicioglu, "Understanding the nature of evolutionary search algorithms", Additional technical report for the project of 110Y309-Tubitak,2013.
%
% (Basic/simple version can be found in matlab file exchange)
%%--------------------------------------------------------------

function [g_minimizer,g_min,y_con]=cds_mincon2(fhandle,nonhandle,method,N_so,Dim,low_lmt,up_lmt,max_epoch,method_init)
 global epk fth lambda
 global pfactor1 t1 t2
 fth 
 pfactor1=0;
% method=4;
% N_so=40;% size
% max_epoch=6000;
% method_init='rnd';
num_feval=0;% function evaluating number
% N_so ; size of population.
% Dim ; size of problem dimension (1,2,3,...,d), where each clan (i.e. sub-so) includes d-individuals.
%Initialization
% control of habitat lmts

if numel(low_lmt)==1,
    low_lmt=low_lmt*ones(1,Dim);
    up_lmt=up_lmt*ones(1,Dim);
end

%% generate initial individuals, clans and so.
switch method_init
    case 'rnd'
        so=genpop_rnd(N_so,Dim,low_lmt,up_lmt);
    case 'kent'
        so=genpop_kent(N_so,Dim,low_lmt,up_lmt);
        % success of clans/so
        % fit_so=feval(fnc,so,mydata);
end

for i=1:N_so
    zn(i)=Fun(fhandle,nonhandle,so(i,:));
    fit_so(i)=zn(i);
end
num_feval=N_so; f_draw=inf*ones(1,max_epoch);fyvary=f_draw;H_b_x=f_draw;
%%%%fit_so=feval(Fun,so);
kk=1;
for epk=1:max_epoch
    % SETTING OF ALGORITHMIC CONTROL PARAMETERS
    % Trial-pattern generation strategy for morphogenesis; 'one-by-one morphogenesis'.
    % p1=0.0*rand;  % i.e.,  0.0 <= p1 <= 0.0
    % p2=0.0*rand;  % i.e.,  0.0 <= p2 <= 0.0
    
    % Trial-pattern generation strategy for morphogenesis; 'one-or-more morphogenesis'. (DEFAULT)
    p1=0.3*rand;  % i.e.,  0.0 <= p1 <= 0.3
    p2=0.3*rand;  % i.e.,  0.0 <= p2 <= 0.3
    %-------------------------------------------------------------------
    [direction,msg]=generate_direction(method(randi(numel(method))),so,N_so,fit_so);
    map=generate_map_of_active_individuals(N_so,Dim,p1,p2);
    %-------------------------------------------------------------------
    % Recommended Methods for generation of Scale-Factor; R
    % R=4*randn;  % brownian walk
    % R=4*randg;  % brownian walk
    % R=lognrnd(rand,5*rand);  % brownian walk
    R=1./gamrnd(1,0.5);   % pseudo-stable walk
    % R=1/normrnd(0,5);    % pseudo-stable walk
    %-------------------------------------------------------------------
    % bio-interaction (morphogenesis)
    stopover=so+(R.*map).*(direction-so);
    % Boundary Control
    stopover=update(stopover,low_lmt,up_lmt);
    % Selection-II
    %     fit_stopover=feval(fnc,stopover,mydata);
    fit_stopover=100*ones(1,N_so);
    for i=1:N_so
       fit_stopover(i)=Fun(fhandle,nonhandle,stopover(i,:));
    end
    num_feval=num_feval+N_so;% fun eval number
    ind=fit_stopover<fit_so;
    fit_so(ind)=fit_stopover(ind);
    so(ind,:)=stopover(ind,:);
    
    
    % update results
    [g_min,indexbest]=min(fit_so);
    g_minimizer=so(indexbest,:);
    f_min=fhandle(g_minimizer);
    % export results
    assignin('base','g_min',g_min);
    assignin('base','g_minimizer',g_minimizer);
    if epk==kk*100
        fprintf('%s  | %5.0f ---> %10.10f   %10.10f\n',msg,epk,f_min,g_min)
        kk=kk+1;
    end    
end

g_minimizer
y_con=feval(nonhandle,g_minimizer)

%% 
function pop=genpop_rnd(N,D,low,up)
pop=ones(N,D);
for i=1:N
    for j=1:D
        pop(i,j)=rand*(up(j)-low(j))+low(j);
    end
end
%% Kent map for initial population
function pop=genpop_kent(N,D,low,up)
a=0.4; 
for k=1:N
    x(k,:)=rand(1,D).*(up-low)+low;
end
for i=1:N
    for j=1:D
        if x(i,j)<=a;
            x(i,j)=x(i,j)/a;
        else
            x(i,j)=(1-x(i,j))/(1-a);
        end
    end
end
% here we may not set the x in [low, up], because this a constrained optimization. 
pop=x;


function p=update(p,low,up)
[popsize,dim]=size(p);
for i=1:popsize
    for j=1:dim
        % first (standard)-method
        if p(i,j)<low(j), if rand<rand, p(i,j)=rand*(up(j)-low(j))+low(j); else p(i,j)=low(j); end, end
        if p(i,j)>up(j),  if rand<rand, p(i,j)=rand*(up(j)-low(j))+low(j); else p(i,j)=up(j); end, end
        
        %{
       %  second-method
        if rand<rand,
            if p(i,j)<low(j) || p(i,j)>up(j), p(i,j)=rand*(up(j)-low(j))+low(j); end
        else
            if p(i,j)<low(j), p(i,j)=low(j); end
            if p(i,j)>up(j),  p(i,j)=up(j); end
        end
        %}
    end
end
function [direction,msg]=generate_direction(method,so,N_so,fit_so)
switch method
    case 1,
        % BIJECTIVE DSA  (B-DSA) (i.e., go-to-rnd DSA);
        % philophy: evolve the so (i.e.,population) towards to "permuted-so (i.e., random directions)"
        direction=so(randperm(N_so),:); msg=' B-DSA';
    case 2,
        % SURJECTIVE DSA (S-DSA) (i.e., go-to-good DSA)
        % philophy: evolve the so (i.e.,population) towards to "me of the random top-best" lutions
        ind=ones(N_so,1);
        [null,B]=sort(fit_so);
        for i=1:N_so, ind(i)=B(randi(ceil(rand*N_so),1)); end;
        direction=so(ind,:);  msg=' S-DSA';
    case 3,
        % ELITIST DSA #1 (E1-DSA) (i.e., go-to-best DSA)
        % philophy: evolve the so (i.e.,population) towards to "one of the random top-best" lution
        [null,jind]=sort(fit_so); ibest=jind(ceil(rand*N_so)); msg='E1-DSA';
        direction=repmat(so(ibest,:),[N_so 1]);
    case 4,
        % ELITIST DSA #2 (E2-DSA) (i.e., go-to-best DSA)
        % philophy: evolve the so (i.e.,population) towards to "the best" lution
        [null,ibest]=min(fit_so); msg='E2-DSA';
        direction=repmat(so(ibest,:),[N_so 1]);
end
return
function map=generate_map_of_active_individuals(N_so,Dim,p1,p2)
%% strategy-selection of active/passive individuals
map=zeros(N_so,Dim);
if rand<rand,
    if rand<p1,
        % Random-mutation #1 strategy
        for i=1:N_so
            map(i,:)=rand(1,Dim) < rand;
        end
    else
        % Differential-mutation strategy
        for i=1:N_so
            map(i,randi(Dim))=1;
        end
    end
else
    % Random-mutation #2 strategy
    for i=1:N_so
        map(i,randi(Dim,1,ceil(p2*Dim)))=0;
    end
end
return
function [ns]=findlmts(n,ns,Lb,Ub)
%% Make sure the fireflies are within the bounds/lmts
for i=1:n,
    % Apply the lower bound
    ns_tmp=ns(i,:);
    I=ns_tmp<Lb;
    ns_tmp(I)=Lb(I);
    
    % Apply the upper bounds
    J=ns_tmp>Ub;
    ns_tmp(J)=Ub(J);
    % Update this new move
    ns(i,:)=ns_tmp;
end

% -----------------------------------------
% d-dimensional objective function
function z=Fun(fhandle,nonhandle,u)
%% Objective
z=fhandle(u);
% Apply nonlinear constraints by the penalty method
% Z=f+sum_k=1^N lam_k g_k^2 *H(g_k) where lam_k >> 1
z=z+getnonlinear(nonhandle,u);

%% %% Apply Dynamic Penalty for inequality constraints as a penalty function
function [Z,H_b_x]=getnonlinear(nonhandle,u)
Z=0;H_b_x=0;beta=2;epsilon=10^(-3);
global pfactor1
%% %%%%global f_init_max
% Get nonlinear constraints
g = nonhandle(u);
for k=1:length(g)
    pf=pfactor1+DPFactor(g(k),epsilon);
    Z=Z+pf*g(k)^beta*getH(g(k),epsilon);
    H_b_x=H_b_x+g(k)^beta*getH(g(k),epsilon);
end

%% get dynamic penalty factor
function pf=DPFactor(gk,epsilon)
global epk t1 t2
global max_epoch
%% % only dely on gen
%  t1=5;t2=8;
%  pf=10.^((t2-t1)./(1+exp(20*(-epk+max_epoch/5)/max_epoch))+t1);
% Test if inequalities hold
% H(g) which is mething like an index function
%% % p dely on fval and gval at point x_k
% gv=g(find(g>epsilon));
% vmax=max(gv);
vmax=gk;
%t1=3;t2=6;% this setting is better for most of 30D 
%t1=3;t2=6;

% pf=10^t1/(1+exp(t1*(-epk+max_epoch/t2)/max_epoch));
if vmax>epsilon
   pf=10^t1/(1+exp(t1*(-epk+max_epoch/t2)/max_epoch));
%    pf=(1+exp(t1*(-epk+max_epoch/2)/max_epoch));
else
    pf=epsilon;
end
% t1=4;t2=10; % this setting is better for 30D
% if vmax>10^-5 
%     %pf=10^t1*sum(gv)/(1+exp(t1*(-epk+max_epoch/2*t2)/max_epoch));
%     pf=10^t1*vmax/(1+exp(t1*(-epk+max_epoch/2*t2)/max_epoch));
% 
% else
%     pf=0;
% end

function H=getH(g,epsilon)
if g<=epsilon,
    H=0;
else
    H=1;
end

% Test if equalities hold
function H=geteqH(g)
if g==0,
    H=0;
else
    H=1;
end
%%
%% ==== End of Algorithm implementation ======