function mcmc = metropolis(dfunc, chain, theta_init, jumpdist, space, burnin, seed)
% function to apply Metropolis Algorithm and returns a Markov Chain vector
% Input :
%     dfunc : function_handle
%       self-defined density function which take only one parameter (i.e. r.v. X)
%     chain : double, optional
%       the length of Markov Chain to be generated (default 5000)
%     theta_init : double, optional
%       initial value for the chain of θ (default 0.5)
%     jumpdist : use 'makedist' function to specify the jump distribution, optional
%       the distribution that proposed jump (Δθ) follows (default makedist('Normal','mu',0,'sigma',0.1))
%     space : 2D vector ([min, max]), optional
%       the accepted span of the parameter θ (default [-Inf,Inf])
%     burnin : double, optional
%       the length of the chain to be burned from the beginning (default 0)
%     seed : optional
%       random seed (default 'shuffle')
% Output :
%     mcmc : vector
%       generated markov chain vector [burnin + 1 : end]
arguments
    dfunc function_handle
    chain double = 5000;
    theta_init double = 0.5;
    jumpdist = makedist('Normal','mu',0,'sigma',0.1);
    space = [-Inf, Inf];
    burnin double = 0;
    seed = 'shuffle';   
end
if ~isa(dfunc,'function_handle')
    error('dfunc must be a function handle')
end

theta_cur = theta_init;
theta_freq = zeros(chain,1);
theta_freq(1) = theta_init;

rng(seed)

for i = 2:chain
    Delta_theta = random(jumpdist);
    theta_pro = theta_cur + Delta_theta;

    if theta_pro < space(1) || theta_pro > space(2)
        pmoving = 0;
    elseif dfunc(theta_pro) == 0
        pmoving = 1;
    else
        pmoving = min([1,dfunc(theta_pro)/dfunc(theta_cur)]);
    end

    if rand(1) <= pmoving
        theta_freq(i) = theta_pro;
        theta_cur = theta_pro;
    else
        theta_freq(i) = theta_cur;
    end
end

mcmc = theta_freq(burnin + 1 : end);
end