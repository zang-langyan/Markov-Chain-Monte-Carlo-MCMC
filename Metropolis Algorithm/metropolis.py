import numpy as np
from scipy.special import beta
import scipy.stats

class MCMC:
    """Markov Chain Monte Carlo method class used by MCMC Algorithms functions (e.g. metropolis)

    Used to initialize instances for MCMC Algorithm functions

    Return an MCMC object

    Parameters
    ----------
    dfunc : function
        self-defined density function which take only one parameter (i.e. r.v. X)
    chain : int, optional
        the length of Markov Chain to be generated (default 5000)
    theta_init : float, optional
        initial value for the chain of θ (default 0.5)
    jumpdist : scipy.stats.rv_continuous or scipy.stats.rv_discrete, optional
        the distribution that proposed jump (Δθ) follows (default scipy.stats.norm(0,0.2))
    space : list, optional
        the accepted span of the parameter θ (default [-np.inf,np.inf])
    burnin : int, optional
        the length of the chain to be burned from the beginning (default 0)
    seed : int or None, optional
        random seed (default None)

    Usage:
    -----
    >>> import metropolis
    >>> density = lambda x: scipy.stats.gamma(2,loc=4,scale=5).pdf(x)
    >>> d = metropolis.MCMC(density, chain = 10000, jumpdist=scipy.stats.norm(loc=0,scale=2), space = [0,np.inf])
    >>> d
    <metropolis.MCMC object ...>
    """
    def __init__(self, 
                 dfunc, 
                 chain = 5000, 
                 theta_init = 0.5, 
                 jumpdist = scipy.stats.norm(0,0.2), 
                 space = [-np.inf,np.inf], 
                 burnin = 0, 
                 seed = None) -> None:
        self.dfunc = dfunc
        self.chain = chain
        self.theta_init = theta_init
        self.jump = jumpdist
        self.space = space
        self.burnin = burnin
        self.seed = seed

    def metropolis(self, *arg, **kwarg):
        """
        function applying Metropolis Algorithm to MCMC objects

        Parameters
        ----------
        dfunc : function, *arg, optional (update)
            self-defined density function which take only one parameter (i.e. r.v. X)
        chain : int, optional (update)
            the length of Markov Chain to be generated
        theta_init : float, optional (update)
            initial value for the chain of θ (default 0.5)
        jumpdist : scipy.stats.rv_continuous or scipy.stats.rv_discrete, optional (update)
            the distribution that proposed jump (Δθ) follows (default scipy.stats.norm(0,0.2))
        space : list, optional (update)
            the accepted span of the parameter θ (default [-np.inf,np.inf])
        burnin : int, optional (update)
            the length of the chain to be burned from the beginning (default 0)
        seed : int or None, optional (update)
            random seed (default None)

        Examples
        --------
        >>> import metropolis
        >>> from scipy.special import beta
        >>> density = lambda x: x**14 * (1-x)**6 / beta(15,7) # Beta(15,7) distribution density function
        >>> d = metropolis.MCMC(density, space = [0,1], burnin = 5, seed = 72)
        >>> result = d.metropolis(chain = 50, theta_init = 0.1)
        >>> result
        [0.8176960093922774, 0.6965658096789994, 0.7918980615882665, 0.7742394795862185, 0.7742394795862185, 0.6215414421599866, 0.6215414421599866, 0.6468248283946754, 0.7523307146724492, 0.7523307146724492, 0.7680822171469903, 0.6717376155310788, 0.6717376155310788, 0.8189878864825765, 0.7533257832372013, 0.7648206889254298, 0.7648206889254298, 0.7648206889254298, 0.7648206889254298, 0.7648206889254298, 0.7967947965010628, 0.707139903476412, 0.707139903476412, 0.707139903476412, 0.707139903476412, 0.7949590848197712, 0.48018105229278457, 0.48018105229278457, 0.6163978253871556, 0.6163978253871556, 0.6163978253871556, 0.6241841537823003, 0.6241841537823003, 0.6241841537823003, 0.6241841537823003, 0.6922332333961135, 0.8204027203103916, 0.7521267229207833, 0.7521267229207833, 0.808566620237602, 0.808566620237602, 0.6460288022318678, 0.5452360032045938, 0.5452360032045938, 0.5934357613729606]
        """
        # update attributes
        for args in arg:
            self.dfunc = args
        
        for kw, value in kwarg.items():
            if kw == 'chain':
                self.chain = value
            elif kw == 'theta_init':
                self.theta_init = value
            elif kw == 'jumpdist':
                self.jump = value
            elif kw == 'space':
                self.space = value
            elif kw == 'burnin':
                self.burnin = value
            elif kw == 'seed':
                self.seed = value
            else:
                raise Exception(f'keyword argument "{kw}" not supported',)

        # check if dfunc callable
        if not callable(self.dfunc):
            raise Exception("dfunc must be a function. recreate the object with a valid density function")
        
        # Metropolis Algorithm
        theta_cur = self.theta_init
        theta_freq = [self.theta_init]
        
        rng = np.random.default_rng(self.seed)

        while True:
            Delta_theta = self.jump.rvs(random_state=rng)
            theta_pro = theta_cur + Delta_theta

            if theta_pro < self.space[0] or theta_pro > self.space[1]:
                pmoving = 0
            elif self.dfunc(theta_cur) == 0:
                pmoving = 1
            else:
                pmoving = min(1,self.dfunc(theta_pro)/self.dfunc(theta_cur))
            
            # np.random.rand()
            if scipy.stats.uniform().rvs(random_state=rng) <= pmoving:
                theta_freq.append(theta_pro)
                theta_cur = theta_pro
            else:
                theta_freq.append(theta_cur)

            if len(theta_freq) >= self.chain:
                break

        return theta_freq[self.burnin:]

def main():
    # density = lambda x: x**14 * (1-x)**6 / beta(15,7)
    density = lambda x: scipy.stats.gamma(2,loc=4,scale=5).pdf(x)
    d = MCMC(density, chain = 10, jumpdist=scipy.stats.norm(loc=0,scale=2), space = [0,np.inf])
    print(d)
    result = d.metropolis(chain=100, seed = 42)
    print(result)
    
if __name__ == "__main__":
    main()