from typing import Callable
import autolik
import numpy as np
from numpy import _ArrayFloat_co
import scipy.stats

class Proposed:
    # Option = ["normal", "hamiltonian"]

    def __init__(self, dim, seed):
        self.dim = dim
        self.seed = seed

    def jump(self):
        rng = np.random.default_rng(self.seed)
        if self.dim == 1:
            return rng.normal(0,1,1)
        else:
            return rng.multivariate_normal(np.zeros(self.dim),np.ones((self.dim,self.dim)),1)


def HMC(U:Callable, eps, L, current_q:_ArrayFloat_co, seed):
    q = current_q
    p = Proposed(len(q), seed).jump()
    current_p = p

    # make a half step for momentum at beginning
    p -= eps * autolik.grad(U)(q) / 2

    # alternate full steps for position and momentum

    for _ in range(L-1):
        # make a full step for the position
        q += eps * p
        # make a full step for the momentum, except the end of the trajectory
        p -= eps * autolik.grad(U)(q)

    # make a full step for the position at the end of the trajectory
    q += eps * p
    # make a half step for the momentum at the end of the trajectory
    p -= eps * autolik.grad(U)(q) / 2

    p = -p

    # evaluate potiential and kinetic energies at the start and the end of trajectory
    current_U = U(current_q)
    current_K = sum(current_p**2) / 2
    proposed_U = U(q)
    proposed_K = sum(p**2) / 2

    # accept or reject the state at end of trajectory, returning either the position at the end of the trajectory or the initial position

    if np.random.default_rng(seed).uniform() < np.exp(current_U - proposed_U + current_K - proposed_K):
        return q # accept
    else:
        return current_q

