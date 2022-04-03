本程辑包仅为学习和练习之用

This repository is only for study and practice purposes

# Metropolis 算法

假设随机变量 $X$ 服从一个概率密度函数 $P$. 在贝叶斯分析中, 是 $\theta$ 服从一个后验概率分布 $p(\theta|D)$

- 第 0 步: 初始化随机变量 $\theta_{cur}$ (贝叶斯分析). (也可以是其他的随机变量, 如, 随机变量 $X$)
- 第 1 步: 生成一个随机游走的跳跃值 $\Delta \theta \sim N(\mu = 0, \sigma)$ (注: 1. 可以是其他合适的分布, 不一定非得是正态分布. 2. 估计多元分布时, 则从一个多元分布中生成跳跃值, 如, $\Delta \bf{\theta} \sim N(\bf{\mu} = 0, \Sigma)$). 得到提议的下一个随机游走的点 $\theta_{pro} = \theta_{cur} + \Delta \theta$
- 第 2 步: 计算移动到提议的下一个随机游走的点的概率: $p_{move} = \min (1, \dfrac{P(\theta_{pro})}{P(\theta_{cur})})$. 对于贝叶斯分析来说, 由于后验分布概率与似然函数和先验概率之积成比例 $p(\theta|D) \propto p(D|\theta)p(\theta)$, 则有 $p_{move} = \min (1, \dfrac{p(D|\theta_{pro})p(\theta_{pro})}{p(D|\theta_{cur})p(\theta_{cur})})$. (注: 提议的随机游走的点 $theta_{pro}$ 落在可接受范围之外, 则令 $p_{move} = 0$)
- 第 3 步: 检查提议的随机游走是否被接受: 生成一个付出均匀分布 $U(0,1)$ 的随机数. 如果其值小于 $p_{move}$, 则接受随机游走; 否则, 待在原位.
- 第 4 步: 重复第1步至第3步, 直到随机游走生成的马尔可夫链足以表示(representativeness)目标分布 (后验分布)

Check Kruschke's book ("***Doing Bayesian Data Analysis***") for the explanation of MCMC representativeness (**7.5.1**), accuracy (**7.5.2**) and efficiency (**7.5.3**)

# Metropolis Algorithm

Suppose a random variable $X$ follows a probability density (or mass) function $P$.  In Bayesian Analysis, it is $\theta$ follows the posterior distribution function $p(\theta|D)$

- step 0: Initialize the random variable $\theta_{cur}$ (posterior case). (For the other cases, it could be any random variable, e.g. r.v. $X$)
- step 1: Generate a random proposed jump $\Delta \theta \sim N(\mu = 0, \sigma)$ (Note: 1. it could be other kind suitable distributions. 2. for the multivariate distributions, generate the jump from a multivariate distribution, e.g. $\Delta \bf{\theta} \sim N(\bf{\mu} = 0, \Sigma)$). $\theta_{pro} = \theta_{cur} + \Delta \theta$
- step 2: Computing the moving probability for the random walk: $p_{move} = \min (1, \dfrac{P(\theta_{pro})}{P(\theta_{cur})})$. For the posterior distribution, since it is proportional to the product of prior and likelihood $p(\theta|D) \propto p(D|\theta)p(\theta)$, we have $p_{move} = \min (1, \dfrac{p(D|\theta_{pro})p(\theta_{pro})}{p(D|\theta_{cur})p(\theta_{cur})})$. (Note: if the proposed r.v. $theta$ is outside of its acceptable span, then we set $p_{move}$ to $0$)
- step 3: Check if the random walk jump is accepted: generate a random number from $U(0,1)$ uniform distribution. If it is less than $p_{move}$, then we accept the jump, otherwise, we stay.
- step 4: Repeating step 1 to step 3, until the random walks from above generate a markov chain which is representative of our target distribution (posterior distribution)

Check Kruschke's book ("***Doing Bayesian Data Analysis***") for the explanation of MCMC representativeness (**7.5.1**), accuracy (**7.5.2**) and efficiency (**7.5.3**)

# 内容
本程辑包提供了R、Julia、Matlab以及Python的代码，除`.R`外，其他文件中均有详细的帮助文档，可有效得帮助理解各个代码

# Contents
This repository provides R, Julia, Matlab and Python codes. Full documentation can be found in all but `.R` file, which hopefully would help users better understand the codes.