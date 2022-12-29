# meth-ising

*Work in progress*

Fitting 5-parameter 1D Ising models for regionwise information-theoretic properties of DNA methylation states. Implemented in Rust.

## Theory

*Notes summarizing Supplementary Theory/Materials (from Jenkinson et al., 2017 & 2018) for 5-parameter 1D Ising model.*

### 1. 1D Ising model for DNA methylation

We can define an energy function for the observed methylation pattern $x$. $x_n$'s denote the methylation state (0 or 1) of $n$th CpG in methylation pattern $x$.

$$U(x) = -\sum_{n=1}^{N}a_n(2x_n-1)-\sum_{n=2}^{N}c_n(2x_n-1)(2x_{n-1}-1)$$

The first term evaluates the average methylation level, and the second term reflects the local correlation of DNA methylation states. Also note that $2x_n-1$ is 1 when $x_n$ is methylated, otherwise -1.

This formulation requires too many parameters; $a_n$ for all $N$ CpGs, and $c_n$ for all $N-1$ adjacent CpG pairs (total 2 * $N$-1). To address this issue, the authors parameterize $a_n$ and $c_n$ with region-specific parameters $\alpha$, $\beta$ and $\gamma$.

$$a_n=\alpha + \beta\rho_n$$
$$c_n=\frac{\gamma}{d_n}$$

where $\rho_n$ is the density (0~1) of CpG sites within $\pm500$ nt window centered at $n$ th CpG, and $d_n$ is the distance (in nt) between $n$ th and $n-1$ th CpG.

### 2. Computing the probability of methylation pattern $x$ given 5 parameters $\alpha', \alpha, \alpha'', \beta, \gamma$


All the followings are applied for each genomic region $\mathcal{R_k}$, which comprises $R$ CpG sites 1, 2, 3, ..., $R$.

By assuming that the methylation states of the first (1) and the last ( $R$ ) CpGs are independent, the probability of observing a methylation pattern in region $\mathcal{R_k}$ is as follows:

$$P_X(x_1, x_2, ..., x_R) = P_X(x_2, x_3, ..., x_{R-1})P_x(x_1, x_R) \\ \simeq P_x(x_2, x_3, ..., x_{R-1})P_x(x_1)P_x(x_R)$$

Besides, defining energy $U(x)$ for pattern $x$ (See 1.) allows us to compute the probability of observing $x$ using the Boltzmann-Gibbs distribution.

$$P(x) = \frac{1}{Z}\exp\{-U(x)\}$$

Therefore,

$$P_x(x_2,x_3,...,x_{R-1}|x_1,x_R) \propto \exp \{ \sum_{r=2}^{R-1} (\alpha_k + \beta_k\rho_r)(2x_r-1) + \sum_{r=2}^{R}\frac{\gamma_k}{d_r}(2x_r-1)(2x_{r-1} - 1) \} $$

And if we set

$$\alpha'_k=\frac{1}{2}\ln\frac{Pr[X_1=1]}{1-Pr[X_1=1]}$$

$$\alpha''_k=\frac{1}{2}\ln\frac{Pr[X_R=1]}{1-Pr[X_R=1]}$$

We obtain

$$P_x(x_1,x_2,...,x_{R-1},X_R) = \frac{1}{Z} \exp \{ \alpha'_k(2x_1-1) + \sum_{r=2}^{R-1} (\alpha_k + \beta_k\rho_r)(2x_r-1) + \alpha''_k(2x_R-1) + \sum_{r=2}^{R}\frac{\gamma_k}{d_r}(2x_r-1)(2x_{r-1} - 1) \} $$

Where the partition function $Z$ is defined as the sum of $P$'s for every possible pattern $u$.

$$Z=\sum_{u}\exp \{ \alpha'_k(2u_1-1) + \sum_{r=2}^{R-1}(\alpha_k+\beta_k\rho_r)(2u_r-1) + \alpha''_k(2u_R-1) + \sum_{r=2}^{R}\frac{\gamma_k}{d_r}(2u_r-1)(2u_{r-1}-1) \}$$

### 3. Computing partition function $Z$

The number of every possible DNA methylation pattern consisting of $n$ CpG sites are $2^n$, which makes naive computation of $Z$ infeasible even with the moderate number of CpGs ($n$). But the computation of $Z$ can be efficiently implemented with dynamic programming. Details are below:

We should first reformulate $P_x$ as multiplicative form

$$P_x(x_1,x_2,...,x_R) = \frac{1}{Z} \prod_{r=1}^{R-1}\phi_r(x_r, x_{r+1})$$

where

$$\phi_1(x_1, x_2) = \exp \{\alpha'_k(2x_1-1) + (\alpha_k+\beta_k\rho_2)(2x_2 - 1) + \frac{\gamma_k}{d_2}(2x_1-1)(2x_2-1) \}$$

$$\phi_1(x_r, x_{r+1}) = \exp \{(\alpha_k+\beta_k\rho_{r+1})(2x_r - 1) + \frac{\gamma_k}{d_{r+1}}(2x_r-1)(2x_{r+1}-1) \}$$

$$\phi_{R-1}(x_{R-1}, x_R) = \exp \{\alpha''_R(2x_{R-1}-1) + \frac{\gamma_k}{d_R}(2x_{R-1}-1)(2x_R-1) \}$$

Then the partition function $Z$ is

$$Z= \sum_{u_1 = 0}^1 \sum_{u_2 = 0}^1 ... \sum_{u_R = 0}^1 \prod_{r=1}^{R-1} \phi_r (u_r, u_{r+1}) $$

which can be computed with the following DP formulation

$$Z_R(x) = 1 \enspace \text{for} \enspace x = 0,1$$

$$Z_r(x) = \phi_r(x, 0)Z_{r+1}(0) + \phi_r(x, 1)Z_{r+1}(1), \enspace \text{for} \enspace x=0,1 \enspace r=R-1, R-2, ..., 1$$

$$Z = Z_1(0) + Z_1(1)$$

### 4. Marginal PMFs

To optimize $\alpha', \alpha, \alpha'', \beta$ and $\gamma$ using observed methylation patterns from bisulfite sequencing data, we should compute marginal PMFs since in most cases we cannot obtain full observations for $x_n$'s due to the short length of bisulfite sequencing reads (mostly <300bp, where the size of $\mathcal{R}_k$ is 3,000bp). The following shows how to compute marginal PMFs given partial observation of methylation pattern $x_{q:q+s}$ and parameters $\alpha', \alpha, \alpha'', \beta$ and $\gamma$.

In this case, we reformulate $P_x$ as a form of inhomogeneous Markov chain.

$$P_x(x_1, x_2, ..., x_R) = z_1(x_1) \prod_{r=1}^{R-1}z_{r+1}(x_{r+1}|x_r)$$

where the initial probability is given as

$$z_1(x_1) = \frac{Z_1(x_1)}{Z}$$

with the transition probability

$$z_{r+1}(x_{r+1}|x+r) = \frac{\phi_r(x_r, x_{r+1})Z_{r+1}(x_{r+1})}{Z_r(x_r)} \enspace \text{for} \enspace r=1,2...R-1$$

We can now compute probabilities for partial observations

$$P_x(x_{1:q+s}) = z_1(x_1)\prod_{r=1}^{q+s-1} z_{r+1} (x_{r+1}|x_r)$$

$$= z_1(x_1) \prod_{r=1}^{q-1} z_{r+1} (x_{r+1}|x_r) \prod_{r=1}^{q+s-1} z_{r+1} (x_{r+1}|x_r)$$

By marginalizing out for 1~$q$ th CpGs for all possible methylation patterns consisting of 1~$q$ th CpGs, we have

$$P_x(x_{q:q+s}) = w_q(x_q) \prod_{r=q}^{q+s-1} z_{r+1} (x_{r+1}|x_r)$$

where

$$w_q(x_q) = \sum_{x_1} \sum_{x_2} ... \sum_{x_{q-1}} z_1(x_1) \prod_{r=1}^{q-1} z_{r+1} (x_{r+1}|x_r)$$



## See also

- Jenkinson, G., Pujadas, E., Goutsias, J., & Feinberg, A. P. (2017). Potential energy landscapes identify the information-theoretic nature of the epigenome. Nature genetics, 49(5), 719-729.

- Jenkinson, G., Abante, J., Feinberg, A. P., & Goutsias, J. (2018). An information-theoretic approach to the modeling and analysis of whole-genome bisulfite sequencing data. BMC bioinformatics, 19(1), 1-23.

- Jenkinson, G., Abante, J., Koldobskiy, M. A., Feinberg, A. P., & Goutsias, J. (2019). Ranking genomic features using an information-theoretic measure of epigenetic discordance. BMC bioinformatics, 20(1), 1-17.

- InformME (https://github.com/GarrettJenkinson/informME)