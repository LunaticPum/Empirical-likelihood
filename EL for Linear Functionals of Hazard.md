## 1. Empirical Likelihood, Poisson Version

假设 $X_1,\ldots,X_n$ 为一组独立同分布的随机变量，表示生存分析中的**生存时间**，且这些变量服从连续分布，其分布函数为 $F_0(t)$，相应的累积风险函数为 $\Lambda_0(t)$ 。进一步，假设 $C_1,\ldots,C_n$ 是另一组独立同分布的随机变量，表示生存分析中的**删失时间**，这些变量同样服从连续分布，其分布函数为 $G_0(t)$ ，且删失时间与生存时间相互独立。基于这些假设，可以得到一组带删失的观测量 $(T_i,\delta_i)$ ，其中：
$$T_i=\min(X_i,C_i)\quad\text{and}\quad\delta_i=\mathbb{1}(X_i\leq C_i)\quad \text{for}\,\,\,i=1,2,\ldots,n$$
其中，$T_i$ 是观测到的生存时间，$\delta_i$ 是删失指示变量（$\delta_i = 0$ 表示观测数据为删失）。

对于带删失的观测数据的似然，可通过经验似然来表示为：
$$\prod_{i=1}^nP(T=t_i\,, \delta=\delta_i)=\prod_{i=1}^np_i$$
为了简化推导，令 $$U_1(t)=P(T>t, \delta=1)\,,\,U_0(t)=P(T>t,\delta=0)$$其中，$U_1(t)$ 表示生存时间大于 $t$ 且未删失的概率（时间 $t$ 之后个体仍然存活且观测到时间发生的概率），$U_0(t)$ 表示生存时间大于 $t$ 且删失的概率。于是有
$$-P(T=t_i\,,\,\delta=1)=\Delta U_1(t_i)=U_1(t_i+)-U_1(t_i-)$$
于是，相应的经验似然可写成
$$\prod_{i=1}^n\{-\Delta U_1(t_i)\}^{\delta_i}\{-\Delta U_0(t_i)\}^{1-\delta_i}$$