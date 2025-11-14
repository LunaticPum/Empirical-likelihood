# 1. 右删失数据的经验似然
假设 $X_1,\ldots,X_n$ 为一组独立同分布的随机变量，表示生存分析中的**生存时间**，且这些变量服从连续分布，其累积分布函数为 $F_0(t)$，相应的*累积风险函数*为 $\Lambda_0(t)$ 。进一步，假设 $C_1,\ldots,C_n$ 是另一组独立同分布的随机变量，表示生存分析中的**删失时间**，这些变量同样服从连续分布，其累积分布函数为 $G_0(t)$ ，且删失时间与生存时间相互独立。基于这些假设，可以得到一组带删失的观测量 $(T_i,\delta_i)$ ，其中：
$$T_i=\min(X_i,C_i)\quad\text{and}\quad\delta_i=\mathbb{1}(X_i\leq C_i)\quad \text{for}\,\,\,i=1,2,\ldots,n$$
其中，$T_i$ 是观测到的生存时间，$\delta_i$ 是删失指示变量（$\delta_i = 0$ 表示观测数据为删失）。

对于带删失的观测数据，其似然函数可通过经验似然来表示为：
$$\prod_{i=1}^nP(T=t_i\,, \delta=\delta_i)=\prod_{i=1}^np_i$$
为了简化推导，令 ：$$U_1(t)=P(T>t, \delta=1)\,,\,U_0(t)=P(T>t,\delta=0)$$其中，$U_1(t)$ 表示生存时间大于 $t$ 且未删失的概率（即在时间 $t$ 之后个体仍然存活且观测到时间发生的概率），$U_0(t)$ 表示生存时间大于 $t$ 且删失的概率。由此，有：
$$-P(T=t_i\,,\,\delta=1)=\Delta U_1(t_i)=U_1(t_i+)-U_1(t_i-)$$
因此，相应的经验似然函数可表示为：
$$\prod_{i=1}^n\{-\Delta U_1(t_i)\}^{\delta_i}\{-\Delta U_0(t_i)\}^{1-\delta_i}$$

接下来，我们引入生存时间和删失时间的累积分布函数，可以得到：
$$\begin{align}
	U_1(t)&=P(T>t\,,\,\delta=1)=\int_t^\infty1-G(s)\,\mathrm{d}F(s) \\
	U_0(t)&=P(T>t\,,\,\delta=0)=\int_t^\infty1-F(s)\,\mathrm{d}G(s)
\end{align}$$
因此，经验似然函数可以写作：
$$\prod_{i=1}^n\{(1-G(t_i))\,\mathrm{d}F(t_i)\}^{\delta_i}\{(1-F(t_i))\,\mathrm{d}G(t_i)\}^{1-\delta_i}$$
在大部分生存分析中，我们通常只关注与生存时间相关的分布 $F$ 。因此，与 $F$ 无关的部分可以视作常量，从而将与 $F$ 相关的经验似然函数简化为：
$$EL(F)\propto\prod_{i=1}^n\{\mathrm{d}F(t_i)\}^{\delta_i}\{(1-F(t_i))\}^{1-\delta_i}$$
其中，$\propto$ 表示成比例关系。相应的对数表达式为：
$$\log EL(F)=\sum_{i=1}^n\delta_i\log\mathrm{d}F(t_i)+(1-\delta_i)\log\left[1-F(t_i)\right]$$
### 1.1 风险函数的引入
我们可以观察到，对于病人个体，其生存时间的累积风险函数 $\Lambda(t)$ 和累积分布函数 $F(t)$ 存在一一对应的关系，因此，能将与 $F$ 相关的经验似然表示为与 $\Lambda$ 相关的经验似然，为此，我们考虑以下两种经验似然的表示形式：
#### 泊松似然
由于生存函数 $S(t) = 1 - F(t)$ 与累积风险函数 $\Lambda(t)$ 之间的关系式为
$$1-F(t)=\exp(-\Lambda(t))$$
因此，与 $\Lambda$ 相关的对数经验似然表示为：
$$\begin{align}
&\quad\,\,\delta_i\log\mathrm{d}F(t_i)+(1-\delta_i)\log\left[1-F(t_i)\right]\\
&= \delta_i\log\mathrm{d}[1-\exp(-\Lambda(t_i))]+(1-\delta_i)(-\Lambda(t_i))\\
&=  \delta_i\log[\exp(-\Lambda(t_i))\mathrm{d}\Lambda(t_i)]+(1-\delta_i)(-\Lambda(t_i))\\
&= \delta_i[-\Lambda(t_i)+\Lambda(t_i)+\log\mathrm{d}\Lambda(t_i)]-\Lambda(t_i)\\
&= \delta_i\log\Delta\Lambda(t_i)-\Lambda(t_i)
\end{align}$$
其中，$\Delta\Lambda(t_i)=\Lambda(t_i+)-\Lambda(t_i-)$ 表示累积风险函数 $\Lambda(t_i)$ 在时间点 $t_i$ 的”**跳跃**“。

然而，在实际使用中，通过 Kaplan-Meier 估计等方法得到的生存相关的累积分布函数 $F(t)$ 的估计是离散的。因此，累积风险函数 $\Lambda(t)$ 也应当是离散的函数形式。基于此，我们不能简单利用连续形式的 $\Lambda(t)$ 来描述生存函数的变化，而必须考虑到 $\Lambda(t)$ 的**跳跃特性**，即 $\Lambda(t)$ 是一个阶梯函数。*在每个事件发生时，累积风险函数会发生跳跃*，且 $\Lambda(t)$ 在某一时刻 $t_i$ 的值等于在 $t_i$ 之前所有事件跳跃的总和，具体表示为：
$$\Lambda(t_i)=\sum_j\Delta\Lambda(t_j)\mathbb{1}[t_j\leq t_i]$$
其中 $\Delta\Lambda(t_j)$ 是在 $t_j$ 时刻的跳跃大小，且只有在事件发生（即未删失）时，累积风险函数才会发生跳跃，即 $\Delta\Lambda(t_j)\neq 0$ 。因此，与 $\Lambda(t)$ 相关的对数经验似然表示为：
$$\log EL(\Lambda)=\sum_{i=1}^n\left(\delta_i\log\Delta\Lambda(t_i)-\sum_j\Delta\Lambda(t_j)\mathbb{1}[t_j\leq t_i]\right)$$

### 1.2 经验似然估计
根据上述过程，我们能得到删失数据的关于累积风险函数的经验似然有如下两个形式：
$$\begin{align}
	EL(\Lambda)&=\prod_{i=1}^n\left\{\Delta\Lambda(t_i)\right\}^{\delta_i}\exp\left\{\sum_j\Delta\Lambda(t_j)\mathbb{1}[t_j\leq t_i]\right\} \\
	
\end{align}$$