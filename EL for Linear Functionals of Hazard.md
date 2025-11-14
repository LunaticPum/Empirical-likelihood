# 1. 右删失数据的经验似然估计
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
根据上述过程，我们可以得到删失数据中关于累积风险函数的经验似然的两种表示形式：
$$\begin{align}
	EL(\Lambda)&=\prod_{i=1}^n\left\{\Delta\Lambda(T_i)\right\}^{\delta_i}\exp\left\{-\sum_j\Delta\Lambda(T_j)\mathbb{1}[T_j\leq T_i]\right\} \\
	EL(\Lambda)&=\prod_{i=1}^n\left\{\Delta\Lambda(T_i)\right\}^{\delta_i}\exp\left\{-\Lambda(T_i)\right\} 
\end{align}$$
其中，第一行表达式是基于假设累积风险函数 $\Lambda(T)$ 为纯离散函数（阶梯函数）的情况下得到的。

>[!tip|pz]+ **最大观测值删失对风险函数估计的影响**
> 当最大观测值被删失时，风险函数的经验似然估计（NPMLE）和Kaplan-Meier估计都会面临类似的问题：*对于删失后的时间段，无法给出唯一的估计值*。这是因为删失数据提供了关于该时间点之后的信息的不完全性，导致无法确定删失时间点之后的生存风险或生存概率。

由于数据中存在删失，意味着我们*对病人的生存时间进行的是不完全观测*，因此生存时间的观测是离散的，而其相应的累积风险函数也必然具有**跳跃性**。因此，我们只能选择离散形式的风险函数作为经验似然的标准，并以离散的风险函数作为我们的经验似然估计目标。

令 $w_i=\Delta\Lambda(T_i)\,,\, i=1,2,\ldots,n$ , 其中 $w_n = \delta_n$ ，因为累积风险函数的最后一跳必须为 1，以确保所有时刻的累积风险跳跃总和为 1。于是，经验似然函数可表示为：
$$EL(\Lambda)=\prod_{i=1}^n\left\{w_i\right\}^{\delta_i}\exp\left\{-\sum_{j=1}^nw_j\mathbb{1}[T_j\leq T_i]\right\}$$
接下来，设置经验似然估计的约束条件，即 $\Lambda(T)$ 需要满足以下积分约束：
$$\int_0^\infty g_1(t)\mathrm{d}\Lambda(t)=\theta_1$$
$$\int_0^\infty g_2(t)\mathrm{d}\Lambda(t)=\theta_2$$
$$\cdots\cdots$$
$$\int_0^\infty g_p(t)\mathrm{d}\Lambda(t)=\theta_p$$
其中，$g_i(t)\,,\,i=1,2\ldots,p$ 是满足一定矩条件的函数，$\theta_i$ 是给定的约束条件。在离散估计下，这些约束可表示为：
$$\sum_{i=1}^n g_1(T_i)w_i=\theta_1$$
$$\sum_{i=1}^n g_2(T_i)w_i=\theta_2$$
$$\cdots\cdots$$
$$\sum_{i=1}^n g_p(T_i)w_i=\theta_p$$
如果在没有约束条件的情况下进行估计，最大化经验似然得到的风险函数为$$w_i = \Delta \hat{\Lambda}_{NA}(T_i) = \frac{\delta_i}{R_i}$$其中 $R_i = \sum_k \mathbb{1}[T_k \geq T_i]$ ，该结果被称为 Nelson-Aalen 估计器。

而在有约束条件下进行估计时，我们可以通过 Nelson-Aalen 估计器得到参数 $\theta$ 的 NPMLE 估计，即： 
$$\hat{\theta}_{NPMLE}=\int g(t)\mathrm{d}\hat{\Lambda}_{NA}(t)$$
### 1.3 最大化风险经验似然
前面我们提到，离散风险函数的最后一个跳跃必须为 1，因此，在最大化经验似然估计时，最后一个跳跃应满足 $w_n=\Delta\hat{\Lambda}_{NA}(T_n)=\delta_n$ 。对于其他跳跃的估计，我们可以基于离散估计下的约束条件，通过以下定理来获得其他跳跃的估计值。

>[!info|wbk wtb tm]+ 
> *定理 1.1*&emsp;如果所有约束都具有唯一的风险函数解，那么在约束条件下，最大化对数经验似然 $\log EL(\Lambda)$ 的解在所有跳跃满足以下式子时获得，即：
> $$\begin{align}
> 	w_i &= \frac{\delta_i}{R_i+n\vec{\lambda}^\top \vec{G}(T_i)\delta_i} \\
> 	&= \frac{\delta_i}{R_i}\times \frac{1}{1+\vec{\lambda}^\top(\delta_i\vec{G}(T_i)\,/\,(R_i/n))} \\
> 	&= \Delta\hat{\Lambda}_{NA}(T_i)\frac{1}{1+\vec{\lambda}^\top \vec{Z}_i}\,,\quad i=1,2,\ldots,n
> \end{align}$$
> 其中
> $$\vec{G}(T_i)=\{g_1(T_i),\ldots,g_p(T_i)\}^\top\,,\quad\vec{Z}_i=\frac{\delta_i\vec{G}(T_i)}{R_i/n}=\left\{\frac{\delta_i g_1(T_i)}{R_i/n},\ldots,\frac{\delta_i g_p(T_i)}{R_i/n}\right\}$$
> 而 $\vec{\lambda}=\{\lambda_1,\ldots,\lambda_p\}$ 为下列等式的解
> $$\sum_{i=1}^{n-1}\frac{1}{n}\frac{Z_{ki}}{1+\lambda^\top\vec{Z}_i}+g_k(T_n)\delta_n=\theta_k\,,\quad k=1,2,\ldots,p$$
> 
> *Ps. 上述过程实际上是通过拉格朗日法求解在约束条件下使目标函数 $\log EL(\Lambda)$ 取最大值的解，$\lambda$ 为拉格朗日乘子。*

>[!info|wbk wtb tm]+ 
> *例子 1.1*&emsp;在本例中，我们检验一个关于生存概率的假设。具体来说，目的是检验在 3 年时的累积风险是否等于给定值 $\theta = -\log(0.4)$。为此，我们定义约束 $g(t) = \mathbb{1}[t \leq 1095.75]$，表示生存时间小于或等于 3 年，并假设原假设为 $H_0 : \int g(t) \mathrm{d}\Lambda_0(t) = \theta = -\log(0.4)$。根据生存函数的关系 $S(t) = \exp(-\Lambda_0(t))$，这个检验等价于检验 3 年时的生存概率 $S(3) = \exp(\log(0.4)) = 0.6$ 是否成立，即检验 3 年时的生存概率是否为 60%。