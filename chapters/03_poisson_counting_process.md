# 3. The Poisson Counting Process

A *Poisson counting process* $N(t)$ is a continuous time stochastic process that assumes
values on the nonnegative integers, and which can be defined in the following ways. Let
$p_n(t) = \text{Prob}\, \{N(t) = n\},\ t \geq 0,\ n = 0,\ 1,\ 2,\ \ldots$. We define
$N(t)$ by a collection of differential equations for $p_n(t)$:

```{math}
:label: eq-3-1
\begin{aligned}
\frac{d}{dt}\ p_0(t) &= -\lambda p_0(t), p_0(0) = 1 \\
\frac{d}{dt}\ p_1(t) &= -\lambda p_1(t) + \lambda p_0(t) \\
\frac{d}{dt}\ p_j(t) &= -\lambda p_j(t) + \lambda p_{j-1}(t),\ j \geq 1 \\
\end{aligned}
```

where $\lambda > 0$, and $p_j(0) = 0$ for $j \geq 1$. These equations can be arranged
conveniently in the form of the infinite matrix equation.

$$
\frac{d}{dt}\ \begin{pmatrix} p_0(t) \\ p_1(t) \\ p_2(t) \\ \vdots \end{pmatrix}\ =\
\begin{pmatrix} -\lambda & 0 & 0 & \ldots \\ \lambda & -\lambda & 0 & \ldots \\ 0 & \lambda & -\lambda & \ldots \end{pmatrix}\ \begin{pmatrix} p_0(t) \\ p_1(t) \\ p_2(t) \\ \vdots \end{pmatrix}
$$

The positive parameter $\lambda$ is called the *rate* of the process. The solution of this
recursive system of differential equations is readily found to be

$$
\begin{aligned}
p_0(t) &= e^{-\lambda t} \\
p_1(t) &= \lambda t e^{-\lambda t} \\
\ \vdots \\
p_j(t) &= \frac{(\lambda t)^j}{j!}\ e^{-\lambda t},\, j \geq 1 \\
\end{aligned}
$$

Evidently, we could define the Poisson counting process by simply starting with the
definitions of probabilities

$$
p_n(t) = \text{Prob } \{N(t) = n\} = \frac{(\lambda t)^n}{n!}\ e^{-\lambda t},\ n = 0,\ 1,\ \ldots
$$

The mean of the Poisson counting process is readily calculated as

$$
\begin{aligned}
EN(t) &= \sum_{n=0}^{\infty}\ n p_n(t) \\
&= \sum_{n=0}^{\infty}\ n\ \frac{(\lambda t)^n}{n!}\ e^{-\lambda t} = \lambda t\ \sum_{n=0}^{\infty}\ \frac{(\lambda t)^{n-1}}{(n-1)!}\ e^{-\lambda t} \\
\end{aligned}
$$

or

$$
EN(t) = \lambda t.
$$

This formula motivates the interpretation of $\lambda$ as the *rate* or *arrival rate* of
the process. Similar calculations show that

$$
EN(t)^2 = \lambda^2\, t^2 + \lambda t,
$$

which implies that

$$
E\, \big(N(t) - EN(t)\big)^2 = \lambda t.
$$

An example of a realization of a Poisson counting process is given in figure (—). The
random points in time $t_1,\ t_2,\ \ldots$ at which $N(t)$ jumps are called *arrival
times*.

```{note}
*(Figure omitted — hand-drawn in the original manuscript.)*
```

Notice that the event that the first arrival time $t_1 \leq T_1$ is the event that
$N(T_1) \geq 1$, which is the complement of the event that $N(T_1) = 0$. Therefore, we
have that

$$
\text{Prob } \{t_1 \leq T_1\} = 1 - e^{-\lambda T_1},
$$

which states that the first arrival time is exponentially distributed.

From the above, it follows that a Poisson process $N(t)$ is characterized by the
properties

(i) $N(0) = 0$

(ii) For $t_1 > t_2$,

$$
\begin{aligned}
\text{Prob } &\big\{N(t_1) - N(t_2) = k \big\} =\ \frac{(\lambda (t_1 - t_2))^k}{k!}\ e^{-\lambda (t_1 - t_2)}, \\
k &= 0,\ 1,\ \ldots \\
\end{aligned}
$$

(iii) If $t_1 > t_2 > t_3 > t_4$, then $(N(t_1) - N(t_2))$ is statistically independent of
$N(t_3) - N(t_4)$, and $E(N(t_1) - N(t_2))\, (N(t_3) - N(t_4)) = \lambda^2 (t_1 - t_2)\, (t_3 - t_4)$,
so that the expectation of the products of the increments $(N(t_1) - N(t_2))$ and
$(N(t_3) - N(t_4))$ is the product of their expectations.

Property (iii) is the "independent increments" property of a Poisson process.

For $t_1 > t_2$, we have

$$
\begin{aligned}
E\, \big\{N(t_1) - N(t_2)\big\} &= \lambda (t_1 - t_2) \\
E\, \big\{N(t_1) - N(t_2)\}^2 &= \lambda^2 (t_1 - t_2)^2 + \lambda (t_1 - t_2). \\
\end{aligned}
$$

Assuming again that $t_1 > t_2 > t_3 > t_4$, and using
$N(t_1) - N(t_3) = (N(t_1) - N(t_2)) + (N(t_2) - N(t_3))$, and
$N(t_2) - N(t_4) = (N(t_2) - N(t_3)) + (N(t_3) - N(t_4))$, we obtain

$$
E\, \big\{ (N(t_1) - N(t_3))\, (N(t_2) - N(t_4)) = \lambda^2 (t_1 - t_3)\ (t_2 - t_4) + \lambda (t_2 - t_3).
$$

It follows that

$$
\text{cov}\ (N(t_1) - N(t_3))\ (N(t_2) - N(t_4)) = \lambda (t_2 - t_3),
$$

so that the covariance between two increments is proportional to their overlap in time
$(t_2 - t_3)$.

Using the preceding formulas, we can calculate the autocorrelation function
$R(t_1,\, t_2)$ and the autocovariance function $C(t_1,\, t_2)$ for the Poisson process.
We find

$$
R(t_1,\, t_2) = \begin{cases} \lambda t_2 &+ \lambda^2 t_1\, t_2,\ t_1 \geq t_2 \\ \lambda t_1 &+ \lambda^2 t_1\, t_2,\ t_1 \leq t_2 \end{cases}
$$

or

$$
R(t_1,\, t_2) = \lambda\, \min\, (t_1,\, t_2) + \lambda^2 t_1\, t_2.
$$

Note that $\min\, (t_1,\, t_2)$ is the length of overlap of $N(t_2)$ and $N(t_1)$. It
follows from (—) that

$$
C\, (t_1,\, t_2) = \lambda\, \min\, (t_1,\, t_2).
$$

From (—) it follows that $N(t)$ is mean square continuous. (Why?) From {eq}`eq-3-1`, we
also have that

$$
\frac{\partial R(t_1,\, t_2)}{\partial t_1}\ =\ \begin{cases} 0 &+ \lambda^2 t_2,\ t_1 \geq t_2 \\ \lambda &+ \lambda^2 t_2,\ t_2 \geq t_1 \end{cases}
$$

or

$$
\frac{\partial R(t_1,\, t_2)}{\partial t_1}\ =\ \lambda u(t_2 - t_1) + \lambda^2 t_2
$$

where $u(t)$ is the Heaviside unit step function defined by $u(t) = 1$ for $t \geq 0,\
u(t) = 0$ for $t < 0$. It follows from (—) that

$$
\frac{\partial^2 R(t_1,\, t_2)}{\partial t_1 \partial t_2}\ =\ \lambda \delta(t_2 - t_1) + \lambda^2
$$

where $\delta(\ \ )$ is the Dirac delta generalized function defined by

$$
\int_{-\infty}^{\infty} g(\tau)\, \delta(\tau)\, d\tau = g(0)
$$

for all test functions $g(t)$ that are continuous at $t = 0$, and that go to zero
sufficiently quickly as $|t| \to \infty$.

Equation (—) states that $\partial^2 R(t_1,\, t_2)/\partial t_1 \partial t_2$ does not
exist as an ordinary function, so that $R(t_1,\, t_2)$ is not twice differentiable. It
follows that $N(t)$ is not mean square differentiable, so that $dN(t)/dt$ does not exist
as an ordinary stochastic process.

Although $dN(t)/dt$ does not exist in the mean square sense, it does exist as a
"generalized stochastic process." The sample paths of the Poisson process can be
represented as

$$
N(t) = \sum_{i=1}^{\infty}\, u(t - t_i)
$$

where $\{ t_1,\, t_2,\, \ldots\}$ are the Poisson arrival times and $u(t)$ is again the
Heaviside unit step function. Taking the generalized derivative of the sample path $N(t)$
and using the definition $\frac{d}{dt}\, u(t - t_i) = \delta(t - t_i)$, we obtain

$$
z(t) \equiv\ \frac{dN(t)}{dt} = \sum_{i=1}^{\infty}\ \delta(t - t_i),
$$

so that $dN(t)/dt$ is a sum of $\delta$-impulses at the random arrival times $t_i$. From
our preceding results, it follows that

$$
Ez(t) = \lambda
$$

$$
Ez(t_1) z(t_2) = \lambda^2 + \lambda \delta(t_1 - t_2)
$$

The zero mean generalized stochastic process $z(t) - \lambda$, which has autocorrelation
function $\lambda \delta(t_1 - t_2)$ is an example of a *white noise*. Any stochastic
process with an autocorrelation function proportional to $\delta(t_1 - t_2)$ is known as a
white noise.

The process $z(t) = dN(t)/dt$ does not exist as an ordinary stochastic process, but can be
regarded as a particular kind of limit point of a process that does exist as an ordinary
stochastic process. In particular, for $\epsilon > 0$, define

$$
y(t)\ =\ \frac{N(t + \epsilon) - N(t)}{\epsilon}
$$

where $N(t)$ is a Poisson counter with rate $\lambda$. It follows that

$$
y(t) = k/\epsilon
$$

where $k = N(t + \epsilon) - N(t) =$ number of arrivals of the $N(t)$ process in the
interval $(t,\, t + \epsilon)$. Therefore,

$$
\text{Prob } \big\{y(t) = k/\epsilon\big\}\ =\ \frac{e^{-\lambda \epsilon} (\lambda \epsilon)^k}{k!}
$$

Using this Poisson probability distribution, and the above results on moments of a Poisson
process, we can deduce that

$$
Ey(t) = \lambda
$$

$$
R(t_1,\, t_2) = \begin{cases} \lambda^2 & |t_1 - t_2| > \epsilon \\ \lambda^2 &+ \frac{\lambda}{\epsilon}\ -\ \frac{\lambda |t_1 - t_2|}{\epsilon^2}\, ,\ \text{ for }\ |t_1 - t_2| \leq \epsilon \end{cases}
$$

The function $\lambda/\epsilon - \lambda |t_1 - t_2|/\epsilon^2$ for
$|t_1 - t_2| \leq \epsilon$ is plotted in figure (—), and inscribes a triangle of area
$\lambda$ above the horizontal axis. It is known that the limit as $\epsilon \to 0$, of the
ordinary function $\max\ (0,\, \frac{\lambda}{\epsilon}\ -\ \frac{\lambda |t_1 - t_2|}{\epsilon^2})$

```{note}
*(Figure omitted — hand-drawn in the original manuscript.)*
```

defines the Dirac delta generalized function with mass $\lambda$,
$\lambda \delta(t_1 - t_2)$.

It is a "spike" of "mass" $\lambda$ at $t_1 - t_2 = 0$, and is equal to zero for
$t_1 - t_2 \neq 0$.

As an example of one use of the generalized stochastic process $z(t) = dN(t)/dt$, let
$L(\tau),\ \tau \in [0,\, \infty)$ be a continuous and square integrable function, i.e.,
$\int_0^{\infty} L(\tau)^2\, d\tau < +\infty$. Then consider a stochastic process defined by
the distributed lag

$$
Y(t) = \int_0^{\infty} L(\tau) z(t - \tau)\, d\tau.
$$

Substituting $z(t) = \sum_{i=1}^{\infty}\, \delta(t - t_i)$ and using the definition of the
delta generalized function, we find

$$
\begin{aligned}
Y(t) &= \int_0^{\infty} L(\tau) \sum_{i=1}^{\infty} \delta(t - t_i - \tau)\, d\tau \\
&= \sum_{i=1}^{\infty} \int_0^{\infty} L(\tau) \delta(t - t_i - \tau)\, d\tau \\
Y(t) &= \sum_{i=1}^{\infty} L(t - t_i),\ t \geq 0 \\
\end{aligned}
$$

The process $Y(t)$ defined by (—) is an ordinary stochastic process, consisting of a sum
of the function $L$ shifted by the random arrival times $t_i$. Such a process $Y(t)$ is
called *shot noise*. In Section (—), we shall show how to calculate its first and second
order moments.
