# 7. Stochastic Differential Equations Driven by a Wiener Process

We wish to consider a class of stochastic differential equations driven by a Wiener process,

$$
dx(t) = f\left(x(t),\, t\right) dt + g\left(x(t),\,t\right) dW(t).
$$

We shall derive rules for manipulating such equations by regarding (—) as the limit as $\lambda \to \infty$ of

$$
dx(t) = f\left(x(t),\, t\right) dt + \frac{1}{\sqrt\lambda}\ g\left(x(t),\, t\right)\ \left(dN_1(t) - dN_2(t)\right),
$$

where $N_1$ and $N_2$ are two independent Poisson counters with identical rates $\lambda/2$.

The first rule we seek is the counterpart of rule (—). Where $\Psi(x)$ is a function of $x$, $\Psi(x)$ obeys the differential equation

$$
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, , \ f\left(x(t),\, t\right)\right\rangle dt + 1/2\, \left\langle \frac{\partial^2\Psi}{\partial x^2}\, g,\, g(x,\, t)\right\rangle dt \\
&+ \left\langle \frac{\partial\Psi}{\partial x}\, ,\ g(x,\, t)\right\rangle dW(t).
\end{aligned}
$$

To generate this rule, we use (—) for $\lambda > 0$ and using rule (—) to obtain

```{math}
:label: eq-7-1
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle dt \\
&+ \left[ \Psi\left(x + \frac{1}{\sqrt\lambda}\ g(x,\, t)\right) - \Psi(x) \right] dN_1(t) \\
&+ \left[ \Psi\left(x - \frac{1}{\sqrt\lambda}\ g(x,\, t)\right) - \Psi(x) \right] dN_2(t).
\end{aligned}
```

Obtain for $\Psi(x \pm\, \frac{1}{\sqrt\lambda}\ g(x,\, t))$ its Taylor series expansion about $x$,

$$
\begin{aligned}
\Psi\left(x \pm \frac{1}{\sqrt\lambda}\ g\left(x,\, t\right)\right) &= \Psi(x) \pm\ \left\langle \frac{\partial\Psi}{\partial x}\, ,\ \frac{1}{\sqrt\lambda}\ g(x,\, t)\right\rangle \\
&+ \frac{1}{2}\ \left(\frac{1}{\sqrt\lambda}\right)^2\ \left\langle \frac{\partial^2\Psi}{\partial x^2}\ g,\ g(x,\, t)\right\rangle + 0\, \left(\frac{1}{\lambda^{3/2}}\right)
\end{aligned}
$$

where $0\, (\frac{1}{\lambda^{3/2}}) \to 0$ as $\lambda \to \infty$. Substituting the Taylor expansions in {eq}`eq-7-1`, we obtain

```{math}
:label: eq-7-2
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle dt + \left\langle \frac{\partial\Psi}{\partial x}\, ,\ g(x,\, t)\right\rangle\ \frac{1}{\sqrt\lambda}\ (dN_1 - dN_2) \\
&+ \frac{1}{2}\ \left\langle \frac{\partial^2\Psi}{\partial x^2}\, g,\, g(x,\, t)\right\rangle\ \frac{1}{\lambda}\ (dN_1 + dN_2) + 0\, \left(\frac{1}{\lambda^{3/2}}\right)
\end{aligned}
```

Now consider the process $z(t)$ governed by

$$
dz(t) = \frac{1}{\lambda}\ (dN_1 + dN_2)
$$

Applying our rules, we find that

$$
\frac{d}{dt}\ Ez(t) = 1
$$

so that

$$
Ez(t) = t.
$$

We also find that

$$
\frac{d}{dt}\ Ez(t)^2 = 2Ez + 2/\lambda,
$$

so that as $\lambda \to \infty$

$$
\frac{d}{dt}\ Ez(t)^2 = 2Ez(t) = 2t.
$$

Therefore

$$
Ez(t)^2 = t^2.
$$

It follows that

$$
E\left(z(t) - Ez(t)\right)^2 = 0,
$$

which implies that as $\lambda \to \infty$,

$$
z(t) = t
$$

with probability 1.

Thus returning to {eq}`eq-7-2` and taking limits as $\lambda \to \infty$, we obtain

$$
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle dt + 1/2\, \left\langle \frac{\partial^2\Psi}{\partial x^2}\, g,\, g(x,\, t)\right\rangle dt \\
&+ \left\langle \frac{\partial\Psi}{\partial x}\, ,\ g(x,\, t)\right\rangle dW(t).
\end{aligned}
$$

This is known as *Ito's rule* for the stochastic differential equation

$$
dx = f(x,\, t)dt + g(x,\, t)dW(t).
$$

Our second rule is

$$
\frac{d}{dt}\ E\Psi(x) = E\ \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle + \frac{1}{2}\, E\ \left\langle \frac{\partial^2\Psi}{\partial x^2}\, g,\, g\right\rangle.
$$

This rule can be derived by the same limiting process. For $\lambda > 0$, we have

$$
\begin{aligned}
d\Psi(x) &= \left\langle \frac{\partial\Psi}{\partial x}\, , \ f(x,\, t)\right\rangle dt + \left( \Psi(x + \frac{1}{\sqrt\lambda}\ g(x,\, t)\right) \\
&- \Psi(x)\bigr) dN_1 + \left(\Psi(x - \frac{1}{\sqrt\lambda}\ g\left(x,\, t)\right) - \Psi(x)\right)dN_2
\end{aligned}
$$

which implies that

$$
\begin{aligned}
\frac{d}{dt}\ E\Psi(x) &= E\, \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle + E\, \left[\Psi(x + \frac{1}{\sqrt\lambda}\ g(x,\, t)) \right. \\
&+ \left. \Psi\left( x - \frac{1}{\sqrt\lambda}\ g(x,\, t) \right) - 2\ \Psi(x) \right]\lambda.
\end{aligned}
$$

Now take Taylor series expansions of $\Psi\, (x \pm\, \frac{1}{\sqrt\lambda}\ g(x,\, t))$ about $x$ to get

$$
\begin{aligned}
\frac{d}{dt}\ E\Psi(x) &= E\, \left\langle \frac{\partial\Psi}{\partial x}\, ,\ f(x,\, t)\right\rangle + E\, \left\langle \frac{\partial\Psi}{\partial x}\, ,\ \frac{1}{\sqrt\lambda}\ \left(g(x,\, t) - g(x,\, t)\right)\right\rangle \\
&+ \frac{1}{2}\ \left\langle \frac{\partial^2\Psi}{\partial x^2}\, ,\ g^2\right\rangle\ \frac{1}{\lambda} \cdot \lambda + 0(\lambda^{-1/2}).
\end{aligned}
$$

Taking the limit as $\lambda \to \infty$ gives the desired result (—).

The next result that we desire is for $\tau > 0$,

$$
\frac{d}{d\tau}\ Ex(t) x(t + \tau) = Ex\, (t) f\left(x(t+\tau),\ t+\tau\right).
$$

To obtain this, we take limits as $\lambda \to \infty$ in the formula

$$
\frac{d}{d\tau}\ Ex(t) x(t+\tau) = Ex\, (t)f\left(x(t+\tau),\ t+\tau\right) + Ex\, (t)\ \left[\frac{1}{\sqrt\lambda}\ \left( \frac{\lambda}{2}\ - \ \frac{\lambda}{2}\right)\right].
$$

As an example of the use of these formulas, we take the linear stochastic differential equation

$$
dx(t) = -ax(t)dt + bdW(t),\ x(0) = 0,\ a,\ b > 0.
$$

Applying our formulas, we find that

$$
\begin{aligned}
dx(t)^2 &= \left(-2ax(t)^2 + b^2\right) dt + 2bx(t)\, dW(t). \\
\frac{d}{dt}\ Ex(t)^2 &= -2a Ex(t)^2 + b^2 \\
\frac{d}{dt}\ Ex(t) &= -a Ex(t) \\
\frac{d}{d\tau}\ Ex(t) x(t + \tau) &= -aEx(t) x(t + \tau),\ \tau > 0.
\end{aligned}
$$
