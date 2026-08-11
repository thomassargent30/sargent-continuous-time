# 4. The Concept of "Physical Realizability"

We have asserted that the white noise "$dN(t)/dt$" does not exist as an ordinary
stochastic process, it being so erratic that its variance is infinite. We have shown how
this process can be regarded as a limiting point of a sequence of ordinary stochastic
processes. It is useful to say a few words about the meaning of the concept of "physical
realizability" of a stochastic process. We will use the concept of physical realizability
to be interchangeable with the notion of an ordinary stochastic process.

We use the following definition:

**Definition 7.** A stochastic process is said to be *ordinary* or *physically realizable*
if its "realizations" or "sample paths" can be represented as ordinary function of time.

Loosely speaking, this means that it is in principle possible to "draw" each sample path
as a function of time. (It may, however, sometimes take a long time to do so, since for
example, the sample path of the Wiener process is not of bounded variation.)

We also use the following definition:

**Definition 8.** A *generalized stochastic process* is a stochastic process whose sample
paths cannot be represented as ordinary functions, but only as limit points of sequences
of ordinary functions. The sample paths of a "generalized stochastic process" can only be
represented as "generalized functions." The autocorrelation function of a generalized
stochastic process will itself be a generalized function.

Loosely speaking, sample paths of a generalized stochastic process cannot be "drawn," but
can only be represented in terms of ideal "pulses" of zero width but positive "mass." For
this reason, a generalized stochastic process is said not to be "physically realizable."

A generalized stochastic process is so erratic that it cannot be drawn. However, often
there exists a moving average of a generalized stochastic process that is physically realizable.
The shot-noise process $Y(t)$ of {doc}`03_poisson_counting_process`, formed by integrating
the generalized white noise $dN/dt$ against a square-integrable kernel, is the prototype: the
noise lives only under the integral sign, while $Y(t)$ itself is an ordinary process.
In our economic models to be constructed below, generalized stochastic processes will
appear only under integral signs. The economic models to be used always imply that the
observable economic variables are physically realizable. This is exactly the role white
noise plays in the moving-average (Wold) representation of {doc}`08_spectral_densities`, and
the admissibility of a physically unrealizable white-noise input recurs in the estimation
problem of {doc}`21_phillips_continuous_time_estimation`.
