rolling_shape
=============

![Rolling pentagon](https://github.com/mjwillatt/rolling_shape/blob/main/figures/rotated_output.gif?raw=true)

Description
-----------

This program produces images of a shape rolling without slipping on a curve.
The curve is determined from the shape and the rolling without slipping
constraint. To find the trajectory we need to solve a first order differential
equation, which is achieved numerically using separation of variables and
Simpson's rule.

The Python3 modules NumPy, SciPy, Matplotlib and argparse are required.

Example usage
-------------

The following command produces 1000 images (nframes=1000) of a regular pentagon
(n=5) completing 10 revolutions (ncycles=10), rolling from right to left
(reverse),

./rolling_shape.py --n 5 --nframes 1000 --ncycles 10 --reverse

Using the following command, these images are rotated 180 degrees (convert) and
combined into a gif (ffmpeg), 

./make_gif.sh

The resulting gif is displayed at the top of this README.

Background theory
-----------------

Suppose we have a shape described by a closed parametric curve $r = r(\theta)$,
where $\theta$ is measured anticlockwise from the negative $y$ axis. The shape
undergoes rigid body motion: clockwise rotation about the origin through an
angle $\phi$ followed by translation by $t$ along the $x$ axis. At time $t$,
the shape is described parametrically by $r = r(\theta + \phi)$, where the
moving origin is the point $(t, 0)$. Suppose the point on the shape directly
below the moving origin is in contact with a curve. This point has Cartesian
coordinates $(t, -r(\phi))$ and velocity $(1 - r(\phi)\dot{\phi}, 0)$. If the
motion corresponds to rolling without slipping, the velocity must be zero,

$$
\begin{equation}
r(\phi) \dot{\phi} = 1.
\end{equation}
$$

Solving this differential equation yields $\phi = \phi(t)$, which completely
determines the trajectory of the shape.

The following analysis applies to any polygon, but let us take the special case
of a square. $r(\theta)$ is determined piecewise as $r(\theta) = w\sec(\theta)$
for $-\pi/4 \le \theta \lt \pi/4$, where $w$ is half the width of the square.
Integrating the differential equation leads to

$$
\begin{equation}
\textrm{gd}^{-1}(\phi) = t / w,
\end{equation}
$$

where the inverse of the Gudermannian function is $\textrm{gd}^{-1}(\phi) =
\textrm{arcsinh}(\tan \phi)$. Therefore,

$$
\begin{align}
r(t) &= w \sec(\textrm{gd}(t/w)) \\
&= w \cosh(t/w).
\end{align}
$$

The curve on which the square rolls $(t, -w\cosh(t/w))$ is a catenary.
