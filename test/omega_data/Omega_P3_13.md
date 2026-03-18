# Membrane index of $\mathrm{Tot}_{\mathbb{P}^3}\mathcal{O}(-1)\oplus \mathcal{O}(-3)$

The files `Omega_P3_13_num.dat` and `Omega_P3_13_den.dat` contain the numerator
and denominator of the membrane index $\Omega_{d[L]}$ of
```math
Z = \mathrm{Tot}_{\mathbb{P}^3}\mathcal{O}(-1)\oplus \mathcal{O}(-3)
```
in degree $1 \leq d\leq 3$ respectively. Here, the $d$th line of the files contains
the numerator resp. denominator of $\Omega_{d[L]}$ which is stored as a polynomial
in the characters `q0,q1,q2,q3,q4` of a certain torus $T=(\mathbb{C}^{\times})^5$
acting on $Z$ and elementary symmetric polynomials `e1,e2,e3,e4` in the variables
`q0,q1,q2,q3`. The action of the torus $T$ is described as follows: Let us present
$Z$ as the quotient
```math
Z = \mathbb{C}^6 \setminus V(x_0 x_1 x_2 x_3) \big/ \mathbb{C}^{\times}
```
where $\mathbb{C}^{\times}$ acts on affine six-space via
```math
t\cdot(x_0,x_1,x_2,x_3,y_0,y_1) = (t x_0,t x_1,t x_2,t x_3, t^{-1} y_0, t^{-3} y_1).
```
Now let $T=(\mathbb{C}^{\times})^5$ act on the factors of $\mathbb{C}^6$ with characters
$q_0^{-2},\ldots,q_4^{-2}$ and $(q_0\cdots q_4)^{2}$ respectively. This descends to a
$T$-action on $Z$. With this for instance the tangent weights at the fixed point
$[1:0:0:0]\in \mathbb{P}^3 \subseteq Z$ read
```math
q_0^{-2}q_1^2, \quad q_0^{-2}q_2^2, \quad q_0^{-2} q_3^2, \quad q_0^2 q_4^2, \quad
q_0^{4} (q_1\cdots q_4)^{-2} .
```
The conjectural expressions for $\Omega_{d[L]}$ stored in
the files is expanded in exactly these $T$-characters $q_0,\ldots,q_4$.
