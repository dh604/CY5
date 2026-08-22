# Details of the Pipeline

```@meta
DocTestSetup = quote
    using Oscar, GKMtools, GW_CY5
end
```

On this page, we present a detailed example illustrating the [Pipeline](pipeline.md) and how it works internally.
The main function [`get_Omega_beta`](@ref) and its two most important helpers [`get_GW_beta`](@ref) and [`get_Omega_prediction`](@ref) have already been introduced on the page [Pipeline](pipeline.md).
We will therefore focus on the internal features that make these functions work.

Throughout, we assume that the user has already loaded all required packages using the following line.

```julia-repl
julia> using Oscar, GKMtools, GW_CY5
```

## Loading the space

Throughout, we will be dealing with the Calabi-Yau 5-fold

```math
Z = \mathcal{A}_2\times\mathbb{C}^3
```
and its equivariantly Calabi-Yau $T$-action given by its description as a _strip geometry_, see [main_paper; Section 2.1](@cite).
Let $G$ be the GKM graph of $(Z,T)$ and let $\beta\in H_2(Z)$ be the curve class associated to the leftmost edge in [main_paper; picture (18)](@cite).
We load $G$ using the function [`Ar_times_C3`](@ref).

```jldoctest pipeline_big_example
julia> G = Ar_times_C3(2)
GKM graph with 3 nodes, valency 5 and axial function:
2,1 -> 1,1 => (1, -1, 0, 0, 0)
3,1 -> 2,1 => (1, -2, 0, 0, 0)
Standalone flags:
1,1.2 => (0, -1, 0, -1, -1)
1,1.3 => (1, 0, 0, 0, 0)
1,1.4 => (0, 0, 0, 1, 0)
1,1.5 => (0, 0, 0, 0, 1)
2,1.3 => (0, -1, 0, -1, -1)
2,1.4 => (0, 0, 0, 1, 0)
2,1.5 => (0, 0, 0, 0, 1)
3,1.2 => (0, -1, 0, -1, -1)
3,1.3 => (-1, 3, 0, 0, 0)
3,1.4 => (0, 0, 0, 1, 0)
3,1.5 => (0, 0, 0, 0, 1)
```

We also load $\beta$.

```jldoctest pipeline_big_example
julia> beta = curve_class(G, Edge(1, 2))
(1, 0)
```

## Setting conjectures for $\Omega_\beta$

Let us see what the predictions for $\Omega_\beta$ and $\Omega_{2\beta}$ are up to genus one.
As discussed in [Pipeline](@ref), this is given to us by the function [`get_Omega_prediction`](@ref) as follows.

```jldoctest pipeline_big_example
julia> R, (t1, t2, t3, t4, t5, u) = polynomial_ring(QQ, ["t1", "t2", "t3", "t4", "t5", "u"]);

julia> max_genus = 1;

julia> get_Omega_prediction(G, [t1, t2, t3, t4, t5], u, beta, max_genus)
(-1//12*t3^2*t4*u^2 - 1//12*t3^2*t5*u^2 - 1//12*t3*t4^2*u^2 - 1//4*t3*t4*t5*u^2 - 1//12*t3*t5^2*u^2 - t3 - 1//12*t4^2*t5*u^2 - 1//12*t4*t5^2*u^2 - t4 - t5)//(t3*t4*t5*u^2)

julia> get_Omega_prediction(G, [t1, t2, t3, t4, t5], u, 2*beta, max_genus)
0
```

### Where does this come from?

Internally, [`get_Omega_prediction`](@ref) does the following.

```jldoctest pipeline_big_example
julia> get_attribute(G, :example_type)
:gkm_5d_strip_from_3d_CY
```

This tells the function that our space $G$ is of the type _strip geometry_.
It then calls the *internal* function `GW_CY5.gkm_5d_strip_prediction` which produces the specific
expression of $\Omega_\beta$ that is conjectured to hold for strip geometries:

```jldoctest pipeline_big_example
julia> GW_CY5.gkm_5d_strip_prediction(G, [t1, t2, t3, t4, t5], u, beta, max_genus)
(-1//12*t3^2*t4*u^2 - 1//12*t3^2*t5*u^2 - 1//12*t3*t4^2*u^2 - 1//4*t3*t4*t5*u^2 - 1//12*t3*t5^2*u^2 - t3 - 1//12*t4^2*t5*u^2 - 1//12*t4*t5^2*u^2 - t4 - t5)//(t3*t4*t5*u^2)

julia> GW_CY5.gkm_5d_strip_prediction(G, [t1, t2, t3, t4, t5], u, 2*beta, max_genus)
0
```

Similarly, all other spaces listed in [Spaces](spaces.md) have the attribute `:example_type`.
For each example type, we have implemented a dedicated prediction function that is called by [`get_Omega_prediction`](@ref), based on this attribute.

### General vs. manual prediction functions

Most conjectured formulae for $\Omega_\beta$ have a sufficiently compact form that we could implement directly so that the degree is an arbitrary parameter (_"general prediction"_).
This applies to all spaces except those listed in [Further Examples](spaces.md#Further-examples).
For those, the only way to arrive at the conjectural expression of $\Omega_\beta$ was by degree-wise computations in PT theory, where the size of the output grows rapidly as the degree increases.
For these examples, we saved the resulting expressions in the folder `test/omega_data/` (_"manual prediction"_).
These expressions are then loaded into Julia by the internal predictor function of the example in question.

For example, the internal prediction function for the space $\text{Tot}_{\mathbb{P}^2}(\mathcal{O}(-1)\oplus\mathcal{O}(-1)\oplus\mathcal{O}(-1))$ given by
[`gkm_5d_P2_111`](@ref) is `GW_CY5.gkm_5d_P2_111_prediction`.
It reads the files `test/omega_data/Omega_P2_111_num.dat` and `test/omega_data/Omega_P2_111_den.dat` to produce the following output in degree one and genus up to 2.

```jldoctest pipeline_big_example
julia> Tot_P2_111 = gkm_5d_P2_111();

julia> get_attribute(Tot_P2_111, :example_type)
:P2_111

julia> b = curve_class(Tot_P2_111, Edge(1, 2));

julia> GW_CY5.gkm_5d_P2_111_prediction(Tot_P2_111, [t1, t2, t3, t4, t5], u, b, 2)
1//32*t1^2*u^2 - 1//32*t1*t2*u^2 - 1//32*t1*t3*u^2 + 1//32*t2^2*u^2 - 1//32*t2*t3*u^2 + 1//32*t3^2*u^2 - 1//8
```

Indeed, this matches the output of [`get_Omega_prediction`](@ref):

```jldoctest pipeline_big_example
julia> get_Omega_prediction(Tot_P2_111, [t1, t2, t3, t4, t5], u, b, 2)
1//32*t1^2*u^2 - 1//32*t1*t2*u^2 - 1//32*t1*t3*u^2 + 1//32*t2^2*u^2 - 1//32*t2*t3*u^2 + 1//32*t3^2*u^2 - 1//8
```

We now return to the running example $Z = \mathcal{A}_2\times\mathbb{C}^3$ of this page.

## Getting $\Omega_\beta$ from GW invariants

Next, we explain how [`get_Omega_beta`](@ref) computes the membrane indices $\Omega_\beta$ from the truncated Gromov-Witten generating series produced by [`get_GW_beta`](@ref).
In our running example $Z = \mathcal{A}_2\times\mathbb{C}^3$, this produces the following in degree up to 2 and maximum genus 1.

```jldoctest pipeline_big_example
julia> GW = get_GW_beta(G, [beta, 2*beta], max_genus; show_bar=false);
Calculating b=(1, 0), g=0
Calculating b=(1, 0), g=1
Calculating b=(2, 0), g=0
Calculating b=(2, 0), g=1

julia> GW[beta]
(1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)

julia> GW[2*beta]
(1//24*t2^2*t4*u^2 + 1//24*t2^2*t5*u^2 + 1//24*t2*t4^2*u^2 + 1//24*t2*t4*t5*u^2 + 1//24*t2*t5^2*u^2 - 1//8*t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)
```

Recall the defining formula for $\Omega_\beta$:

```math
\sum_{\beta\neq 0}\sum_{g\ge 0} Q^\beta u^{2g-2}  GW_{g,\beta}(X,T)
=
\sum_{\beta\neq 0}
\sum_{k>0} \frac{1}{k} Q^{k\beta} \Omega_\beta(\epsilon,ku)
```
Extracting the coefficient of $Q^\beta$ gives
```math
\sum_{g\ge 0} u^{2g-2}GW_{g,\beta}(X,T)
=
\sum_{k\mid \beta}
\frac{1}{k} \Omega_{\beta/k}(\epsilon,ku)
```
On the right, the sum is over all positive integers $k$ dividing $\beta$.
Therefore, $\Omega_\beta$ can be computed if we know
- the left hand side,
- and $\Omega_{\beta/k}$ for all $k> 1$ dividing $\beta$.

Thus, we can compute $\Omega_\beta$ recursively, which is implemented as follows:

### Determining all relevant curve classes

If we want to compute $\Omega_\beta$ for all $\beta$ in some set $S$ of curve classes, the above
recursion requires us to compute $\Omega_\beta$ and $\sum_{g=0}^{\text{gMax}} u^{2g-2}GW_{g,\beta}(X,T)$
for all $\beta$ that have some positive integer multiple in $S$.

The internal function `GW_CY5.downward_close_cc` takes a list $S$ of curve classes and returns all curve classes with some
positive integer multiple contained in $S$:

```jldoctest pipeline_big_example
julia> GW_CY5.downward_close_ccs([6*beta, 10*beta])
6-element Vector{AbstractAlgebra.Generic.FreeModuleElem{ZZRingElem}}:
 (6, 0)
 (3, 0)
 (2, 0)
 (1, 0)
 (10, 0)
 (5, 0)
```

Indeed, for $S=\{6\beta, 10\beta\}$, the set of curve classes with positive integer multiple in $S$ is
$\{\beta, 2\beta, 3\beta, 5\beta, 6\beta, 10\beta\}$.
The function [`get_Omega_beta`](@ref) calls [`get_GW_beta`](@ref), which in turn uses `GW_CY5.downward_close_cc`
to determine all relevant curve classes.

### Solving the recursion for $\Omega_\beta$

Once [`get_GW_beta`](@ref) returned the truncated Gromov-Witten series in all relevant curve classes,
[`get_Omega_beta`](@ref) proceeds by computing $\Omega_{\beta}$ using the above equation.
This is done using the internal function `GW_CY5.cc_mobius`.
First it sorts the curve classes so that for any integer $k>1$, the curve class $k\beta$ is considered after $\beta$.
Second, it computes $\Omega_\beta$ in the given order.
This ensures that when $\Omega_\beta$ is computed, all relevant $\Omega_{\beta/k}$ have already been computed before.

Let us see an example of this.
Recall that we have previously computed the Gromov-Witten series of $Z=\mathcal{A}_2\times\mathbb{C}^3$ in degrees 1 and 2
up to genus 1:

```jldoctest pipeline_big_example
julia> GW[beta]
(1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)

julia> GW[2*beta]
(1//24*t2^2*t4*u^2 + 1//24*t2^2*t5*u^2 + 1//24*t2*t4^2*u^2 + 1//24*t2*t4*t5*u^2 + 1//24*t2*t5^2*u^2 - 1//8*t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)
```

Let us see how the internal function `GW_CY5.cc_mobius` turns this into $\Omega_\beta$.

```jldoctest pipeline_big_example
julia> Omega_from_cc_mobius = GW_CY5.cc_mobius(GW);

julia> Omega_from_cc_mobius[beta]
(1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)

julia> Omega_from_cc_mobius[2*beta]
0
```
Indeed, this is the output of [`get_Omega_beta`](@ref):

```jldoctest pipeline_big_example
julia> Omega = get_Omega_beta(G, [beta, 2*beta], max_genus; show_bar=false);
Calculating b=(1, 0), g=0
Calculating b=(1, 0), g=1
Calculating b=(2, 0), g=0
Calculating b=(2, 0), g=1

julia> Omega[beta]
(1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)

julia> Omega[2*beta]
0

julia> Omega[beta] == Omega_from_cc_mobius[beta]
true

julia> Omega[2*beta] == Omega_from_cc_mobius[2*beta]
true
```

## Testing conjectures for $\Omega_\beta$

Now that we know how [`get_Omega_beta`](@ref) computes the true values of $\Omega_\beta$ (up to some genus) from
Gromov-Witten invariants, let us see how it tests the true values vs. the conjectured values when
its optional argument `check_predictions=true` is set.

### Applying the equivariantly CY substitution

Note that the true value `Omega[beta]` above does not match the predicted value of $\Omega_\beta$ from above yet, which was as follows:

```jldoctest pipeline_big_example
julia> p = get_Omega_prediction(G, [t1, t2, t3, t4, t5], u, beta, max_genus)
(-1//12*t3^2*t4*u^2 - 1//12*t3^2*t5*u^2 - 1//12*t3*t4^2*u^2 - 1//4*t3*t4*t5*u^2 - 1//12*t3*t5^2*u^2 - t3 - 1//12*t4^2*t5*u^2 - 1//12*t4*t5^2*u^2 - t4 - t5)//(t3*t4*t5*u^2)

julia> p == Omega[beta]
false
```

The reason is that the following attribute is set.
```jldoctest pipeline_big_example
julia> has_attribute(G, :equiCY_substitution)
true
```
This means that we still need to impose the Calabi-Yau condition on the equivariant parameters of the predicted value.
To know what the right substitution is, [`get_Omega_beta`](@ref) internally does the following.

```jldoctest pipeline_big_example
julia> get_attribute(G, :equiCY_substitution)
5-element Vector{QQMPolyRingElem}:
 t1
 t2
 -t2 - t4 - t5
 t4
 t5
```
This means that, to impose the equivariantly Calabi-Yau condition, we need to substitute
$\epsilon_3 = -\epsilon_2-\epsilon_4-\epsilon_5$.
Indeed, the difference between the true value `Omega[beta]` and the predicted value `p` before imposing the Calabi-Yau condition
can be seen to be divisible by $\epsilon_2+\epsilon_3+\epsilon_4+\epsilon_5$:

```jldoctest pipeline_big_example
julia> factor(numerator(p - Omega[beta]))
(-1//12) * (t2 + t3 + t4 + t5) * (t4 + t5) * (t2*t3*u^2 + t3*t4*u^2 + t3*t5*u^2 + t4*t5*u^2 + 12)
```
So let us do the substitution.

```jldoctest pipeline_big_example
julia> t = [t1, t2, t3, t4, t5];

julia> CY_subst = vcat([evaluate(x, t) for x in get_attribute(G, :equiCY_substitution)], [u])
6-element Vector{QQMPolyRingElem}:
 t1
 t2
 -t2 - t4 - t5
 t4
 t5
 u

julia> prediction_equi_CY = evaluate(p, CY_subst)
(1//12*t2^2*t4*u^2 + 1//12*t2^2*t5*u^2 + 1//12*t2*t4^2*u^2 + 1//12*t2*t4*t5*u^2 + 1//12*t2*t5^2*u^2 - t2)//(t2*t4*t5*u^2 + t4^2*t5*u^2 + t4*t5^2*u^2)
```
This is now the predicted value of $\Omega_\beta$ after taking into account that we have a Calabi-Yau relation between the
equivariant parameters.
As expected, it matches the true value of $\Omega_\beta$ computed from Gromov-Witten invariants:

```jldoctest pipeline_big_example
julia> prediction_equi_CY == Omega[beta]
true
```

## Conclusion

Let us see how the entire workflow described above can be executed in a single line of code to test
if our conjectured values of $\Omega_\beta$ and $\Omega_{2\beta}$ are correct up to genus 1
in the example of $Z=\mathcal{A}_2\times\mathbb{C}^3$.

```jldoctest pipeline_big_example
julia> get_Omega_beta(G, [beta, 2*beta], max_genus; check_predictions=true, show_bar=false);
Calculating b=(1, 0), g=0
Calculating b=(1, 0), g=1
Calculating b=(2, 0), g=0
Calculating b=(2, 0), g=1
Prediction holds for (2, 0)
Prediction holds for (1, 0)
All predictions hold.
```

This shows that the conjectured values of $\Omega_\beta$ and $\Omega_{2\beta}$ are correct up to genus 1.
