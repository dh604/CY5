# Examples from the paper

On this page, we mirror the SageMath code presented in the paper.

## Example 1.13

In [main_paper; Example 1.13](@cite) we showcase an example calculation of the degree one PT invariants of the threefold $X_1 = \text{Tot}_{\mathbb{P}^1}(\mathcal{O}(-2))\times \mathbb{C}$ with a Nekrsov-Okounkov insertion coming from an ambient fivefold
$Z = \text{Tot}_{\mathbb{P}^1}(\mathcal{O}(-2))\times \mathbb{C}^3$ with respect to a certain torus action.

After setting up the ring for our calculation
```python
sage> setup_ring('q0,q1,q2,q3')
```
the PT series in curve class $\beta=[\mathbb{P}^1]$ is computed via
```python
sage> (
sage>   V([[1], [], []], [q1/q0, q2*q0^2, q3, 1, 1/(q0*q1*q2*q3)], 3)
sage>   * E([1], [-2, 0, 0, 0], [q1/q0, q2*q0^2, q3, 1, 1/(q0*q1*q2*q3)])
sage>   * V([[1], [], []], [q0/q1, q2*q1^2, q3, 1, 1/(q0*q1*q2*q3)], 3)
sage> )
((q0^2*q1^2*q2^2*q3^2 - q3^2)/(-q3^2 + 1))*q^2 + ((q0^4*q1^4*q2^4*q3^4 - q0^2*q1^2*q2^2*q3^4 + q0^2*q1^2*q2^2*q3^2 - q3^2)/(-q3^2 + 1))*q^4 + ((-q0^6*q1^6*q2^6*q3^6 + q0^4*q1^4*q2^4*q3^6 - q0^4*q1^4*q2^4*q3^4 + q0^2*q1^2*q2^2*q3^4 - q0^2*q1^2*q2^2*q3^2 + q3^2)/(q3^2 - 1))*q^6 + O(q^8)
```

## Example 1.14

Similarly, in [main_paper; Example 1.14](@cite) we compute the degree one PT series of the threefold $X_2 = \mathbb{P}^1 \times \mathbb{C}^2$ with a Nekrsov-Okounkov insertion coming from the ambient fivefold
$Z = \text{Tot}_{\mathbb{P}^1}(\mathcal{O}(-2))\times \mathbb{C}^3$ with respect to a certain torus action.

This PT series can be calculated calling
```python
sage> setup_ring('q0,q1,q4,q5')
sage> (
sage>   V([[1], [], []], [q1/q0, q4, q5, q0^2, 1/(q0*q1*q4*q5)], 3)
sage>   * E([1], [0, 0, -2, 0], [q1/q0, q4, q5, q0^2, 1/(q0*q1*q4*q5)])
sage>   * V([[1], [], []], [q0/q1, q4, q5, q1^2, 1/(q0*q1*q4*q5)], 3)
sage> )
((-q4^2*q5^2)/(q4^2*q5^2 - q4^2 - q5^2 + 1)) + ((q0^2*q1^2*q4^4*q5^4 - q0^2*q1^2*q4^2*q5^2)/(-q4^2*q5^2 + q4^2 + q5^2 - 1))*q^2 + ((q0^4*q1^4*q4^6*q5^6 - q0^4*q1^4*q4^4*q5^4)/(-q4^2*q5^2 + q4^2 + q5^2 - 1))*q^4 + O(q^6)
```
