# Functions provided by PT_CY5.sage



## Example 1.13

After loading the package
```python
sage> load(PT_CY5.sage)
```
and setting up the ring for our calculation
```python
sage> setup_ring('q0,q1,q2,q3')
```

```python
sage> E([1], [-2, 0, 0, 0], [q1/q0, q2*q0^2, q3, 1, 1/(q0*q1*q2*q3)])
((q0^2*q1^2*q2^2*q3^2 - q3^2)/(-q3^2 + 1))*q^2
```

```python
sage> V([[1], [], []], [q1/q0, q2*q0^2, q3, 1, 1/(q0*q1*q2*q3)], 3)
1 + ((q0^4*q1^2*q2^2*q3^2 - q1^2)/(q0^2 - q1^2))*q^2 + ((q0^10*q1^4*q2^4*q3^4 - q0^6*q1^4*q2^2*q3^2 - q0^4*q1^6*q2^2*q3^2 + q1^6)/(q0^6 - q0^4*q1^2 - q0^2*q1^4 + q1^6))*q^4 + O(q^6)
```

## Example 1.14
