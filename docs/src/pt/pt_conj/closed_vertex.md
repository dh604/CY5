# Closed vertex

The so-called closed vertex $X$ is discussed in [main_paper; Section 2.4](@cite). A closed formula for the PT series of $X$ is presented in Conjecture 2.11. Evidence 2.12 for this conjecture can be reproduced by running the script `closed_vertex.sage`. For convenience, we also included a Jupyter notebook `closed_vertex.ipynb` containing multiple comments on our code.

The script will generate the following output:

```
Expansion of the plethystic exponential relating membrane indices and the PT series took 475.55h.
Now we check the formulae against PT series in multi-degree (d1, d2, d3) with d1 + d2 + d3 <= 5 for 6 orders beyond the leading order in q.

[d1, d2, d3] = [1, 0, 0]  -->  True  (0.00h)
[d1, d2, d3] = [2, 0, 0]  -->  True  (0.01h)
[d1, d2, d3] = [1, 1, 0]  -->  True  (0.00h)
[d1, d2, d3] = [3, 0, 0]  -->  True  (0.34h)
[d1, d2, d3] = [2, 1, 0]  -->  True  (0.31h)
[d1, d2, d3] = [1, 1, 1]  -->  True  (1.73h)
[d1, d2, d3] = [4, 0, 0]  -->  True  (8.58h)
[d1, d2, d3] = [3, 1, 0]  -->  True  (20.51h)
[d1, d2, d3] = [2, 2, 0]  -->  True  (14.31h)
[d1, d2, d3] = [2, 1, 1]  -->  True  (33.83h)
[d1, d2, d3] = [5, 0, 0]  -->  True  (122.60h)
[d1, d2, d3] = [4, 1, 0]  -->  True  (375.96h)
[d1, d2, d3] = [3, 2, 0]  -->  True  (274.94h)
[d1, d2, d3] = [3, 1, 1]  -->  True  (333.54h)
[d1, d2, d3] = [2, 2, 1]  -->  True  (353.85h)
```
Each line corresponds to a check of the formula for $\mathsf{PT}_{(d_1,d_2,d_3)}(Z,\mathbb{C}^{\times}_q,\mathsf{T}')$ where $Z = X \times \mathbb{C}^2$. If the conjectured formula holds true we print `True`. In brackets we provide information on how long each check took.
