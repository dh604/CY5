# Local $\mathbb{P}^3$

The purpose of the file `local_P3.sage` is to test our conjectural formulas for membrane indices of the fivefolds
$$ \mathrm{Tot}_{\mathbb{P}^3}\, \mathcal{O}(-a) \oplus \mathcal{O}(a-4), \qquad 0 \leq a \leq 4 $$
against PT calculations. The formulae for membrane indices are partly stated in [main_paper; Section 4.2](@cite). All remaining formulae can be found in `CY5\tests\omega_data\`.

The script `local_P2.sage` reads the files with our conjectural formulae and compares them with the PT series computed via the K-theoretic vertex. The script produces the following output, which we record in our paper as [main_paper; Evidence 4.9](@cite):

```
--------------------
Tot_P3 O(0) + O(-4)
--------------------
Finished reading membrane indices of Tot_P3 O(0) + O(-4). This took 0.01h.
Now we compute the plethystic exponential up to degree 3.

Expansion of conjectural formula took 0.04h.
Now we check the formula up to 7 orders beyond leading order in q.

Tot_P3 O(0) + O(-4), d = 1  -->  True  (0.00h)
Tot_P3 O(0) + O(-4), d = 2  -->  True  (0.57h)
Tot_P3 O(0) + O(-4), d = 3  -->  True  (136.98h)

--------------------
Tot_P3 O(-1) + O(-3)
--------------------
Finished reading membrane indices of Tot_P3 O(-1) + O(-3). This took 0.06h.
Now we compute the plethystic exponential up to degree 2.

Expansion of conjectural formula took 0.97h.
Now we check the formula up to 7 orders beyond leading order in q.

Tot_P3 O(-1) + O(-3), d = 1  -->  True  (0.00h)
Tot_P3 O(-1) + O(-3), d = 2  -->  True  (0.69h)

--------------------
Tot_P3 O(-2) + O(-2)
--------------------
Finished reading membrane indices of Tot_P3 O(-2) + O(-2). This took 0.00h.
Now we compute the plethystic exponential up to degree 3.

Expansion of conjectural formula took 0.00h.
Now we check the formula up to 7 orders beyond leading order in q.

Tot_P3 O(-2) + O(-2), d = 1  -->  True  (0.00h)
Tot_P3 O(-2) + O(-2), d = 2  -->  True  (0.31h)
Tot_P3 O(-2) + O(-2), d = 3  -->  True  (48.94h)

--------------------
Tot_P3 O(-3) + O(-1)
--------------------
Finished reading membrane indices of Tot_P3 O(-3) + O(-1). This took 0.23h.
Now we compute the plethystic exponential up to degree 2.

Expansion of conjectural formula took 0.97h.
Now we check the formula up to 7 orders beyond leading order in q.

Tot_P3 O(-3) + O(-1), d = 1  -->  True  (0.00h)
Tot_P3 O(-3) + O(-1), d = 2  -->  True  (0.68h)

--------------------
Tot_P3 O(-4) + O(0)
--------------------
Finished reading membrane indices of Tot_P3 O(-4) + O(0). This took 0.01h.
Now we compute the plethystic exponential up to degree 3.

Expansion of conjectural formula took 0.06h.
Now we check the formula up to 7 orders beyond leading order in q.

Tot_P3 O(-4) + O(0), d = 1  -->  True  (0.00h)
Tot_P3 O(-4) + O(0), d = 2  -->  True  (0.44h)
Tot_P3 O(-4) + O(0), d = 3  -->  True  (50.58h)
```
