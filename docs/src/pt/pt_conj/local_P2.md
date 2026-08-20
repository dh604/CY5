(# Local $\mathbb{P}^2$

The purpose of the file `local_P2.sage` is to test different claims made in [main_paper](@cite) about PT series of threefolds inside local projective planes.

More precisely, for the fivefold
$$ \mathrm{Tot}_{\mathbb{P}^2}\, \mathcal{O}(-1) \oplus \mathcal{O}(-1) \oplus \mathcal{O}(-1) $$
we show that the PT series of the threefold $\mathrm{Tot}_{\mathbb{P}^2}\, \mathcal{O}(-1)$ inside this ambient space matches with some explicit formulae for membrane indices $\Omega_{dH}$ in degree $d \leq 4$ as discussed in [main_paper; Section 4.1](@cite).

For the fivefold
$$ Z=\mathrm{Tot}_{\mathbb{P}^2}\, \mathcal{O}(-2) \oplus \mathcal{O}(-1) \oplus \mathcal{O}  $$
we test that for all six distinct choices of fibrewise actions of $\mathbb{C}^{\times}_q$ on $Z$ with three-dimensional fixed locus $X$, the associated PT series vanishes. This is a special instance of [main_paper; Conjecture 5.6](@cite).

By running the script `local_P2.sage` one can hence reproduce Evidence 4.4 and parts of Evidence 5.7 in our paper. The other part of Evidence 5.7 is checked [here](local_P1xP1.md). The script produces the following output:

```
----------------------------
Tot_P2 O(-1) + O(-1) + O(-1)
----------------------------
Finished reading Omega. This took 5.29h.
Now we compute the plethystic exponential
Expansion of conjectural formula took 25.64h.
Now we check the formula up to 7 orders beyond leading order in q and degree d <= 4.

Tot_P2 O(-1) + O(-1) + O(-1), d = 1  -->  True  (0.00h)
Tot_P2 O(-1) + O(-1) + O(-1), d = 2  -->  True  (0.09h)
Tot_P2 O(-1) + O(-1) + O(-1), d = 3  -->  True  (16.93h)
Tot_P2 O(-1) + O(-1) + O(-1), d = 4  -->  True  (932.29h)

---------------------------------------------
Tot_P2 O(-2)  in  Tot_P2 O(-2) + O(-1) + O(0)
---------------------------------------------
We test the vanishing of the PT series up to 6 orders beyond leading order in q and degree d <= 4:
Tot_P2 O(-2) + O(-1) + O(0), d = 1  -->  True  (0.29h)
Tot_P2 O(-2) + O(-1) + O(0), d = 2  -->  True  (8.08h)
Tot_P2 O(-2) + O(-1) + O(0), d = 3  -->  True  (439.06h)
Tot_P2 O(-2) + O(-1) + O(0), d = 4  -->  True  (11595.40h)

---------------------------------------------
Tot_P2 O(-2)  in  Tot_P2 O(-2) + O(0) + O(-1)
---------------------------------------------
We test the vanishing of the PT series up to 6 orders beyond leading order in q and degree d <= 4:
Tot_P2 O(-2) + O(0) + O(-1), d = 1  -->  True  (0.07h)
Tot_P2 O(-2) + O(0) + O(-1), d = 2  -->  True  (7.76h)
Tot_P2 O(-2) + O(0) + O(-1), d = 3  -->  True  (417.85h)
Tot_P2 O(-2) + O(0) + O(-1), d = 4  -->  True  (11926.15h)

---------------------------------------------
Tot_P2 O(-1)  in  Tot_P2 O(-1) + O(-2) + O(0)
---------------------------------------------
We test the vanishing of the PT series up to 6 orders beyond leading order in q and degree d <= 4:
Tot_P2 O(-1) + O(-2) + O(0), d = 1  -->  True  (0.06h)
Tot_P2 O(-1) + O(-2) + O(0), d = 2  -->  True  (7.32h)
Tot_P2 O(-1) + O(-2) + O(0), d = 3  -->  True  (385.44h)
Tot_P2 O(-1) + O(-2) + O(0), d = 4  -->  True  (10256.00h)

---------------------------------------------
Tot_P2 O(-1)  in  Tot_P2 O(-1) + O(0) + O(-2)
---------------------------------------------
We test the vanishing of the PT series up to 6 orders beyond leading order in q and degree d <= 4:
Tot_P2 O(-1) + O(0) + O(-2), d = 1  -->  True  (0.06h)
Tot_P2 O(-1) + O(0) + O(-2), d = 2  -->  True  (7.69h)
Tot_P2 O(-1) + O(0) + O(-2), d = 3  -->  True  (393.64h)
Tot_P2 O(-1) + O(0) + O(-2), d = 4  -->  True  (10017.49h)

---------------------------------------------
Tot_P2 O(0)  in  Tot_P2 O(0) + O(-2) + O(-1)
---------------------------------------------
We test the vanishing of the PT series up to 6 orders beyond leading order in q and degree d <= 4:
Tot_P2 O(0) + O(-2) + O(-1), d = 1  -->  True  (0.06h)
Tot_P2 O(0) + O(-2) + O(-1), d = 2  -->  True  (9.24h)
Tot_P2 O(0) + O(-2) + O(-1), d = 3  -->  True  (405.34h)
Tot_P2 O(0) + O(-2) + O(-1), d = 4  -->  True  (9712.55h)

---------------------------------------------
Tot_P2 O(0)  in  Tot_P2 O(0) + O(-1) + O(-2)
---------------------------------------------
We test the vanishing of the PT series up to 6 orders beyond leading order in q and degree d <= 4:
Tot_P2 O(0) + O(-1) + O(-2), d = 1  -->  True  (0.06h)
Tot_P2 O(0) + O(-1) + O(-2), d = 2  -->  True  (8.76h)
Tot_P2 O(0) + O(-1) + O(-2), d = 3  -->  True  (419.47h)
Tot_P2 O(0) + O(-1) + O(-2), d = 4  -->  True  (9979.44h)

```
