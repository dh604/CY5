# Local $\mathbb{P}^1 \times \mathbb{P}^1$

The purpose of the file `local_P1xP1.sage` is to probe the vanishing claimed for the two local surfaces
```math
\mathrm{Tot}_{\mathbb{P}^1 \times \mathbb{P}^1}\, \mathcal{O}(-1,-1) \oplus \mathcal{O}(-1,-1) \oplus \mathcal{O}
```
and
```math
\mathrm{Tot}_{\mathbb{P}^1 \times \mathbb{P}^1}\, \mathcal{O}(-2,0) \oplus \mathcal{O}(0,-2) \oplus \mathcal{O}
```
in [main_paper; Evidence 5.7](@cite).

For the first, the claim is that the PT series of any threefold inside the fivefold vanishes. We test this for the two distinct choices of threefolds that arise as some fixed locus of a fibre-wise $\mathbb{C}^{\times}_q$-action on the fivefold.

For the second local surface the claim is that all membrane indices but $\Omega_{(1,0)}$ and $\Omega_{(0,1)}$ vanish. For the two non-zero indices we also provide a conjectural formula. We check the formula for two choices of a threefold in the fivefold.

The script `local_P1xP1.sage` can be used to probe this claim. It creates the following output:
```
The PT invariants of the following targets all cojecturally vanish.
We test this conjecture up to 6 orders beyond leading order in q in q and bi-degrees (d1,d2) with d1 + d2 <= 4:

--------------------------------------------------------------
Tot_P1xP1 O(-1,-1)  in  Tot_P1xP1 O(-1,-1) + O(-1,-1) + O(0,0)
--------------------------------------------------------------
(d1,d2) = (0,0)  -->  True  (0.00h)
(d1,d2) = (0,1)  -->  True  (0.00h)
(d1,d2) = (0,2)  -->  True  (0.01h)
(d1,d2) = (0,3)  -->  True  (0.72h)
(d1,d2) = (0,4)  -->  True  (18.85h)
(d1,d2) = (1,0)  -->  True  (0.00h)
(d1,d2) = (1,1)  -->  True  (0.00h)
(d1,d2) = (1,2)  -->  True  (1.13h)
(d1,d2) = (1,3)  -->  True  (34.30h)
(d1,d2) = (2,0)  -->  True  (0.01h)
(d1,d2) = (2,1)  -->  True  (0.13h)
(d1,d2) = (2,2)  -->  True  (17.92h)
(d1,d2) = (3,0)  -->  True  (0.17h)
(d1,d2) = (3,1)  -->  True  (1.75h)
(d1,d2) = (4,0)  -->  True  (2.40h)

--------------------------------------------------------------
Tot_P1xP1 O(0,0)  in  Tot_P1xP1 O(0,0) + O(-1,-1) + O(-1,-1)
--------------------------------------------------------------
(d1,d2) = (0,0)  -->  True  (0.00h)
(d1,d2) = (0,1)  -->  True  (0.00h)
(d1,d2) = (0,2)  -->  True  (0.01h)
(d1,d2) = (0,3)  -->  True  (0.17h)
(d1,d2) = (0,4)  -->  True  (2.51h)
(d1,d2) = (1,0)  -->  True  (0.00h)
(d1,d2) = (1,1)  -->  True  (0.00h)
(d1,d2) = (1,2)  -->  True  (0.15h)
(d1,d2) = (1,3)  -->  True  (1.68h)
(d1,d2) = (2,0)  -->  True  (0.01h)
(d1,d2) = (2,1)  -->  True  (0.15h)
(d1,d2) = (2,2)  -->  True  (3.04h)
(d1,d2) = (3,0)  -->  True  (0.17h)
(d1,d2) = (3,1)  -->  True  (1.81h)
(d1,d2) = (4,0)  -->  True  (2.47h)

Start setting up the conjectural formulas for the next PT series to analyse...
Done! This took 11.59h.
Now we check the formula up to 6 orders beyond leading order in q and bi-degrees (d1,d2) with d1 + d2 <= 4:

-----------------------------------------------------------
Tot_P1xP1 O(-2,0)  in  Tot_P1xP1 O(-2,0) + O(0,-2) + O(0,0)
-----------------------------------------------------------
(d1,d2) = (0,0)  -->  True  (0.00h)
(d1,d2) = (0,1)  -->  True  (0.00h)
(d1,d2) = (0,2)  -->  True  (0.01h)
(d1,d2) = (0,3)  -->  True  (0.18h)
(d1,d2) = (0,4)  -->  True  (2.75h)
(d1,d2) = (1,0)  -->  True  (0.00h)
(d1,d2) = (1,1)  -->  True  (0.00h)
(d1,d2) = (1,2)  -->  True  (0.15h)
(d1,d2) = (1,3)  -->  True  (1.78h)
(d1,d2) = (2,0)  -->  True  (0.01h)
(d1,d2) = (2,1)  -->  True  (0.15h)
(d1,d2) = (2,2)  -->  True  (3.06h)
(d1,d2) = (3,0)  -->  True  (0.44h)
(d1,d2) = (3,1)  -->  True  (1.68h)
(d1,d2) = (4,0)  -->  True  (2.23h)

-----------------------------------------------------------
Tot_P1xP1 O(0,0)  in  Tot_P1xP1 O(0,0) + O(0,-2) + O(-2,0)
-----------------------------------------------------------
(d1,d2) = (0,0)  -->  True  (0.00h)
(d1,d2) = (0,1)  -->  True  (0.00h)
(d1,d2) = (0,2)  -->  True  (0.01h)
(d1,d2) = (0,3)  -->  True  (0.53h)
(d1,d2) = (0,4)  -->  True  (8.57h)
(d1,d2) = (1,0)  -->  True  (0.00h)
(d1,d2) = (1,1)  -->  True  (0.00h)
(d1,d2) = (1,2)  -->  True  (0.14h)
(d1,d2) = (1,3)  -->  True  (1.93h)
(d1,d2) = (2,0)  -->  True  (0.01h)
(d1,d2) = (2,1)  -->  True  (0.15h)
(d1,d2) = (2,2)  -->  True  (3.65h)
(d1,d2) = (3,0)  -->  True  (0.21h)
(d1,d2) = (3,1)  -->  True  (2.04h)
(d1,d2) = (4,0)  -->  True  (12.37h)
```
These checks confirm a part of [main_paper; Evidence 5.7](@cite). For the other part see [here](local_P2.md).
