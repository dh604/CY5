# File to test a O(-1)+O(-1) chain of length r.

r = 3
G = minus_one_minus_one_chain(r)
b = curve_class(G, Edge(1, 2))
c = curve_class(G, Edge(2, 3))
d = curve_class(G, Edge(3 ,4))

max_genus = 2

get_Omega_beta(G, [6*b], max_genus; check_predictions=true)
# All predictions hold up to genus 2.

get_Omega_beta(G, [3*b + 3*c], max_genus; check_predictions=true)
# All predictions hold up to genus 2.

get_Omega_beta(G, [4*b + 2*c], max_genus; check_predictions=true)
# All predictions hold up to genus 2.

get_Omega_beta(G, [2*(b+c+d)], max_genus; check_predictions=true)
# All predictions hold up to genus 2.

get_Omega_beta(G, [2*b+c+d], max_genus; check_predictions=true)
# All predictions hold up to genus 2.

get_Omega_beta(G, [b+2*c+d], max_genus; check_predictions=true)
# All predictions hold up to genus 2.