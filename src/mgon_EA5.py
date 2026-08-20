r"""E[A_5] for regular m-gons at 50 digits (route P), and a search for a closed form.
Reported honestly: the PSLQ basis tried is printed whether or not it hits."""
import json
from mpmath import mp, mpf, nstr, cos, sin, pi, pslq, identify
from polygon_exact import T_and_moments, regular_polygon_mp
mp.dps = 50

rows = []
print("m    E[A_3]                                    E[A_5]                                    E[A_6]-(3E[A_5]-5E[A_3])")
for m in list(range(3, 17)) + [20, 24, 30, 40, 60]:
    r = T_and_moments(regular_polygon_mp(m, 50), nmax=7, dps=50, nodes=9)
    a3, a5, a6 = r['E[A_3]'], r['E[A_5]'], r['E[A_6]']
    rows.append(dict(m=m, EA3=nstr(a3, 40), EA4=nstr(r['E[A_4]'], 40), EA5=nstr(a5, 40),
                     EA6=nstr(a6, 40), resid_I=nstr(r['E[A_4]'] - 2 * a3, 5),
                     resid_Ip=nstr(a6 - (3 * a5 - 5 * a3), 5)))
    print(f"{m:3d}  {nstr(a3,30):40s}  {nstr(a5,30):40s}  {nstr(a6-(3*a5-5*a3),5)}", flush=True)

print("\n--- closed-form search for E[A_5](regular m-gon) ---")
print("Alikoski's E[A_3] has the shape (quadratic in cos w)/(36 m^2 sin^2 w), w = 2pi/m.")
print("Trying the analogous ansatz E[A_5] = (a0 + a1 C + a2 C^2 + a3 C^3)/(D m^2 S^2)+... ")
for m in [4, 6]:
    r = T_and_moments(regular_polygon_mp(m, 50), nmax=7, dps=50, nodes=9)
    a5 = r['E[A_5]']
    print(f"  m={m}: E[A_5] = {nstr(a5,40)}   identify -> {identify(a5)}")
mdisk = T_and_moments(regular_polygon_mp(2000, 50), nmax=7, dps=50, nodes=9)
print(f"  m=2000 (disk proxy): E[A_5] = {nstr(mdisk['E[A_5]'],30)}")
print(f"  PSLQ basis [E[A_5], 1, 1/pi^2, 1/pi^4]:",
      pslq([mdisk['E[A_5]'], mpf(1), 1 / pi ** 2, 1 / pi ** 4], maxcoeff=10**10, maxsteps=10**6, tol=mpf(10)**-20))
json.dump(rows, open("../results/C_mgon_EA5.json", "w"), indent=1)
print("wrote ../results/C_mgon_EA5.json")
