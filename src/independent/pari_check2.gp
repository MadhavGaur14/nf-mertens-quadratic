\\ Diagnose the zero-count discrepancy between lfunzeros(zeta_K) and our certified counts.
\\ (a) zeros of zeta and L(chi_D) separately; (b) zeta_K with a finer subdivision (3rd argument of lfunzeros);
\\ (c) list the zeros found by (a) but not by plain lfunzeros(zeta_K); (d) h* from the union (a).
default(parisizemax, 2^31);
\p 60
f(y) = (1 - y)*cos(Pi*y) + sin(Pi*y)/Pi;
hs(Z, LK, T, t) = 2*real(sum(j = 1, #Z, my(g = Z[j], rho = 1/2 + I*g); if(g > 0 && g < T, f(g/T)*exp(I*g*t)/(rho*lfun(LK, rho, 1)), 0)));
run(name, pol, D, T, tl) = {
  my(LK = lfuncreate(pol), Zz = lfunzeros(lfuncreate(1), T), ZL = lfunzeros(lfuncreate(D), T), ZK = lfunzeros(LK, T),
     ZK64 = lfunzeros(LK, T, 64), U, miss);
  U = vecsort(concat(Zz, ZL));
  printf("%s T=%d: zeta %d + L(chi_%d) %d = %d | zeta_K plain %d | zeta_K divz=64 %d\n",
         name, T, #Zz, D, #ZL, #U, #ZK, #ZK64);
  miss = select(g -> vecmin(apply(h -> abs(h - g), ZK)) > 1e-8, U);
  printf("   zeros missed by plain lfunzeros(zeta_K): %s\n", apply(g -> round(g*10^6)/10^6., miss));
  printf("   distance from each missed zero to its nearest neighbour: %s\n",
         apply(g -> my(o = select(h -> h != g, U)); round(vecmin(apply(h -> abs(h - g), o))*10^6)/10^6., miss));
  for(i = 1, #tl, printf("   h* from union (a), t=%s: %.15f\n", tl[i], hs(U, LK, T, eval(tl[i]))));
}
run("Q(sqrt-3)", x^2 + x + 1, -3, 200, ["1560.51100774", "18558.8477306"]);
run("Q(sqrt2)",  x^2 - 2,      8, 150, ["-24.64", "2607.183"]);
quit;
