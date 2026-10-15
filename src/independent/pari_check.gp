\\ Independent verification in PARI/GP (no code shared with src/): the Dedekind zeta function is built directly from the
\\ field polynomial with lfuncreate, its zeros come from lfunzeros, zeta_K'(rho) from lfun(.,rho,1), and
\\   h*_{K,T}(t) = 2 Re Sum_{0<gamma<T} f(gamma/T) e^{i gamma t} / (rho zeta_K'(rho)),  f(y) = (1-y)cos(Pi y) + sin(Pi y)/Pi.
\\ Usage: gp -q pari_check.gp   (prints one line per test)
default(parisizemax, 2^31);
\p 60
f(y) = (1 - y)*cos(Pi*y) + sin(Pi*y)/Pi;
setup(pol, T) = {
  my(L = lfuncreate(pol), Z = lfunzeros(L, T), C);
  C = vector(#Z, j, my(rho = 1/2 + I*Z[j]); [Z[j], 1/(rho*lfun(L, rho, 1))]);
  [L, Z, C];
}
hstar(S, T, t) = 2*real(sum(j = 1, #S[3], my(g = S[3][j][1]); if(g > 0 && g < T, f(g/T)*S[3][j][2]*exp(I*g*t), 0)));
report(name, pol, T, tl) = {
  my(S = setup(pol, T), n = #select(g -> g > 0 && g < T, S[2]));
  printf("%s T=%d: lfunzeros finds %d zeros of zeta_K in (0,T)\n", name, T, n);
  for(i = 1, #tl, printf("   h*(%s) = %.15f\n", tl[i], hstar(S, T, eval(tl[i]))));
}
report("Q(sqrt-3)", x^2 + x + 1, 200, ["1560.51100774", "18558.8477306", "1560.6"]);
report("Q(i)",      x^2 + 1,     600, ["72.85", "-85.15", "158.854", "160.454"]);
report("Q(sqrt2)",  x^2 - 2,     150, ["-24.64", "2607.183", "195.906"]);
report("Q(sqrt3)",  x^2 - 3,     200, ["17.32", "631.336", "1746.876"]);
report("Q(sqrt5)",  x^2 - x - 1, 1000, ["1016.24860331", "722.734868673", "0"]);
quit;
