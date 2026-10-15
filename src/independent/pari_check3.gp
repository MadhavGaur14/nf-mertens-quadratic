\\ Independent PARI/GP verification, remaining witnesses. Zeros of zeta and L(chi_D) computed separately with
\\ lfunzeros (lfunzeros on zeta_K itself misses close zeta/L pairs at the default step; see pari_check2.gp);
\\ zeta_K'(rho) from the Dedekind zeta function built from the field polynomial.
default(parisizemax, 2^31);
\p 60
f(y) = (1 - y)*cos(Pi*y) + sin(Pi*y)/Pi;
run(name, pol, D, T, tl) = {
  my(LK = lfuncreate(pol), U = vecsort(concat(lfunzeros(lfuncreate(1), T), lfunzeros(lfuncreate(D), T))), C);
  C = vector(#U, j, my(rho = 1/2 + I*U[j]); f(U[j]/T)/(rho*lfun(LK, rho, 1)));
  printf("%s T=%d: %d zeros (zeta + L(chi_%d))\n", name, T, #U, D);
  for(i = 1, #tl, my(t = eval(tl[i])); printf("   h*(%s) = %.15f\n", tl[i], 2*real(sum(j = 1, #U, C[j]*exp(I*U[j]*t)))));
}
run("Q(i)",     x^2 + 1,     -4, 600,  ["72.85", "-85.15", "158.854", "160.454"]);
run("Q(sqrt3)", x^2 - 3,     12, 200,  ["17.32", "631.336", "1746.876"]);
run("Q(sqrt5)", x^2 - x - 1,  5, 1000, ["1016.24860331", "722.734868673"]);
quit;
