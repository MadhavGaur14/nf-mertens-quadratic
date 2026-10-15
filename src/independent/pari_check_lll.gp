\\ Independent PARI/GP check of the lattice-reduction witnesses (paper Table 2), T = 1000.
\\ Zeros of zeta and L(chi_D) separately (lfunzeros on zeta_K misses close pairs; see pari_check2.gp) at \p 130,
\\ zeta_K'(rho) from the Dedekind zeta function of the field polynomial. The zeros are also written to a text file
\\ for a digit-by-digit comparison with the Arb-certified 100-digit zeros (compare_zeros.py).
default(parisizemax, 2^31);
\p 130
f(y) = (1 - y)*cos(Pi*y) + sin(Pi*y)/Pi;
run(name, pol, D, T, tl, zfile) = {
  my(LK = lfuncreate(pol), U, C);
  U = vecsort(concat(lfunzeros(lfuncreate(1), T), lfunzeros(lfuncreate(D), T)));
  C = vector(#U, j, my(rho = 1/2 + I*U[j]); f(U[j]/T)/(rho*lfun(LK, rho, 1)));
  printf("%s T=%d: %d zeros (zeta + L(chi_%d))\n", name, T, #U, D);
  for(j = 1, #U, write(zfile, Strprintf("%.110f", U[j])));
  for(i = 1, #tl, my(t = eval(tl[i])); printf("   h*(%s) = %.18f\n", tl[i], 2*real(sum(j = 1, #U, C[j]*exp(I*U[j]*t)))));
}
\\ (gp ends a statement at each newline, so each call must be on one line)
run("Q(sqrt-3)", x^2 + x + 1, -3, 1000, ["165562943829466520552313925650444335142025426569159622719799633443383327010312592", "3431377126899343121084375092595828581527345096545547547377"], "pari_zeros_Qsqrt-3_T1000.txt");
run("Q(sqrt5)", x^2 - x - 1, 5, 1000, ["6701024002746072686211658750800192982264466383201014802804", "-3826144185540628524340499373955955297176444083490968202937", "1016.24860331"], "pari_zeros_Qsqrt5_T1000.txt");
quit;
