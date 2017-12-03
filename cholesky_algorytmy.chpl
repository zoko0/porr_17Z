module cholesky_algorytmy {

  /*
    algorytm wierszowy
    input: kwadrwatowa symetrczna macierzy - A
    output: nadpisana macierz A przez faktoryzacje cholesky - L
    zwraca: true  | jesli A (j, j) > 0.0
            false | jesli A (j, j) <= 0.0
  */

  proc cholesky_wierszowa_skalarna ( A : [] )  where ( A.domain.rank == 2 )  {

    const wskazuje_wiersz_kolumne = A.domain.dim (1);

    for j in wskazuje_wiersz_kolumne do {

      if A (j, j) > 0.0 then {

	       A (j, j)      = sqrt ( A (j, j) );
	       A (j, j+1..) /= A (j, j);

	       forall k in wskazuje_wiersz_kolumne (j+1..) do
	         A (k, k..) -= A(j, k..) * A (j, k);
      }
      else return false;
    }
    return true;
  }

  /*
    algorytm kolumnowy
    input: kwadrwatowa symetrczna macierzy - A
    output: nadpisana macierz A przez faktoryzacje cholesky - L
    zwraca: true  | jesli A (j, j) > 0.0
            false | jesli A (j, j) <= 0.0
  */
  proc cholesky_kolumnowa_skalarna ( A : [] )  where ( A.domain.rank == 2 ) {

    const wskazuje_wiersz_kolumne = A.domain.dim (1);

    for j in wskazuje_wiersz_kolumne do {

      if A (j, j) > 0.0 then {

	       A (j, j)       = sqrt ( A (j, j) );
	       A (j+1.., j ) /= A (j, j);


	       forall k in wskazuje_wiersz_kolumne (j+1..) do
	         A (k.., k) -= A(k.., j) * A (k, j);
         }
         else return false;
      }
    return true;
  }

}
