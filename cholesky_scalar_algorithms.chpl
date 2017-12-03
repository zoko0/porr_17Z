module cholesky_scalar_algorithms {

  /*
    algorytm wierszowy
    input: kwadrwatowa symetrczna macierzy - A
    output: nadpisana macierz A przez faktoryzacje cholesky - L
    zwraca: true  | jesli A (j, j) > 0.0
            false | jesli A (j, j) <= 0.0
  */

  /*
  (TODO) ??? 
    musimy rozkminic co oznacza:
    where ( A.domain.rank == 2 ) na poczatku obu funkcji
    zmienna A_rc_indices
    konkretnie jak dzialają oba algorytmy (kazda linijka od poczatku pierwszego fora)!!

  */

  proc scalar_row_major_outer_product_cholesky ( A : [] )  where ( A.domain.rank == 2 )  {

    const A_rc_indices = A.domain.dim (1);  // indices of either row or column

    for j in A_rc_indices do {

      if A (j, j) > 0.0 then {

	       A (j, j)      = sqrt ( A (j, j) );
	       A (j, j+1..) /= A (j, j);

	       forall k in A_rc_indices (j+1..) do
	         A (k, k..) -= A(j, k..) * A (j, k);
      }
      else return false; // error return if matrix is not positive sdefinite
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
  proc scalar_column_major_outer_product_cholesky ( A : [] )  where ( A.domain.rank == 2 ) {

    const A_rc_indices = A.domain.dim (1);  // row and column indices of A

    for j in A_rc_indices do {

      if A (j, j) > 0.0 then {

	       A (j, j)       = sqrt ( A (j, j) );
	       A (j+1.., j ) /= A (j, j);


	       forall k in A_rc_indices (j+1..) do
	         A (k.., k) -= A(k.., j) * A (k, j);
         }
         else return false;
      }
    return true;
  }

}
