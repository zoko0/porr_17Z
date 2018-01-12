module cholesky_algorytmy_wektorowe {

  use cholesky_algorytmy;

  const empty_range = 1..0;

  /*
    wersja kolumnowa cholesky
  */
  proc blokowy_cholesky ( A : [], wielkosc_bloku : int ) where ( A.domain.rank == 2 )

  {
    assert ( A.domain.dim (1) == A.domain.dim (2) && wielkosc_bloku > 0 );

    const wskazuje_wiersz_kolumne = A.domain.dim (1);
    var   czy_pozytywne_wartosci : bool;

    writeln ( "wielkosc_bloku: ", wielkosc_bloku );

    for (kolumny, kolumny_aktywne, kolumny_next) in iteruje_blok_kolumny ( wskazuje_wiersz_kolumne, wielkosc_bloku ) do {

    	// oblicz choleskiego dla przekatnej maciezry
    	czy_pozytywne_wartosci = cholesky_kolumnowa_bez_zrownoleglenia
    	                        ( A (kolumny_aktywne, kolumny_aktywne) );

    	if czy_pozytywne_wartosci && kolumny_next.length > 0 then {

    	  // compute the remainder of the active block column of L by a
    	  // block triangular solve realizing the equation
    	  //      L (kolumny_next, kolumny_aktywne) =
    	  //                              L (kolumny_next, kolumny_aktywne) *
    	  //                              L (kolumny_aktywne, kolumny_aktywne) ** (-T)

    	  rozwiaz_blok_transponowany ( A (kolumny_aktywne, kolumny_aktywne),
    		         		      A (kolumny_next, kolumny_aktywne) );

    	// make rank wielkosc_bloku (outerproduct) modification to the remaining
    	// block rows and columns of  A, which become the Schur complement

    	 // symetryczny_blok_uzupelnianie (  A (kolumny_next, kolumny_next),
    		//			      A (kolumny_next, kolumny_aktywne),
    		//			      wielkosc_bloku );

    	}
	    else if !czy_pozytywne_wartosci then return false;
      }
    return true;
  }

  /*
    Block Triangular Solve

    // ------------------------------------------------------
    // Solve the block equation
    //      L_przekatna_T = A_offdiag * L_przekatna^{-T}
    //           or
    //      L_przekatna_T^T = L_przekatna^{-1} A_offdiag^T
    // by triangular solve.
    // This code is specialized to a factorization case where
    // L and A are submatrices of a common larger matrix.
    // ------------------------------------------------------
  */
  proc rozwiaz_blok_transponowany ( L_przekatna : [], L_przekatna_T : [] ) {


    const kolumny_aktywne = L_przekatna.domain.dim(1);

    for (i,j) in L_przekatna_T.domain do {
      L_przekatna_T (i,j) -=
	+reduce [k in kolumny_aktywne (.. j-1)] L_przekatna_T (i,k) * L_przekatna (j,k);
      L_przekatna_T (i,j) = L_przekatna_T (i,j) / L_przekatna (j,j);
      }
  }

  /*
    Symmetric Block Outer Product_Modification
  */
  proc symetryczny_blok_uzupelnianie ( A : [] , L : [], block_size ) where ( A.domain.rank == 2 && L.domain.rank == 2) {
    for ( A_top_and_bottom_rows, A_top_rows, A_bottom_rows ) in iterated_block_column_partition (L.domain.dim (1), block_size) do {
    	symetryczna_przekatna_modyfikacja
    	             ( L (A_top_rows, ..),
    		       A (A_top_rows, A_top_rows) );

    	if A_bottom_rows.length > 0 then
    	  symetryczna_przekatna_T_modyfikacja
    	          ( L (A_top_and_bottom_rows, ..),
    		    A (A_bottom_rows, A_top_rows) );
    }
  }


  /*
    Symmetric Block Outer Product Modification for a single diagonal block

    // -----------------------------------------------------------
    // form diagonal block A (K,K) = A (K,K) - L (K,J) L^T (J,K)
    //                             = A (K,K) - L (K,J) L (K,J)^T
    // code is specialized to factorization case where L and A
    // are submatrices of a single larger matrix.
    // -----------------------------------------------------------
  */
 proc symetryczna_przekatna_modyfikacja ( L : [], A : [] ) {
   assert ( A.domain.dim (1) == A.domain.dim (2) && A.domain.dim (1) == L.domain.dim (1) );
   const A_diag_rows   = A.domain.dim (1),
         L_active_cols = L.domain.dim (2);

   forall i in A_diag_rows do
     forall j in A_diag_rows (..i) do
        A (i,j) -= + reduce [k in L_active_cols] L (i,k) * L (j,k);
   }


 /*
  Symmetric Block Outer Product Modification for a single offdiagonal block

  // -------------------------------------------------------------
  // Form a single offdiagonal block
  //       A (I,K) = A (I,K) - L (I,J) L^T (J,K)
  //               = A (I,K) - L (I,J) L (J,K)^T
  // This code is specialized to the triangular factorization case
  // where L and A are submatrices of a common larger matrix.
  // -------------------------------------------------------------
*/
 proc symetryczna_przekatna_T_modyfikacja ( L : [], A : [] ) {

   const L_active_cols  = L.domain.dim (2);

   forall (i,j) in A.domain do
     A (i,j) -= + reduce [k in L_active_cols] L (i,k) * L (j,k);
 }


  /*
    iteruje po wektorach ( kolumny lub wiersze macierzy )
  */
  iter iteruje_blok_czesc ( indeks_zasieg, wielkosc_bloku ) {
    var liczba_krokow_blokowych = ( indeks_zasieg + wielkosc_bloku -1 ) / wielkosc_bloku;

    // poczatek
    var dolny_blok = indeks_zasieg.low;

    // kontynuacja
    for blok_krok in 1 .. liczba_krokow_blokowych -1 do {
      yield dolny_blok .. #wielkosc_bloku;
      dolny_blok += wielkosc_bloku;
    }

    // koncowy blok
    yield dolny_blok .. indeks_zasieg.high;
  }


  /*
    zwraca wiersze w danej kolumnie
  */
  iter iteruje_blok_kolumny ( indeks_zasieg, wielkosc_bloku ) {
    var liczba_krokow_blokowych = ( indeks_zasieg.length + wielkosc_bloku - 1 ) / wielkosc_bloku;

    //poczatek
    var dolny_blok      = indeks_zasieg.low;
    var next_dolny_blok = dolny_blok + wielkosc_bloku;

    //kontynuacja
    for block_step in 1 .. liczba_krokow_blokowych - 1 do {
      yield ( dolny_blok      .. indeks_zasieg.high,
	      dolny_blok      .. #wielkosc_bloku,
	      next_dolny_blok .. indeks_zasieg.high );

      dolny_blok       = next_dolny_blok;
      next_dolny_blok += wielkosc_bloku;
    }

    // koncowy blok
    yield  ( dolny_blok .. indeks_zasieg.high,
	     dolny_blok .. indeks_zasieg.high,
	     empty_range );

  }

}
