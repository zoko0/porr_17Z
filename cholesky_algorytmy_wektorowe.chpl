module cholesky_algorytmy_wektorowe {

  use cholesky_algorytmy;

  const empty_range = 1..0;
  config const wielkosc_bloku = 13;


  /*
    wersja kolumnowa cholesky
  */
  proc blokowy_cholesky ( A : [], wielkosc_bloku : int ) where ( A.domain.rank == 2 )

  {
    assert ( A.domain.dim (1) == A.domain.dim (2) && wielkosc_bloku > 0 );

    const wskazuje_wiersz_kolumne = A.domain.dim (1);
    var   czy_pozytywne_wartosci : bool;

    writeln ( "wielkosc_bloku: ", wielkosc_bloku );

    for (kolumny, kolumny_aktywne, kolumny_next) in symmetric_reduced_matrix_2_by_2_block_partition ( wskazuje_wiersz_kolumne) do {

    	// oblicz choleskiego dla przekatnej maciezry
    	czy_pozytywne_wartosci = cholesky_kolumnowa_skalarna
    	                        ( A (kolumny_aktywne, kolumny_aktywne) );

    	if czy_pozytywne_wartosci && kolumny_next.length > 0 then {

    	  // oblicza pozostala czesc macierzy bloku L
    	  rozwiaz_blok_transponowany ( A (kolumny_aktywne, kolumny_aktywne), A (kolumny_next, kolumny_aktywne) );

    	  // modyfikacja poprzez symetryczne uzupelnianie bloku przy wykorzystaniu komplementacji Schur
    	  symetryczny_blok_uzupelnianie (  A (kolumny_next, kolumny_next),
    					      A (kolumny_next, kolumny_aktywne),
    					      wielkosc_bloku );

    	}
	    else if !czy_pozytywne_wartosci then return false;
    }
    return true;
  }

  /*
    Block Triangular Solve
    rozwiazuje blok o rownaniu
    //      L_przekatna_T = A_offdiag * L_przekatna^{-T}
    //           or
    //      L_przekatna_T^T = L_przekatna^{-1} A_offdiag^T
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
  proc symetryczny_blok_uzupelnianie ( A : [] , L : [], wielkosc_bloku ) where ( A.domain.rank == 2 && L.domain.rank == 2) {
    for ( A_top_and_bottom_rows, A_top_rows, A_bottom_rows ) in symmetric_reduced_matrix_2_by_2_block_partition (L.domain.dim (1)) do {
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
 */
   proc symetryczna_przekatna_T_modyfikacja ( L : [], A : [] ) {

     const L_active_cols  = L.domain.dim (2);

     forall (i,j) in A.domain do
       A (i,j) -= + reduce [k in L_active_cols] L (i,k) * L (j,k);
   }


  /*
    zwraca blok 2x2
  */
  iter symmetric_reduced_matrix_2_by_2_block_partition ( indeks_zasieg ) {

    for dolny_blok in indeks_zasieg by wielkosc_bloku do {
      var next_dolny_blok = dolny_blok + wielkosc_bloku;

      yield ( dolny_blok      .. indeks_zasieg.high,
	      dolny_blok      .. min ( next_dolny_blok - 1, indeks_zasieg.high ),
	      next_dolny_blok .. indeks_zasieg.high );
    }
  }

}
