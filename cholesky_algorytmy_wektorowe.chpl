module cholesky_wektorowy_algorytm {

  use cholesky_algorytmy_skalarne;

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

    for (all_cols, active_cols, later_cols) in iteruje_blok_kolumny ( wskazuje_wiersz_kolumne, wielkosc_bloku ) do {

    	// compute the Cholesky factor of the active diagonal block

    	czy_pozytywne_wartosci = cholesky_kolumnowa_bez_zrownoleglenia
    	                        ( A (active_cols, active_cols) );

    	if czy_pozytywne_wartosci && later_cols.length > 0 then {

    	  // compute the remainder of the active block column of L by a
    	  // block triangular solve realizing the equation
    	  //      L (later_cols, active_cols) =
    	  //                              L (later_cols, active_cols) *
    	  //                              L (active_cols, active_cols) ** (-T)

    	 // transposed_block_triangular_solve ( A (active_cols, active_cols),
    		//         		      A (later_cols, active_cols) );

    	// make rank wielkosc_bloku (outerproduct) modification to the remaining
    	// block rows and columns of  A, which become the Schur complement

    	 // symmetric_block_schur_complement (  A (later_cols, later_cols),
    		//			      A (later_cols, active_cols),
    		//			      wielkosc_bloku );

    	}
	    else if !czy_pozytywne_wartosci then return false;
      }
    return true;
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
