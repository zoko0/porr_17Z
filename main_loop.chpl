odule performance_cholesky_test {

  use Random, Time;
  use cholesky;

  /*****************************************************************************
    Konfigurajca:
    n - wymiar macierzy
    index_base - tylko i wylacznie do wyswietlania, powinno byc rowne 0 (TODO)
    print_matrix_details - czy drukowac macierze (funkcja print_lower_triangle)
    reproductible_output - false - oznacza, ze drukuje czasy
  *****************************************************************************/

  config const n = 1000;
  config const index_base = 0;

  config const print_matrix_details = false;
  config const reproducible_output = false;

  /*****************************************************************************
    Koniec konfiguracji
  *****************************************************************************/

  proc main {

    var Rand = new RandomStream ( real, seed = 314159) ;

    const mat_dom : domain (2) = { index_base .. #n, index_base .. #n };

    var A : [mat_dom] real,
        B : [mat_dom] real,
        L : [mat_dom] real;

    var positive_definite : bool;

    writeln ("Faktoryzacja Cholesky'ego");
    writeln ("   Wymiar macierzy: ", n);
    writeln ("");

    Rand.fillRandom (B);

    A = 0.0;

    forall (i,j) in mat_dom do
      A (i,j) = + reduce (  [k in mat_dom.dim (1) ]
    			    B (i, k) * B (j, k) );

    L = A; // algorytm nadpisuje macierz A


    writeln ("\n\n");
    writeln ("Zrownoleglenie wierszowe: " );

    var clock : Timer;

    clock.start ();

    positive_definite = scalar_row_major_outer_product_cholesky ( L );

    clock.stop ();
    if !reproducible_output then {
      writeln ( "Czas wykonania:    ", clock.elapsed () );
      writeln ( "Predkosc w megaflops: ",
		( (n**3) / 3.0 )  / (10.0**6 * clock.elapsed () ) );
    }

    forall j in mat_dom.dim (1) do
      forall i in j+1 .. mat_dom.dim(1).high do
	      L (i,j) = L (j,i);
    print_lower_triangle ( L );


    if positive_definite then
      check_factorization ( A, L );
    else
      writeln ("Niepowodzenie faktoryzacji");

    writeln ("\n\n");


    L = A;

    /*if print_matrix_details then {
      writeln ("test matrix");
      print_lower_triangle ( L );
    }*/

    writeln ("\n\n");
    writeln ("Zrownoleglenie kolumnowa: ");

    clock.clear ();
    clock.start ();

    positive_definite = scalar_column_major_outer_product_cholesky ( L );

    clock.stop ();

    if !reproducible_output then {
      writeln ( "Czas wykonania:    ", clock.elapsed () );
      writeln ( "Predkosc w megaflops: ",
		( (n**3) / 3.0 )  / (10.0**6 * clock.elapsed () ) );
    }

    print_lower_triangle ( L );

    if positive_definite then
      check_factorization ( A, L );
    else
      writeln ("Niepowodzenie faktoryzacji");

    delete Rand;
  }

  // sprawdza poprawnosc laczac L*L^T i porownujac z A
  proc check_factorization ( A : [], L : [] )
    where ( A.domain.rank == 2 ) {

    assert ( A.domain.dim (1) == A.domain.dim (2)  &&
	     L.domain.dim (1) == A.domain.dim (1)  &&
	     L.domain.dim (2) == A.domain.dim (2)
	     );

    const mat_dom  = A.domain,
          mat_rows = A.domain.dim(1),
          n        = A.domain.dim(1).length;

    const unit_roundoff = 2.0 ** (-53),
          gamma_n1      = (n * unit_roundoff) / (1.0 - n * unit_roundoff);

    var   max_ratio = 0.0;

    var   d : [mat_rows] real;

    for i in mat_rows do
      d (i) = sqrt ( A (i,i) );

    forall (i,j) in mat_dom with (ref max_ratio) do {
      const resid : real  = abs (A (i,j) -
		    + reduce ( [k in mat_dom.dim(1) (..min (i,j))]
			       L (i,k) * L (j,k) ) ) ;
      max_ratio = max ( max_ratio,
			resid * (1 - gamma_n1) /
			( gamma_n1 * d (i) * d (j) ) );
    }

    if max_ratio <= 1.0 then
      writeln ("Faktoryzacja wykonana z powodzeniem");
    else
      writeln ("Blad faktoryzacji. Kod:", max_ratio);
  }

  // wyswietla dolna czesc macierzy (gorna jest taka sama)
  proc print_lower_triangle ( L : [] ) {

    if print_matrix_details then
      for (i_row, i_col) in zip( L.domain.dim(1), L.domain.dim(2) ) do
	writeln (i_row, ":  ", L(i_row, ..i_col) );
  }
}
