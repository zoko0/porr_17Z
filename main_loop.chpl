module cholesky_test_wydajnosci {

  use Random, Time;
  use cholesky_algorytmy_skalarne;
  use cholesky_algorytmy_wektorowe;

  /*
    Konfiguracjca:
    n - wymiar macierzy
    bazowy_indeks - tylko i wylacznie do wyswietlania, powinno byc rowne 0
    wyswietlaj_macierz - czy drukowac macierze (funkcja wyswietl_dolny_trojkat_macierzy)
    nie_drukuj_czasow - false - oznacza, ze drukuje czasy
  */

  config const n = 10;
  config const bazowy_indeks = 0;

  config const wyswietlaj_macierz = false;
  config const nie_drukuj_czasow = false;

  /*
    Koniec konfiguracji
  */

  proc main {

    var Rand = new RandomStream ( real, seed = 314159) ;

    const zakres_macierzy : domain (2) = { bazowy_indeks .. #n, bazowy_indeks .. #n };

    var A : [zakres_macierzy] real,
        B : [zakres_macierzy] real,
        L : [zakres_macierzy] real;

    var czy_pozytywne_wartosci : bool;

    writeln ("Faktoryzacja Cholesky'ego");
    writeln ("Wymiar macierzy: ", n);
    writeln ("");

    Rand.fillRandom (B);

    A = 0.0;

    forall (i,j) in zakres_macierzy do
      A (i,j) = + reduce (  [k in zakres_macierzy.dim (1) ]
    			    B (i, k) * B (j, k) );

    L = A; // algorytm nadpisuje macierz A

    writeln ("\n\n");
    writeln ("Wersja wierszowa, bez zrownoleglenia: " );

    var zegar : Timer;

    zegar.start ();

    czy_pozytywne_wartosci = cholesky_wierszowa_bez_zrownoleglenia ( L );

    zegar.stop ();
    if !nie_drukuj_czasow then {
      writeln ( "Czas wykonania:    ", zegar.elapsed () );
      writeln ( "Predkosc w megaflops: ",	( (n**3) / 3.0 )  / (10.0**6 * zegar.elapsed () ) );
    }

    forall j in zakres_macierzy.dim (1) do
      forall i in j+1 .. zakres_macierzy.dim(1).high do
	      L (i,j) = L (j,i);
    wyswietl_dolny_trojkat_macierzy ( L );


    if czy_pozytywne_wartosci then
      sprawdzenie_poprawnosci ( A, L );
    else
      writeln ("Niepowodzenie faktoryzacji");

    writeln ("\n\n");

    writeln ("\n\n");
    writeln ("Wersja wierszowa, z zrownolegleniem: " );

    var zegar : Timer;

    zegar.start ();

    czy_pozytywne_wartosci = cholesky_wierszowa_skalarna ( L );

    zegar.stop ();
    if !nie_drukuj_czasow then {
      writeln ( "Czas wykonania:    ", zegar.elapsed () );
      writeln ( "Predkosc w megaflops: ",	( (n**3) / 3.0 )  / (10.0**6 * zegar.elapsed () ) );
    }

    forall j in zakres_macierzy.dim (1) do
      forall i in j+1 .. zakres_macierzy.dim(1).high do
	      L (i,j) = L (j,i);
    wyswietl_dolny_trojkat_macierzy ( L );


    if czy_pozytywne_wartosci then
      sprawdzenie_poprawnosci ( A, L );
    else
      writeln ("Niepowodzenie faktoryzacji");

    writeln ("\n\n");


    L = A;

    /*if wyswietlaj_macierz then {
      writeln ("test matrix");
      wyswietl_dolny_trojkat_macierzy ( L );
    }*/

    writeln ("\n\n");
    writeln ("Wersja kolumnowa, bez zrownoleglenia: " );

    zegar.clear ();
    zegar.start ();

    czy_pozytywne_wartosci = cholesky_kolumnowa_bez_zrownoleglenia ( L );

    zegar.stop ();

    if !nie_drukuj_czasow then {
      writeln ( "Czas wykonania:    ", zegar.elapsed () );
      writeln ( "Predkosc w megaflops: ",
    ( (n**3) / 3.0 )  / (10.0**6 * zegar.elapsed () ) );
    }

    wyswietl_dolny_trojkat_macierzy ( L );

    if czy_pozytywne_wartosci then
      sprawdzenie_poprawnosci ( A, L );
    else
      writeln ("Niepowodzenie faktoryzacji");

    delete Rand;
    }


    writeln ("\n\n");
    writeln ("Wersja kolumnowa, z zrownolegleniem: " );

    zegar.clear ();
    zegar.start ();

    czy_pozytywne_wartosci = cholesky_kolumnowa_skalarna ( L );

    zegar.stop ();

    if !nie_drukuj_czasow then {
      writeln ( "Czas wykonania:    ", zegar.elapsed () );
      writeln ( "Predkosc w megaflops: ",
		( (n**3) / 3.0 )  / (10.0**6 * zegar.elapsed () ) );
    }

    wyswietl_dolny_trojkat_macierzy ( L );

    if czy_pozytywne_wartosci then
      sprawdzenie_poprawnosci ( A, L );
    else
      writeln ("Niepowodzenie faktoryzacji");

    delete Rand;
  }


  // sprawdza poprawnosc laczac L*L^T i porownujac z A
  proc sprawdzenie_poprawnosci ( A : [], L : [] )
    where ( A.domain.rank == 2 ) {

    assert ( A.domain.dim (1) == A.domain.dim (2)  &&
	     L.domain.dim (1) == A.domain.dim (1)  &&
	     L.domain.dim (2) == A.domain.dim (2)
	     );

    const zakres_macierzy  = A.domain,
          wiersze_macierzy = A.domain.dim(1),
          n        = A.domain.dim(1).length;

    const jednostka_zaokroglenia = 2.0 ** (-53),
          gamma_n1      = (n * jednostka_zaokroglenia) / (1.0 - n * jednostka_zaokroglenia);

    var   max_ratio = 0.0;

    var   d : [wiersze_macierzy] real;

    for i in wiersze_macierzy do
      d (i) = sqrt ( A (i,i) );

    forall (i,j) in zakres_macierzy with (ref max_ratio) do {
      const resid : real  = abs (A (i,j) -
		    + reduce ( [k in zakres_macierzy.dim(1) (..min (i,j))]
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
  proc wyswietl_dolny_trojkat_macierzy ( L : [] ) {

    if wyswietlaj_macierz then
      for (i_row, i_col) in zip( L.domain.dim(1), L.domain.dim(2) ) do
	writeln (i_row, ":  ", L(i_row, ..i_col) );
  }
}
