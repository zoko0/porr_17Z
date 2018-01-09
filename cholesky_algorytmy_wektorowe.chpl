module cholesky_wektorowy_algorytm {

  /*
    iteruje po wektorach ( kolumny lub wiersze macierzy )
  */
  iter blok_czesc ( indeks_zasieg, wielkosc_bloku ) {
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


}
