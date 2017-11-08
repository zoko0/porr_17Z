module Hello{
  config const message="Hello, world!";
  config const numMessages=100;
  proc main() {
    forall msg in 1..numMessages do
      writeln(message +" (from iteration ", msg, " of ", numMessages, ")");
  }
}
