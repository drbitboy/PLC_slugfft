program br(output);
var
  i,j,ioffset,iadd,ilim : Integer;
  ibr : array [0..511] of Integer;
  pfx : string;
begin
  ibr[0] := 0;
  iadd := 512;
  ioffset := 1;
  i := 1;
  while iadd > 1 do begin
    iadd := iadd shr 1;
    ilim := ioffset shl 1;
    j := i - ioffset;
    while i < ilim do begin
      ibr[i] := ibr[i - ioffset] + iadd;
      i := i + 1;
      j := j + 1;
    end;
    ioffset := ioffset shl 1;
  end;

  pfx := '=(';
  for i := 0 to 511 do begin
    if ((i mod 16)=0) then begin
      writeln('');
    end;
    write(pfx,ibr[i]);
    pfx := ',';
  end;
  writeln(');');
end.
