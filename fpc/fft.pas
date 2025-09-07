(* fast Fourier Transform                                             *)
(* Cooley-Tukey factorization of the discrete Fourier Transform       *)
(* from:                                                              *)
(*   https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm *)
(* Simplifications:                                                   *)
(* - 512-element input data array                                     *)
(* - Real-only input and output                                       *)
(*   - Output result contains magnitude data only; no phase data      *)
type
  Myfloat = Real;

procedure fft(Var indata : array of Myfloat; Var fftmagn : array of Myfloat);
var
  (* Cooley-Tukey Wiki variables                                      *)
  (*                                                                  *)
  (* Argument parameters (above):                                     *)
  (*                                                                  *)
  (*   indata       a, input, real only                               *)
  (*   fftmagn      Magnitudes of A, output, real only                *)
  (*                                                                  *)
  (* Local variables:                                                 *)
  (*                                                                  *)
  Areal  : array [0..511] of Myfloat; (* A, DFT, real components      *)
  Aimag  : array [0..511] of Myfloat; (* A, imag(inary) components    *)
  Wreal  : Myfloat = -1.0;            (* omega, real                  *)
  Wimag  : Myfloat = 0.0;             (* omega, imag(inary)           *)
  WMreal : Myfloat = -1.0;            (* omega-sub-m, real            *)
  WMimag : Myfloat = 0.0;             (* omega-sub-m, imag            *)
  Treal  : Myfloat = 0.0;             (* t, real                      *)
  Timag  : Myfloat = 0.0;             (* t, imag                      *)
  Ureal  : Myfloat = 0.0;             (* u, real (imag not needed)    *)
  m      : Integer = 0;               (* m, 2**s                      *)
  halfm  : Integer = 0;               (* m/2                          *)
  kstep  : Integer = 0;               (* k, steps through A by m      *)
  jincr  : Integer = 0;               (* j, increments from 0 to m/2  *)
  kj     : Integer = 0;               (* k + j                        *)
  kjhm   : Integer = 0;               (* k + j + m/2                  *)

  (* The rest are used in the bit reversal algorithm *)

  ibrs : array [0..511] of Integer; (* Source indices of copied value *)
  itrg : Integer;                   (* Target index into ibrs         *)
  isrc : Integer;                   (* Source index to which to add   *)
  iadd : Integer;                   (* index difference               *)
  ilim : Integer;                   (* loop contol                    *)
begin

  (* Start Cooley-Tukey:  bit-reversal; divide and conquer *)

  (* bit-reverse-copy(a, A); indices stored in array ibrs[0..511]     *)
  ibrs[0] := 0;   (* Seed first index; reversal of index 0 is 0       *)
  iadd := 512;    (* Seed difference to add                           *)
  itrg := 1;       (* Start at ibrs[1]                                *)
  while iadd > 1 do begin   (* Loop until next halved iadd would be 0 *)
    iadd := iadd shr 1;     (* halve iadd i.e. =>256=>...=>1          *)
    ilim := 512 Div iadd;   (* double limit =>2=>4=>...=>512          *)
    isrc := 0;              (* Start source index at ibrs[0]          *)
    while itrg < ilim do begin            (* Loop from itrg to ilim-1 *)
      ibrs[itrg] := ibrs[isrc] + iadd;       (* Save target index     *)
      Areal[itrg] := indata[ibrs[itrg]];     (* Copy real data        *)
      Aimag[itrg] := 0.0;                    (* Zero out imag data    *)
      isrc := isrc + 1;                      (* Next source index     *)
      itrg := itrg + 1;                      (* Next target index     *)
    end;
  end;

  (* Divide and conquer *)

  (* wM <= exp(-2PIi/m) : part 1 of 2; initialization when m is 2 *)
  WMreal := -1.0; (*  cos(-2PI/m) =  cos(-2PI/2) =  cos(-PI) = -1.0 *)
  WMimag := 0.0;  (* isin(-2PI/m) = isin(-2PI/2) = isin(-PI) =  0.0 *)

  (* for s = 1 to log(n) : base-2 logarithm *)
  halfm := 1;
  while halfm < 512 do begin
    (* m <= 2^s *)
    m := halfm shl 1;
    (* for k = 0 to n-1 by m : part 1 of 2 *)
    kstep := 0;
    while kstep < 512 do begin
      (* w <= 1 *)
      Wreal := 1.0;
      Wimag := 0.0;
      (* for j = 0 to m/2-1 *)
      for jincr := 0 to (halfm-1) do begin
        kj := kstep + jincr;
        kjhm := kj + halfm;
        (* t <= w * A[k+j+m/2] *)
        (* u <= A[k+j] : not needed because A[k+j] is last *)
        Treal := (Wreal * Areal[kjhm]) - (Wimag * Aimag[kjhm]);
        Timag := (Wreal * Aimag[kjhm]) + (Wimag * Areal[kjhm]);
        (* A[k+j+m/2] <= u - t : use A....[kj] instead of u *)
        Areal[kjhm] := Areal[kj] - Treal;
        Aimag[kjhm] := Aimag[kj] - Timag;
        (* A[k+j] <= u + t *)
        Areal[kj] := Areal[kj] + Treal;
        Aimag[kj] := Aimag[kj] + Timag;
        (* w <= w * wM : Ureal holds the new value of Wreal until the *)
        (*               calculation of Wimag using Wreal is complete *)
        Ureal := (Wreal * WMreal) - (Wimag * WMimag);
        Wimag := (Wreal * WMimag) + (Wimag * WMreal);
        Wreal := Ureal;
      end;
      (* for k = 0 to n=1 by m : part 2 of 2 *)
      kstep := kstep + m;
    end;
    (* Prepare for next pass through the outer while loop             *)
    (* (i) double halfm (m/2), which will double m on next pass start *)
    halfm := m;
    (* (ii) wM <= exp(-2PIi/m) : part 2 of 2, with doubled m          *)
    WMimag := sqrt((1.0 - WMreal) / 2.0);
    WMreal := sqrt((1.0 + WMreal) / 2.0);
  end;                                         (* End of Cooley-Tukey *)

  (* Convert complex DFT to magnitudes N.B. not part of Cooley-Tukey  *)
  for kj := 0 to 511 do begin
    fftmagn[kj] := sqrt((Areal[kj] * Areal[kj]) + (Aimag[kj] * Aimag[kj]));
  end;
end;

var
  indata : array [0..511] of Myfloat;
  fftmagn : array [0..511] of Myfloat;
  icount : Integer = 0;
  twopi : Myfloat = 6.283185307179586;
  F : Text;
begin
  icount := 0;

(* ****************************************************************** *)
(*                                                                    *)
(* Usage:                                                             *)
(*                                                                    *)
(*   (i) Generate test data internally; write input data to file x.x  *)
(*       ./fft > x.x                                                  *)
(*                                                                    *)
(*   (ii) Read 512 time-domain data from file, one datum per line     *)
(*    ./fft ../default_data.dat                                       *)
(*                                                                    *)
(* Output format:                                                     *)
(* - [offset,value] pairs                                             *)
(* - one pair per line                                                *)
(* - space delimited                                                  *)
(* - offset range is [0,N-1]                                          *)
(* - First N pairs are time-domain input data,                        *)
(* - followed by N pairs of frequency-domain data (FFT)               *)
(* - See example below; ### and text that follow are comments         *)
(*                                                                    *)
(* 0  1.0000000000000000E+000   ### Time-domain data, 1st input value *)
(* 1  1.6872781844821536E+000   ### 2nd input value                   *)
(* ...                          ### etc.                              *)
(* 511 -2.7306462210913929E-001 ### Last input value                  *)
(* 0  2.8421709430404007E-014   ### Result data, first FFT value      *)
(* 1  1.8012239161397711E-014   ### 2nd FFT value                     *)
(* ...                          ### etc.                              *)
(* 511  1.8012239161392905E-014 ### Last FFT value                    *)
(*                                                                    *)
(* ****************************************************************** *)

  if ParamCount > 0 then begin
    (* Get file name from command line parameter 1 *)
    WriteLn(StdErr, 'Reading data from file [', ParamStr(1), ']...');
    Assign(F, ParamStr(1));
    Reset(F);
    while icount < 512 do begin
      ReadLn(F,indata[icount]);
      WriteLn(icount,' ',indata[icount]);
      icount := icount + 1;
    end;
  end
  else begin
    (* No command line parameters; generate data *)
    WriteLn(StdErr, 'Generating test data...');
    while icount < 512 do begin
      indata[icount] := (10.0 * sin(icount * twopi / (512.0 /  8.0)))
                      + ( 1.0 * cos(icount * twopi / (512.0 / 64.0)));
      WriteLn(icount,' ',indata[icount]);
      icount := icount + 1;
    end;
  end;

  fft(indata, fftmagn);
  icount := 0;
  while icount < 512 do begin
    WriteLn(icount,' ',fftmagn[icount]);
    icount := icount + 1;
  end;
end.
