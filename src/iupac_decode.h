#ifndef IUPAC_DECODE_H
#define IUPAC_DECODE_H

inline void decodeIUPAC(char base, int &A, int &C, int &G, int &T, int &nAlleles) {
  switch(base) {
  case 'A': A += 2; nAlleles += 2; break;
  case 'C': C += 2; nAlleles += 2; break;
  case 'G': G += 2; nAlleles += 2; break;
  case 'T': T += 2; nAlleles += 2; break;
  case 'R': A += 1; G += 1; nAlleles += 2; break;
  case 'Y': C += 1; T += 1; nAlleles += 2; break;
  case 'S': C += 1; G += 1; nAlleles += 2; break;
  case 'W': A += 1; T += 1; nAlleles += 2; break;
  case 'K': G += 1; T += 1; nAlleles += 2; break;
  case 'M': A += 1; C += 1; nAlleles += 2; break;
  case 'B':
  case 'D':
  case 'H':
  case 'V':
  case 'N':
  case 'X':
  case '-':
  case '.':
    break;
  default:
    break;
  }
}

#endif
