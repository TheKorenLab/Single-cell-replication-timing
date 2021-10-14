# /usr/bin/awk -f

BEGIN { OFS = "\t" }

{ split($1, full_line, ":");
  id=full_line[1];
  for (i=2; i<=7; i++) {
      id=id":"full_line[i];
  }
}

{ if (length(full_line)==8) {
     barcode=full_line[8];
     $1 = id;
     { print $0"\tCB:Z:"barcode };
} else
      { print $0 }
}
