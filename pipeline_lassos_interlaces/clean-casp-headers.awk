!($1 == "PFRMAT" || $1 == "TARGET" || $1 == "PARENT" || $1 == "AUTHOR" || $1 == "REMARK" || $1 == "METHOD") && $1 != "MODEL" && $1 != "END"
$1 == "MODEL" {if (no != $2) {print $0; no = $2;}}
$1 == "TER" {a=1;}
$1 == "END" {if (a!=1) {print "TER";} a=0;}
END {print "END";}