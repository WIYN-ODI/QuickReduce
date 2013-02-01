awk '
{
    x=($1+$2)/2;
    if (x>260 && x<890) print x, $3+ 0.00555 - (940-x)*3.2e-5;
}
' xxx.plot >odi_r.profile


  