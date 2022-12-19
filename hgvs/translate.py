import hgvs.parser
import sys
import hgvs.dataproviders.uta
import hgvs.assemblymapper

hgvs_p = sys.argv[1]
hp = hgvs.parser.Parser()
var_p = hp.parse_hgvs_variant(hgvs_p)
print(var_p.format(conf={"p_3_letter": False}))
