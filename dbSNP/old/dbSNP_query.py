import myvariant
import json

mv = myvariant.MyVariantInfo()
x=json.dumps(mv.query('rs58991260'), indent=1)

print(x)
