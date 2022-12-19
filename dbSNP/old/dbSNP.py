import time
from Bio import Entrez

Entrez.email = "txdshsnbsdna@gmail.com"
Entrez.api_key = "699902df92858d290ac7526781abb3d32c08"

# dbSNP supported query terms (https://www.ncbi.nlm.nih.gov/snp/docs/entrez_help/) can be build and test online using web query builder (https://www.ncbi.nlm.nih.gov/snp/advanced)
# esearch handle
handle = Entrez.esearch(db="snp",  # search dbSNP
                          term="rs58991260",
                          usehistory="y", #cache result on server for download in batches
                         )
response = handle.read()

print(response)


