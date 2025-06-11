import ssl
from pymongo import MongoClient
from pymongo.server_api import ServerApi

uri = f"mongodb+srv://imhaoyu:987123YH@sequence.b5wpc.mongodb.net/?retryWrites=true&w=majority&appName=sequence"
# Create a new client and connect to the server
client = MongoClient(uri, ssl=True, ssl_cert_reqs=ssl.CERT_NONE)
# Send a ping to confirm a successful connection
try:
    client.admin.command('ping')
    print("Pinged your deployment. You successfully connected to MongoDB!")
except Exception as e:
    print(e)
db=client["sequence_db"]
collection=db["sequences"]

fasta_path="/projectnb2/docking/imhaoyu/25_epitope_pLM/amino_classification/LoRAforESM-2/epitope_dataset/discotope3.0_downloads/trainval_solved/trainval_solved.fasta"

with open(fasta_path,"r") as fasta_file:
    records ={}
    current_id = None
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):
            current_id = line[1:]
            records[current_id] = ""
        else:
            records[current_id] += line

documents = [{"pdb_id":k, "sequence":v} for k,v in records.items()]
collection.insert_many(documents)

record = collection.find_one({"pdb_id":"1a2y_C"})
print(record)