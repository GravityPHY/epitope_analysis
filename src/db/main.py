import ssl
from pymongo import MongoClient
from pymongo.server_api import ServerApi
# Fetching environment variables
#USERNAME = config('MONGODB_USERNAME')
#PASSWORD = config('MONGODB_PASSWORD')

uri = f"mongodb+srv://imhaoyu:987123YH@sequence.b5wpc.mongodb.net/?retryWrites=true&w=majority&appName=sequence"
# Create a new client and connect to the server
client = MongoClient(uri, ssl=True, ssl_cert_reqs=ssl.CERT_NONE)
# Send a ping to confirm a successful connection
try:
    client.admin.command('ping')
    print("Pinged your deployment. You successfully connected to MongoDB!")
except Exception as e:
    print(e)
databases = client.list_database_names()
print(databases)