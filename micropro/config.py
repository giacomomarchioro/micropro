import json
import os

configfile = os.path.join(os.path.dirname(__file__), 'config.json')
print(configfile)
if hasattr(__builtins__, 'raw_input'):
    input = raw_input

if os.path.exists(configfile):
    with open(configfile, 'r') as f:
        myconfig = json.load(f)

else:
    with open(configfile, 'wb') as f:
        myconfig = {}
        myconfig['databasepath'] = input(("Specify the database"
    " path if LabDB is installed in the system: "))
        json.dump(myconfig, f)
