import json
import os

configfile = os.path.join(os.path.dirname(__file__), 'config.json')
print(configfile)
if hasattr(__builtins__, 'raw_input'):
    input = raw_input


with open(configfile, 'r') as f:
    myconfig = json.load(f)

if myconfig['databasepath'] is None:
    myconfig['databasepath'] = input(("Specify the database"
    " path if LabDB is installed in the system: "))
    with open(configfile, 'wb') as f:
        json.dump(myconfig, f)
