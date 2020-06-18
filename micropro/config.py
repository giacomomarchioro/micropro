import json
if hasattr(__builtins__, 'raw_input'):
    input = raw_input

with open('config.json', 'r') as f:
    mydict = json.load(f)
