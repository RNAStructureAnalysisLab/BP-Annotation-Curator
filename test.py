import re

RESIDUE_TEMPLATE = r'([a-zA-Z]+|\d[a-zA-Z]+|\d)(-?\d+)'

string = "08"

groups = re.match(RESIDUE_TEMPLATE, string)

print(groups.group(1))
print(groups.group(2))
