import pandas as pd

orgs = pd.read_table('speclist.txt', sep='=', skiprows=59, header=None, nrows=58800)
locs = pd.read_table('subcell.txt', skiprows=43, header=None)
mhcs = pd.read_table('allelenames.txt', sep=" ", header=None)


# 3rd column of orgs is the org name
org_names = orgs[1].tolist()

# keep only rows with ID or AC in the first column from locs
locs_id_ac = locs[locs[0].str.contains('ID|AC')]
# remove dots and first 7 character from elements in locs_id_ac
locs_id_ac = locs_id_ac[0].str.replace(r'.', '').str[5:].tolist()

# keep 2nd column as list
mhcs = mhcs[1].tolist()


# write lists each to a file
with open('orgs.txt', 'w') as f:
    for name in org_names:
        f.write(name + '\n')
with open('locs.txt', 'w') as f:
    for loc in locs_id_ac:
        f.write(loc + '\n')
with open('mhcs.txt', 'w') as f:
    for mhc in mhcs:
        f.write(mhc + '\n')


