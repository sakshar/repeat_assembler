from Bio import SeqIO

hifi = SeqIO.parse('reads/old/hifi.chrXII_complex.fasta', 'fasta')
nano = SeqIO.parse('reads/old/nano.chrXII_complex.fasta', 'fasta')

hifi_map = dict()
nano_map = dict()
for record in hifi:
    hifi_map[record.id] = len(record.seq)
for record in nano:
    nano_map[record.id] = len(record.seq)

print(len(hifi_map))
print(len(nano_map))

#hifi_sorted = dict(sorted(hifi_map.items(), key=lambda item: item[1], reverse=True))
#nano_sorted = dict(sorted(nano_map.items(), key=lambda item: item[1], reverse=True))
ranges = ['0-10k', '10k-20k', '20k+']
hifi_range = {'0-10k': 0, '10k-20k': 0, '20k+': 0}
nano_range = {'0-10k': 0, '10k-20k': 0, '20k+': 0}

for key in hifi_map.keys():
    if hifi_map[key] < 10000:
        hifi_range[ranges[0]] += 1
    elif hifi_map[key] < 20000:
        hifi_range[ranges[1]] += 1
    else:
        hifi_range[ranges[2]] += 1

for key in nano_map.keys():
    if nano_map[key] < 10000:
        nano_range[ranges[0]] += 1
    elif nano_map[key] < 20000:
        nano_range[ranges[1]] += 1
    else:
        nano_range[ranges[2]] += 1

print('HiFi reads\' range:', hifi_range)
print('ONT reads\' range:', nano_range)
