import sys

f = open(sys.argv[1], 'r')
f1 = open("filtered_"+sys.argv[1], 'w')
count = 0
while True:
    line = f.readline()
    if not line:
        break
    line1 = f.readline()
    #if line[0] != '>':
    count += 1
    line1 = line1.strip()
    i = 0
    faulty = 0
    while i < len(line1):
        if line1[i] not in ['A', 'C', 'G', 'T']:
            faulty = 1
            print(count, "Faulty read found")
            break
        i += 1
    if faulty == 0:
        f1.write(line)
        f1.write(line1+'\n')
