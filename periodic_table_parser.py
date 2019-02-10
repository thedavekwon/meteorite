
def main():
    pt = open("periodic_table.txt", mode='r')
    element_table = {}

    for line in pt:
        words = line.split(',')
        element_table[words[1][1:]] = words[3]



main()