# skripta koja izvlači sve gene iz mreže u .sif formatu 
def parse_line(data_line):
    parts = data_line.split('-')
    identifier1 = parts[0].strip()
    identifier2 = parts[1].strip()
    return identifier1, identifier2

# postaviti iz koje mreže se izvlače geni
if __name__ == "__main__":
    s = set()
    with open("") as f:
        for line in f:
            id1, id2 = parse_line(line)
            s.add(id1)
            s.add(id2)
    print(s)