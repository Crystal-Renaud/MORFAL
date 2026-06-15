import os
import argparse

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--MORFEE', dest='MORFEE', type=str,
                        help='/path/to/MORFEE_DB')
    parser.add_argument('-o', '--MORFEE_ind', dest='MORFEE_ind', type=str,
                        help='/path/to/MORFEE_ind')
    return parser.parse_args()


# create index file
def indexing(MORFEE, MORFEE_ind):
    with open(f'{MORFEE}', 'r') as mor_file:
        with open(f'{MORFEE_ind}', 'w') as file:
            for line_morf in mor_file:
                if line_morf.startswith('#'):
                    en_tete = line_morf.split("\t")
                    file.write(f"{en_tete[0]}\t{en_tete[1]}\n")
                else : 
                    line_morfee = line_morf.split("\t")
                    chrom_morf = line_morfee[0]
                    start = int(line_morfee[1])
                    file.write(f"{chrom_morf}\t{start}\n")


def main():

    # get arguments
    args = get_arguments()
    MORFEE = os.path.realpath(args.MORFEE)
    MORFEE_ind = os.path.realpath(args.MORFEE_ind)


    indexing(MORFEE=MORFEE, MORFEE_ind=MORFEE_ind)

if __name__ == "__main__":
    main()