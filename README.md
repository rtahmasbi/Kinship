# Kinship
Calculate Kinship matrix


# Download
    git clone https://github.com/rtahmasbi/Kinship

# Compile
    make


# Examples
## Example 1
In this example, all the individuals are in one file.

    kinship --file_pedigree pedigree1.txt --out test1

## Example 2
You can also calculate kinship from different generation files. For example, to calculate kinship form `GeneEvolve` output info files for the last 6 generations (generations 5 to 10), you can run the following commands:

    # step 1: create a list file for pedigree info:
    rm file_pedigree_list.txt; 
    for i in $(seq 5 10); do echo -e out1.info.pop1.gen$i.txt >> file_pedigree_list.txt; done
    # step 2: run kinship
    kinship --file_pedigree_list file_pedigree_list.txt
  
For more info on `GeneEvolve`, please check https://github.com/rtahmasbi/GeneEvolve.
  
