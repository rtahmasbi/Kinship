# Kinship
Calculate Kinship matrix


# Download
  git clone https://github.com/rtahmasbi/Kinship

# Compile
  make


# Examples
## Example 1
  kinship --file_pedigree pedigree1.txt --out test1

## Example 2
  rm file_pedigree_list.txt; 
  for i in $(seq 5 10); do echo -e out1.info.pop1.gen$i.txt >> file_pedigree_list.txt; done
  kinship --file_pedigree_list file_pedigree_list.txt
  
  
  
