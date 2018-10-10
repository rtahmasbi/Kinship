
#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <ostream>
#include <iostream>
#include <ctime>
#include <vector>
#include <string> // str1.compare(str2)
#include <sstream> //std::istringstream
#include <cmath> //pow
#include <algorithm> //find

#define VER "1.0"



class Human
{
public:
    std::string ID;
    std::string ID_Father;
    std::string ID_Mother;
    long int idx_fa;
    long int idx_mo;

    Human() {};

    Human(std::string id, std::string id_fa, std::string id_mo) {
        ID = id;
        ID_Father = id_fa;
        ID_Mother = id_mo;
        idx_fa=-1;
        idx_mo=-1;
    };

    friend std::ostream &operator << (std::ostream &output, const Human &h)
    {
        output << "ID=" << h.ID << ", ID_Father=" << h.ID_Father << ", ID_Mother=" << h.ID_Mother;
        return output;
    }

    friend bool operator== (const Human &h1, const Human &h2)
    {
        return h1.ID == h2.ID;
    }

};





long int ras_read_fam_file(std::string file_fam, std::vector<Human> &humans, long int gen);
bool ras_read_pedigree_file(std::string file_pedigree, std::vector<Human> &humans);
void ras_help(void);
int ras_read_in_gen_file(std::string in_gen_info, std::vector<std::string> &vec_gen_file_name);
double ras_similarity(const std::vector<Human> &humans, long int i, long int j);
bool write_similarities(std::string file_out, const std::vector<Human> &humans, long int pedigree_st, long int pedigree_en);
bool ras_set_index(std::vector<Human> &humans);



int main(int argc, char ** argv)
{
    int start_time = time(0);

    std::cout << "Rasool Tahmasbi; ver " << VER << std::endl;

    //proc params
    std::string file_pedigree_list="";
    std::string file_pedigree="";
    std::string file_out="out";
    int generation_info_type=1; // 1=pedigree, 2=fam files

    for (int i=1; i<argc; i++)
    {
        std::string s(argv[i]);
        if (s=="--help")
        {
            ras_help();
            return 0;
        }
        else if (s=="--file_pedigree")
        {
            file_pedigree=std::string(argv[++i]);
            generation_info_type=1;
        }
        else if (s=="--file_pedigree_list")
        {
            file_pedigree_list=std::string(argv[++i]);
            generation_info_type=2;
        }
        else if (s=="--out")
        {
            file_out=std::string(argv[++i]);
        }
        else
        {
            std::cout << "Error: unknown parameter [" << s << "]" << std::endl;
        }
    }


    /////////////////////////////////////////////////////////
    // read human info from pedigree or fam files
    std::vector<Human> humans;
    humans.clear();
    long int pedigree_st=0, pedigree_en=0;


    if(generation_info_type==1)
    {
        if(!ras_read_pedigree_file(file_pedigree, humans))
            return -1;
    }

    if(generation_info_type==2)
    {
        std::vector<std::string> vec_gen_file_name;
        int ngen=0;
        ngen = ras_read_in_gen_file(file_pedigree_list, vec_gen_file_name);
        std::cout << " Number of generations = " << ngen << std::endl;
        for (int igen=0; igen<ngen; igen++)
        {
            if(igen==ngen-1) pedigree_st=humans.size(); // compute just for the last generation
            if (!ras_read_fam_file(vec_gen_file_name[igen], humans, igen+1))
                return -1;
        }
    }

    std::cout << " Size of pedigree for all the specified generations = " << humans.size() << std::endl;
    pedigree_en=humans.size();

    /////////////////////////////////////////////////////////
    std::cout << " Start finding indexes" << std::endl;
    ras_set_index(humans);


    /////////////////////////////////////////////////////////
    // compute kinship
    std::cout << " Start computing kinship" << std::endl;
    write_similarities(file_out, humans, pedigree_st, pedigree_en);


    /////////////////////////////////////////////////////////
    // end
    int time_tot = time(0) - start_time;
    std::cout << std::endl << " Program Successfully finished." << std::endl;
    std::cout << "Total Run completed in " << time_tot / 3600 << " hours, " << (time_tot % 3600) / 60 << " mins, " << time_tot % 60 << " seconds." << std::endl;

    return 0;
}


bool write_similarities(std::string file_out, const std::vector<Human> &humans, long int pedigree_st, long int pedigree_en)
{
    std::string sep=" ";
    std::ofstream myfile((file_out+".kin").c_str());
    long int tot=(pedigree_en-pedigree_st)*(pedigree_en-pedigree_st+1)/2;
    long int r=0;

    if (myfile.is_open())
    {
        for (long int i=pedigree_st; i<pedigree_en; i++)
        {
            for (long int j=pedigree_st; j<=i; j++)
            {
                myfile << humans[i].ID << sep << humans[j].ID << sep << ras_similarity(humans,i,j) << std::endl;
                if (r%10000==0) std::cout << "\r " << r << " of " << tot << " wrote ..." << std::flush;
                r++;
            }
        }
        std::cout << "\r " << "All wrote successfully." << std::flush << std::endl;
        myfile.close();
    }
}


double ras_similarity(const std::vector<Human> &humans, int i, int j)
{
    if (i==-1 || j==-1)
        return 0;

    if (i==j)
    {
        long int idx_hi_f=humans[i].idx_fa;
        long int idx_hi_m=humans[i].idx_mo;

        return 1+ras_similarity(humans,idx_hi_f,idx_hi_m);
    }
    else
    {
        long int idx_hi_f=humans[i].idx_fa;
        long int idx_hi_m=humans[i].idx_mo;
        long int idx_hj_f=humans[j].idx_fa;
        long int idx_hj_m=humans[j].idx_mo;

        double d_i_jf=0, d_i_jm=0, d_j_if=0, d_j_im=0;
        if (idx_hj_f>-1)
            d_i_jf=ras_similarity(humans,i,idx_hj_f);
        if (idx_hj_m>-1)
            d_i_jm=ras_similarity(humans,i,idx_hj_m);
        if (idx_hi_f>-1)
            d_j_if=ras_similarity(humans,j,idx_hi_f);
        if (idx_hi_m>-1)
            d_j_im=ras_similarity(humans,j,idx_hi_m);

        return std::max((d_i_jf+d_i_jm)/2,(d_j_if+d_j_im)/2);

    }
}


bool ras_set_index(std::vector<Human> &humans)
{
    std::string id;
    long int pos;
    for (long int i=0; i<humans.size(); i++)
    {
        if (i%1000==0)
            std::cout << "\r   " << "running index " << i << " of " << humans.size() << " ..." << std::flush;
        // ID_Father
        id=humans[i].ID_Father;
        pos = std::find(humans.begin(), humans.end(), Human(id,"NA","NA")) - humans.begin();
        if(pos >= humans.size()) pos=-1;
        humans[i].idx_fa=pos;
        // ID_Mother
        id=humans[i].ID_Mother;
        pos = std::find(humans.begin(), humans.end(), Human(id,"NA","NA")) - humans.begin();
        if(pos >= humans.size()) pos=-1;
        humans[i].idx_mo=pos;
    }
    std::cout << "\r " << "ras_set_index done." << std::flush << std::endl;
}





// read file_pedigree_list
int ras_read_in_gen_file(std::string file_pedigree_list, std::vector<std::string> &vec_gen_file_name)
{
    std::ifstream ifile(file_pedigree_list.c_str());
    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+ file_pedigree_list +"] to read." << std::endl;
        return 0;
    }

    std::string line;

    while (std::getline(ifile, line))
    {
        vec_gen_file_name.push_back(line);
    }
    return vec_gen_file_name.size();
}





int ras_read_fam_file(std::string file_fam, std::vector<Human> &humans, int gen)
{
    char sep=' ';


    std::string file_name=file_fam;
    std::ifstream ifile(file_name.c_str());

    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
        return 0;
    }

    std::string line;

    // discard the first line which is header
    std::getline(ifile, line);

    while (std::getline(ifile, line))
    {
        std::istringstream iss(line);
        std::string token;

        std::getline(iss, token, sep); // ID
        std::string ID=std::to_string(gen) + "$" + token;
        std::getline(iss, token, sep); // ID_Father
        std::string ID_Father=std::to_string(gen-1) + "$" + token;
        std::getline(iss, token, sep); // ID_Mother
        std::string ID_Mother=std::to_string(gen-1) + "$" + token;

        Human h;
        h.ID=ID;
        h.ID_Father=ID_Father;
        h.ID_Mother=ID_Mother;
        humans.push_back(h);
    }
}




bool ras_read_pedigree_file(std::string file_pedigree, std::vector<Human> &humans)
{
    char sep=' ';


    std::string file_name=file_pedigree;
    std::ifstream ifile(file_name.c_str());

    if(!ifile)
    {
        std::cout << "Error: can not open the file ["+ file_name +"] to read." << std::endl;
        return false;
    }

    std::string line;

    // discard the first line which is header
    std::getline(ifile, line);

    while (std::getline(ifile, line))
    {
        std::istringstream iss(line);
        std::string token;

        std::getline(iss, token, sep); // ID
        std::string ID=token;
        std::getline(iss, token, sep); // ID_Father
        std::string ID_Father=token;
        std::getline(iss, token, sep); // ID_Mother
        std::string ID_Mother=token;

        Human h;
        h.ID=ID;
        h.ID_Father=ID_Father;
        h.ID_Mother=ID_Mother;
        humans.push_back(h);
    }
    return true;
}



void ras_help(void)
{
    std::cout << std::endl;
    std::cout << "Help" << std::endl;
    std::cout << " kinship --file_pedigree pedigree1.txt" << std::endl;
    std::cout << "         [pedigree1.txt] is a file with 3 columns: ID, ID_Father, ID_Mother" << std::endl;
    std::cout << std::endl;
    std::cout << " kinship --file_pedigree_list file_pedigree_list.txt" << std::endl;
    std::cout << "         [file_pedigree_list.txt] is the list of pedigree files, where each file's name is in one line: gen0, gen1, ..." << std::endl;
    std::cout << "         each file has (3+) columns: ID, ID_Father, ID_Mother" << std::endl;

}
