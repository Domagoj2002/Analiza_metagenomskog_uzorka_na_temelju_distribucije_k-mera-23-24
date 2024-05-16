#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <memory>
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <functional> // for std::reference_wrapper
#include "bioparser/fasta_parser.hpp" 
using namespace std;

 class Sequence  //prema  uputama za korištenje bioparser modula https://github.com/rvaser/bioparser
  { 
  public:
    Sequence(const char *name, std::uint32_t name_len, const char *data,
             std::uint32_t data_len)
    {
      name_.assign(name, name_len);
      data_.assign(data, data_len);
    }

  public:
    std::string name_, data_;
  };

string get_data_directory() { //za dohvaćanje puta do direktorija data
    filesystem::path current_path = filesystem::current_path();
    filesystem::path data_directory = current_path.parent_path().parent_path() / "data";
    return data_directory.string();
}

uint64_t sequence_length(string name) {  //fja za dohvaćanje duljine očitanja iz metapodataka fasta datoteke
    size_t last_slash = name.rfind('/');
    size_t last_underscore = name.rfind('_');
    //cout << "izracuanata distribucija1"<<endl;
    if (last_slash != string::npos && last_underscore != string::npos) {
        string start_str = name.substr(last_slash + 1, last_underscore - last_slash - 1);
        string end_str = name.substr(last_underscore + 1);
        //cout << "izracuanata distribucija2"<<endl;
        uint64_t start = 0;
        uint64_t end = 0;
        try {
            start = stoull(start_str);
            end = stoull(end_str);
        } catch (const std::invalid_argument& ia) {
            printf("Neispravan format u fasta datoteci: %s\n", ia.what());
            return 0;
        } catch (const std::out_of_range& oor) {
            printf("Predugačak broj u fasta datoteci: %s\n", oor.what());
            return 0;
        }
        return end - start;
    }
    printf("Neispravan format u fasta datoteci!");
    return 0;
}



vector<double> get_distribution_vector(const unique_ptr<Sequence>& fragment, int k, bool referent) {
    string name = fragment->name_;
    uint64_t length;
    if(referent==false){
        length = sequence_length(name);       
    }
    else{
        length = fragment->data_.size();
    }
    length = length - k + 1;
    //cout<<length<<endl;
    map<string, double> kmer_counts;
    string sequence = fragment->data_;
    vector<double> distribution_vector;

    // Initialization of the map
    for (int i = 0; i < pow(4, k); i++) {
        string kmer = "";
        int x = i;
        for (int j = 0; j < k; j++) {
            char base = "ACGT"[x % 4];
            kmer = base + kmer;
            x /= 4;
        }
        kmer_counts[kmer] = 0;
    }

    // Counting k-mers
    for (uint64_t i = 0; i < length; i++) {
        string kmer = sequence.substr(i, k);
        kmer_counts[kmer]++;
    }
    
    // Calculating distribution
    for (auto& pair : kmer_counts) {
        pair.second /= length;
    }

    // Adding components to the vector in the same order as the addition
    for (int i = 0; i < pow(4, k); i++) {
        string kmer = "";
        int x = i;
        for (int j = 0; j < k; j++) {
            char base = "ACGT"[x % 4];
            kmer = base + kmer;
            x /= 4;
        }
        distribution_vector.push_back(kmer_counts[kmer]);
    }

    return distribution_vector;
}



map<string, double> get_distribution_vector_map(const unique_ptr<Sequence>& fragment, int k) { //verzija fje sa mapama
    string name = fragment->name_;
    uint64_t length = sequence_length(name);
    length=length-k+1;
    map<string, double> kmer_counts;
    string sequence = fragment->data_;  

    // Inicijalizacija mape
    for (int i = 0; i < pow(4, k); i++) {
        string kmer = "";
        int x = i;
        for (int j = 0; j < k; j++) {
            char base = "ACGT"[x % 4];
            kmer = base + kmer;
            x /= 4;
        }
        kmer_counts[kmer] = 0;
    }

    // Brojanje k-mera
    for (uint64_t i = 0; i < length; i++) {
        string kmer = sequence.substr(i, k);
        kmer_counts[kmer]++;
    }

    // Izračunavanje distribucije
    for (auto& pair : kmer_counts) {
        pair.second /= length;
    }

    return kmer_counts;
}

void write_csv(const vector<reference_wrapper<vector<double>>>& references_to_distribution_vectors, const vector<string>& names, const string& filename) {
    if (references_to_distribution_vectors.size() != names.size()) {
        throw runtime_error("Veličina vektora referenci i vektora imena se ne podudara.");
    }
    
    ofstream file(filename);
    if (!file.is_open()) {
    throw runtime_error("Ne mogu otvoriti datoteku: " + filename);
}
    
    // Write the header
    for (size_t i = 0; i < names.size(); i++) {
        file << names[i];
        if (i != names.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
    
    // Write the data
    size_t length = references_to_distribution_vectors[0].get().size();
    for (size_t i = 0; i < length; i++) {
        for (size_t j = 0; j < references_to_distribution_vectors.size(); j++) {
            if (i < references_to_distribution_vectors[j].get().size()) {
                file << references_to_distribution_vectors[j].get()[i];
            }
            if (j != references_to_distribution_vectors.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }
    
    file.close();
}


double get_cosinus_similarity(std::map<std::string, double> vector1, std::map<std::string, double> vector2) {
    double dot_product = 0.0;
    double norm1 = 0.0;
    double norm2 = 0.0;

    for (const auto& pair : vector1) {
        const auto& key = pair.first;
        const auto& value1 = pair.second;

        double value2 = 0.0;
        if (vector2.count(key) > 0) {
            value2 = vector2[key];
        }

        dot_product += value1 * value2;
        norm1 += value1 * value1;
    }

    for (const auto& pair : vector2) {
        const auto& value2 = pair.second;
        norm2 += value2 * value2;
    }

    norm1 = std::sqrt(norm1);
    norm2 = std::sqrt(norm2);

    if (norm1 == 0.0 || norm2 == 0.0) {
        throw std::invalid_argument("Norma jednog od vektora je 0, kosinusna sličnost nije definirana.");
    }

    return dot_product / (norm1 * norm2);
}


int main(){
          // najprije ide file sa cjelokupnim genomom, a potom fragmenti
    // auto genome_parser =
    //       bioparser::Parser<Sequence>::Create<FastaParser>( //izbrisano bioparser::
    //           argv[optind]);
    // auto genomes = genome_parser->Parse(-1);
    //---------------------------------------------------
    // auto fragment_parser =
    //          bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
    //              "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/NCTC4450.fasta");
    // auto fragments = fragment_parser->Parse(-1); //fragments je vektor unique_ptr<Sequence> objekata
    // cout << fragments[0]->name_ << endl;
    // cout << fragments[0]->data_ << endl;
    // cout <<"prva 3 znaka ocitanja su"<< fragments[0]->data_.substr(0,3) << endl;
    // cout << "duljina jednog očitanja"<< fragments[0]->data_.size()<< endl;
    // cout << "Broj fragmenata: " << fragments.size() << endl;
    // cout << sequence_length(">m151004_144909_00127_c100873272550000001823191402121654_s1_p0/8/4883_9223") << endl;

    // map<string, double> mapa;
    // mapa=get_distribution_vector_map(ref(fragments[0]),3); //koristim referencu da izbjegnem kopiranje podataka

    //     // Ispis sadržaja mape
    // cout << "Sadržaj mape kmer_counts:" << endl;
    // for (const auto& pair : mapa) {
    //     cout << pair.first << ": " << pair.second << endl;
    // }

    // // for (auto& fragment : fragments) { //unique pointer se ne može kopirati pa koristiš referencu
    // //     printf(fragment->data_);
    // // }

    // printf("Sve funkcionira\n");
     //---------------------------------------------------

     //----parsiranje referentnih genoma------------------
    auto ref_genome_parser1 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/bacillus_cereus_reference.fasta");
    auto ref_genome_parser2 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/escherichia_coli_reference.fasta");
    auto ref_genome_parser3 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/haemophilus_influenzae_reference.fasta");
    auto ref_genome_parser4 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/helicobacter_pylori_reference.fasta");
    auto ref_genome_parser5 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/legionella_spiritensis_reference.fasta");
    auto ref_genome_parser6 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/proteus_vulgaris_reference.fasta");
    auto ref_genome_parser7 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/pseudomonas_aeruginosa_reference.fasta");
    auto ref_genome_parser8 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/salmonella_enterica_reference.fasta");
    auto ref_genome_parser9 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/streptococcus_pneumoniae_reference.fasta");
    auto ref_genome_parser10 =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/tissierella_praeacuta_reference.fasta");
    
    auto ref1 = ref_genome_parser1->Parse(-1); //ref1...10 su vektori unique_ptr<Sequence> objekata
    auto ref2 = ref_genome_parser2->Parse(-1); 
    auto ref3 = ref_genome_parser3->Parse(-1); 
    auto ref4 = ref_genome_parser4->Parse(-1); 
    auto ref5 = ref_genome_parser5->Parse(-1); 
    auto ref6 = ref_genome_parser6->Parse(-1); 
    auto ref7 = ref_genome_parser7->Parse(-1); 
    auto ref8 = ref_genome_parser8->Parse(-1); 
    auto ref9 = ref_genome_parser9->Parse(-1); 
    auto ref10 = ref_genome_parser10->Parse(-1); 

    //----kreiranje distribucijskih vektora referentnih genoma------------------

    vector<double> ref_distribution_vector1;
    vector<double> ref_distribution_vector2;
    vector<double> ref_distribution_vector3;
    vector<double> ref_distribution_vector4;
    vector<double> ref_distribution_vector5;
    vector<double> ref_distribution_vector6;
    vector<double> ref_distribution_vector7;
    vector<double> ref_distribution_vector8;
    vector<double> ref_distribution_vector9;
    vector<double> ref_distribution_vector10;

    int k = 3;

    ref_distribution_vector1 = get_distribution_vector(ref(ref1[0]),k,true); //koristim referencu da izbjegnem kopiranje podataka)
    ref_distribution_vector2 = get_distribution_vector(ref(ref2[0]),k,true);
    ref_distribution_vector3 = get_distribution_vector(ref(ref3[0]),k,true);
    ref_distribution_vector4 = get_distribution_vector(ref(ref4[0]),k,true);
    ref_distribution_vector5 = get_distribution_vector(ref(ref5[0]),k,true);
    ref_distribution_vector6 = get_distribution_vector(ref(ref6[0]),k,true);
    ref_distribution_vector7 = get_distribution_vector(ref(ref7[0]),k,true);
    ref_distribution_vector8 = get_distribution_vector(ref(ref8[0]),k,true);
    ref_distribution_vector9 = get_distribution_vector(ref(ref9[0]),k,true);
    ref_distribution_vector10 = get_distribution_vector(ref(ref10[0]),k,true);


    // Stvaranje vektora referenci
vector<reference_wrapper<vector<double>>> references_to_distribution_vectors;

// Dodavanje referenci na distribucijske vektore
references_to_distribution_vectors.push_back(ref_distribution_vector1);
references_to_distribution_vectors.push_back(ref_distribution_vector2);
references_to_distribution_vectors.push_back(ref_distribution_vector3);
references_to_distribution_vectors.push_back(ref_distribution_vector4);
references_to_distribution_vectors.push_back(ref_distribution_vector5);
references_to_distribution_vectors.push_back(ref_distribution_vector6);
references_to_distribution_vectors.push_back(ref_distribution_vector7);
references_to_distribution_vectors.push_back(ref_distribution_vector8);
references_to_distribution_vectors.push_back(ref_distribution_vector9);
references_to_distribution_vectors.push_back(ref_distribution_vector10);

// Nazivi bakterija
vector<string> names = {
        "b. cereus", "e. coli", "h. influenzae", "h. pylori", "l. spiritensis",
        "p. vulgaris", "p. aeruginosa", "s. enterica", "s. pneumoniae", "t. praeacuta"
    };


//ispis u csv datoteku
string filename = "reference_vectors.csv";
string absolute_path = get_data_directory() + "/" + filename;
write_csv(references_to_distribution_vectors, names, absolute_path);

//----parsiranje očitanja-----------------------------------------
auto fragment_parser =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/NCTC4450.fasta");
auto fragments = fragment_parser->Parse(-1); //fragments je vektor unique_ptr<Sequence> objekata

//----kreiranje distribucijskih vektora očitanja------------------


    return 0;
}