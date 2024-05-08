#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <memory>
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


uint64_t sequence_length(string name) {  //fja za dohvaćanje duljine očitanja iz metapodataka fasta datoteke
    size_t last_slash = name.rfind('/');
    size_t last_underscore = name.rfind('_');
    if (last_slash != string::npos && last_underscore != string::npos) {
        string start_str = name.substr(last_slash + 1, last_underscore - last_slash - 1);
        string end_str = name.substr(last_underscore + 1);
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

map<string, double> get_distribution_vector(const unique_ptr<Sequence>& fragment, int k) {
    string name = fragment->name_;
    uint64_t length = sequence_length(name);
    length=length-k+1;
    map<string, double> kmer_counts;
    string sequence = fragment->data_;  // Pretpostavljam da postoji metoda za dohvaćanje sekvence

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



int main(){
          // najprije ide file sa cjelokupnim genomom, a potom fragmenti
    // auto genome_parser =
    //       bioparser::Parser<Sequence>::Create<FastaParser>( //izbrisano bioparser::
    //           argv[optind]);
    // auto genomes = genome_parser->Parse(-1);
    auto fragment_parser =
             bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                 "/home/domagoj/Desktop/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24/data/NCTC4450.fasta");
    auto fragments = fragment_parser->Parse(-1); //fragments je vektor unique_ptr<Sequence> objekata
    cout << fragments[0]->name_ << endl;
    cout << fragments[0]->data_ << endl;
    cout <<"prva 3 znaka ocitanja su"<< fragments[0]->data_.substr(0,3) << endl;
    cout << "duljina jednog očitanja"<< fragments[0]->data_.size()<< endl;
    cout << "Broj fragmenata: " << fragments.size() << endl;
    cout << sequence_length(">m151004_144909_00127_c100873272550000001823191402121654_s1_p0/8/4883_9223") << endl;

    map<string, double> mapa;
    mapa=get_distribution_vector(ref(fragments[0]),3); //koristim referencu da izbjegnem kopiranje podataka

        // Ispis sadržaja mape
    cout << "Sadržaj mape kmer_counts:" << endl;
    for (const auto& pair : mapa) {
        cout << pair.first << ": " << pair.second << endl;
    }

    // for (auto& fragment : fragments) { //unique pointer se ne može kopirati pa koristiš referencu
    //     printf(fragment->data_);
    // }

    printf("Sve funckonira\n");
    return 0;
}