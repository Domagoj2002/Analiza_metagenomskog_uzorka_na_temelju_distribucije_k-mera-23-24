#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <memory>
#include <fstream>
#include <stdexcept>
#include <regex>
#include <filesystem>
#include <random>
#include <functional> // for std::reference_wrapper
#include "bioparser/fasta_parser.hpp" 


using namespace std;
namespace fs = std::filesystem;

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
        length = (int) length;       
    }
    else{
        length = fragment->data_.size();
    }
    length = length - k + 1;
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
        printf("Velicina vektora referenci je %d,a velicina vektora imena je %d\n",(int)references_to_distribution_vectors.size(),(int)names.size());
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
    size_t length = references_to_distribution_vectors[0].get().size();//to je broj redaka

    for (size_t i = 0; i < length; i++) { //redak
        for (size_t j = 0; j < references_to_distribution_vectors.size(); j++) { //stupac           
                file << references_to_distribution_vectors[j].get()[i];         
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

// Funkcija za izračunavanje srednje vrijednosti
double calculateMean(const vector<int>& data, int start, int count) {
    double sum = 0.0;
    for (int i = start; i < start + count; ++i) {
        sum += data[i];
    }
    return sum / count;
}

// Funkcija za izračunavanje standardne devijacije
double calculateStdDev(const vector<int>& data, int start, int count, double mean) {
    double sum = 0.0;
    for (int i = start; i < start + count; ++i) {
        sum += pow(data[i] - mean, 2);
    }
    return sqrt(sum / count);
}

// Funkcija za pretvaranje fastq datoteke u fasta pomoću awk
void convert_fastq_to_fasta(const string& fastq_file, const string& data_directory) {
    string fasta_file = regex_replace(fastq_file, regex("\\.fastq$"), ".fasta");
    string command = "awk '{if(NR%4==1) {printf(\">%s\\n\", substr($0, 2));} else if(NR%4==2) print;}' " + data_directory +"/Reads/"+ fastq_file + " > " + data_directory +"/Reads/"+ fasta_file;
    system(command.c_str());
}


// Funkcija za parsiranje imena bakterija iz naziva datoteka
std::string parse_bacteria_name(const std::string& filename) {
    // Provjeravamo sadrži li naziv datoteke očekivani format
    size_t fasta_pos = filename.find(".fasta");
    if (fasta_pos != std::string::npos) {
        // Pronalazimo prvu poziciju podvlake
        size_t first_underscore_pos = filename.find('_');
        if (first_underscore_pos != std::string::npos) {
            // Provjeravamo postoji li "reference" u nazivu
            size_t ref_pos = filename.find("_reference");
            size_t second_underscore_pos = (ref_pos != std::string::npos) ? ref_pos : fasta_pos;
            
            // Formiramo rezultat s prvim slovom prve riječi i cijelom drugom riječi
            std::string result;
            result += std::tolower(filename[0]);
            result += ' ';
            result += filename.substr(first_underscore_pos + 1, second_underscore_pos - first_underscore_pos - 1);
            return result;
        }
    }
    // Ako format nije odgovarajući, vraćamo neizmijenjeni naziv datoteke
    return filename;
}


// Funkcija za kreiranje metagenomskog uzorka
void create_metagenomic_sample(const vector<string>& read_files, 
                               const vector<int>& read_counts, 
                               const string& output_file,const string& data_directory) {
    cout << output_file <<endl;
    ofstream outfile(output_file);
    if (!outfile.is_open()) {
    throw runtime_error("Ne mogu otvoriti datoteku: " + output_file);
}
    for (size_t i = 0; i < read_files.size(); ++i) {      
         cout << read_files[i]<< endl;
    }
    for (size_t i = 0; i < read_files.size(); ++i) { 
 
        string input_file = read_files[i];
        int count = read_counts[i];
        if (input_file.find(".fastq") != string::npos) {
            convert_fastq_to_fasta(input_file, data_directory);
            input_file = regex_replace(input_file, regex("\\.fastq$"), ".fasta");
        }
        string command = "awk 'BEGIN {RS=\">\"; ORS=\"\"} NR>1 {print \">\"$0}' " + data_directory + "/Reads/" + input_file + " | head -n " + to_string(count * 2) + " >> " + output_file;
        system(command.c_str());
    }
    outfile.close();
}

// Funkcija za parsiranje referentnih genoma
vector<vector<unique_ptr<Sequence>>> parse_reference_genomes(const vector<string>& reference_genomes, int k_mer_length, const vector<string>& names) {
    vector<vector<unique_ptr<Sequence>>> refs;
    for (const auto& genome_file : reference_genomes) {
        auto ref_genome_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(genome_file);
        auto ref = ref_genome_parser->Parse(-1);
        refs.push_back(std::move(ref));
    }

    vector<vector<double>> ref_distribution_vectors(refs.size());
    vector<reference_wrapper<vector<double>>> references_to_distribution_vectors;

    for (size_t i = 0; i < refs.size(); ++i) {
        ref_distribution_vectors[i] = get_distribution_vector(ref(refs[i][0]), k_mer_length, true);
        references_to_distribution_vectors.push_back(ref_distribution_vectors[i]);
    }

        //ispis u csv datoteku
    string filename = "reference_vectors.csv";
    string absolute_path = get_data_directory() + "/" + filename;
    write_csv(references_to_distribution_vectors, names, absolute_path); 

    return refs;
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
    // Definiranje putanje do CONFIGURATION_file.txt
    filesystem::path current_path = filesystem::current_path();
    filesystem::path path = current_path.parent_path().parent_path();
    string config_file = path.string() + "/CONFIGURATION_file.txt";
    string data_directory = get_data_directory();

    
    // Definiranje potrebnih varijabli
    vector<string> reference_genomes;
    vector<string> read_files;
    vector<int> read_counts;
    string line;
    int k;

    // Učitavanje konfiguracijske datoteke
    cout << "config file:" << config_file<< endl;
bool read_references = false, read_reads = false, k_mer = false;
    ifstream file(config_file);
    if (file.is_open()) {
        while (getline(file, line)) {
            // Preskoči prazne redove i one koji počinju s '#'
            if (line.empty() || line[0] == '#') {
                // Provjeri je li linija komentar za referentne genome ili očitanja
                if (line.find("NAZIVE DATOTEKA KOJE SADRŽE REFERENTNE GENOME") != std::string::npos) {
                    read_references = true;
                    read_reads = false;
                    k_mer = false;
                } else if (line.find("NAZIVE DATOTEKA KOJE SADRŽE OČITANJA") != std::string::npos) {
                    read_reads = true;
                    read_references = false;
                    k_mer = false;
                } else if (line.find("UNESITE ŽELJENU DULJINU K-MERA")!= std::string::npos){
                    read_reads = false;
                    read_references = false;
                    k_mer = true;
                }
                continue;
            }

            if (read_references) {
                // Učitavanje referentnih genoma
                reference_genomes.push_back(line);
            } else if (read_reads) {
                // Učitavanje očitanja i broja očitanja
                std::istringstream iss(line);
                std::string file_name;
                int count;
                iss >> file_name >> count;
                read_files.push_back(file_name);
                read_counts.push_back(count);
            } else if (k_mer){
                // Učitavanje željene duljine k-mera
                std::istringstream iss(line);
                iss >> k;
            }
        }
        file.close();
    } else {
        std::cerr << "Ne mogu otvoriti konfiguracijsku datoteku." << std::endl;
    }

    // Ispis učitanih podataka za provjeru
    std::cout << "Referentni genomi:" << std::endl;
    for (const auto& genome : reference_genomes) {
        std::cout << genome << std::endl;
    }
    std::cout << "Očitanja i broj očitanja:" << std::endl;
    for (size_t i = 0; i < read_files.size(); ++i) {
        std::cout << read_files[i] << " " << read_counts[i] << std::endl;
    }
    std::cout << "Željena duljina k-mera: " << k << std::endl;

    // Parsiranje imena bakterija za referentne genome
    vector<string> names;
    for (const auto& genome_file : reference_genomes) {
        names.push_back(parse_bacteria_name(genome_file));
    }

    // Kreiranje metagenomskog uzorka ako datoteka ne postoji
    string metagenomic_sample_file = data_directory + "/" + "metagenomic_sample.fasta";
    if (!fs::exists(metagenomic_sample_file)) {
        for (size_t i = 0; i < read_files.size(); ++i) {      
         cout << read_files[i]<< endl;
          cout << read_files.size()<< endl;
    }
        create_metagenomic_sample(read_files, read_counts, metagenomic_sample_file, data_directory);
    } else {
        cout << "Datoteka metagenomic_sample.fasta već postoji. Preskačem stvaranje nove datoteke." << endl;
    }

    // Parsiranje imena bakterija za metagenomski uzorak i broja očitanja
    vector<string> names_metagenomic_sample;
    for (const auto& read_file : read_files) {
        names_metagenomic_sample.push_back(parse_bacteria_name(read_file));
    }

    // Poklapanje redoslijeda imena bakterija i broja očitanja
    //vector<int> frags_num(read_counts);

// Parsiranje referentnih genoma i stvaranje distribucijskih vektora
    vector<vector<unique_ptr<Sequence>>> refs;
    for (const auto& genome_file : reference_genomes) {
        auto ref_genome_parser = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(data_directory + "/Referent_genomes/"+ genome_file);
        auto ref = ref_genome_parser->Parse(-1);
        refs.push_back(move(ref));
    }

    vector<vector<double>> ref_distribution_vectors(refs.size());
    vector<reference_wrapper<vector<double>>> references_to_distribution_vectors;

    for (size_t i = 0; i < refs.size(); ++i) {
        ref_distribution_vectors[i] = get_distribution_vector(ref(refs[i][0]), k, true);
        references_to_distribution_vectors.push_back(ref_distribution_vectors[i]);
    }

    // Ispis distribucijskih vektora u CSV datoteku
    string filename = "reference_vectors.csv";
    string absolute_path = data_directory + "/" + filename;
    write_csv(references_to_distribution_vectors, names, absolute_path);



    //----parsiranje očitanja-----------------------------------------
    auto fragment_parser =
                bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(
                    data_directory + "/" + "metagenomic_sample.fasta");
    auto fragments = fragment_parser->Parse(-1); //fragments je vektor unique_ptr<Sequence> objekata

    //----kreiranje distribucijskih vektora očitanja------------------

    // Stvaranje vektora referenci
    vector<reference_wrapper<vector<double>>> references_to_fragments_distribution_vectors;
    vector<shared_ptr<vector<double>>> real_vectors; 
    vector<int> frag_sizes; //pohranjuje duljinu svakog očitanja u uzorku

    for(int i=0;i<fragments.size();i++){
        if(fragments[i]->data_.size()<k){ //preskačem prekratka očitanja (manja od zadanog k)
            printf("Na poziciji %d u metagenomskom uzorku je prekratko ocitanje\n",i*2+2);
            int zbroj=0;
            for(int k=0;k<read_counts.size();k++){
                zbroj += read_counts[k]; 
                if(i<=zbroj-1){
                    read_counts[k]=read_counts[k]-1; //zbog izbacivanja smanjim broj
                    break;
                }
            }
        }
        else{
        frag_sizes.push_back(fragments[i]->data_.size()); //zapišem duljinu svakog očitanja
        // Stvaranje shared_ptr novog vektora za svaki fragment
        shared_ptr<vector<double>> fragment_distribution_vector = make_shared<vector<double>>(get_distribution_vector(fragments[i], k, false));

        //dodajem u vanjski vektor da bude vidljivo i izvan else bloka
        real_vectors.push_back(fragment_distribution_vector);
        // Pohrana referenci na nove vektore
        references_to_fragments_distribution_vectors.push_back(ref(*fragment_distribution_vector));
        }
    }

    //header csv datoteke
/*     vector<string> names_frags;
    for(int i = 1; i <= references_to_fragments_distribution_vectors.size(); i++) {
        names_frags.push_back("readr" + to_string(i));
    } */
    //header csv datoteke
    vector<string> names_frags;
    // Popunjavanje vektora names_frags
    for (size_t i = 0; i < read_counts.size(); ++i) {
        for (int j = 0; j < read_counts[i]; ++j) {
            names_frags.push_back(names_metagenomic_sample[i]);
        }
    }

    //ispis u csv datoteku
    string filename2 = "fragments_vectors.csv";
    string absolute_path2 = get_data_directory() + "/" + filename2;
    write_csv(references_to_fragments_distribution_vectors, names_frags, absolute_path2);



    //----kreiranje očitanja bez greške---------------------------

    //vector<int> frags_num = { 100, 91, 82, 106, 119, 97, 78, 111, 144, 70 }; //83 umjesto 82 za k<9 treba napisati


    vector<reference_wrapper<vector<double>>> references_to_NO_ERROR_fragments_distribution_vectors;
    vector<shared_ptr<vector<double>>> real_vectors2;

    int size_index = 0;
    random_device rd;
    mt19937 gen(rd());

    for (size_t ref_idx = 0; ref_idx < refs.size(); ++ref_idx) {
        for (int i = 0; i < read_counts[ref_idx]; ++i) {
            int frag_size = frag_sizes[size_index++];
            int max_start_idx = refs[ref_idx][0]->data_.size() - frag_size;
            uniform_int_distribution<> dis(0, max_start_idx); //iz jednolike razdiobe odabirem brojeve
            int start_idx = dis(gen);

            string fragment = refs[ref_idx][0]->data_.substr(start_idx, frag_size);
            //cout << "Generated fragment: " << fragment << endl;

            // Simulate creating a Sequence object from the fragment
            auto seq = make_shared<Sequence>("frag", 4, fragment.c_str(), fragment.size());

            // Convert shared_ptr to unique_ptr
            auto seq_unique = make_unique<Sequence>(*seq);

            // Stvaranje shared_ptr novog vektora za svaki fragment
            shared_ptr<vector<double>> fragment_distribution_vector = make_shared<vector<double>>(get_distribution_vector(seq_unique, k, true));

            //dodajem u vanjski vektor da bude vidljivo i izvan else bloka
            real_vectors2.push_back(fragment_distribution_vector);
            // Pohrana referenci na nove vektore
            references_to_NO_ERROR_fragments_distribution_vectors.push_back(ref(*fragment_distribution_vector));
        }
    }

    //ispis u csv datoteku
    string filename3 = "NO_ERROR_fragments_vectors.csv";
    string absolute_path3 = get_data_directory() + "/" + filename3;
    write_csv(references_to_NO_ERROR_fragments_distribution_vectors, names_frags, absolute_path3);



    //----ispis osnovnih statistika---------------------------



    // Pretpostavljamo da frag_sizes ima točan broj elemenata
    if (frag_sizes.size() != accumulate(read_counts.begin(), read_counts.end(), 0)) {
        cerr << "Vektor frag_sizes nema očekivani broj elemenata!" << endl;
        return 1;
    }

    int start = 0;
    cout << endl;
    // Ispis zaglavlja tablice
    cout << "                  ";
    for (const auto& name : names_metagenomic_sample) {
        cout << setw(15) << name;
    }
    cout << endl;

    // Ispis srednjih vrijednosti
    cout << "average_length    ";
    for (size_t i = 0; i < names_metagenomic_sample.size(); ++i) {
        double mean = calculateMean(frag_sizes, start, read_counts[i]);
        cout << setw(15) << mean;
        start += read_counts[i];
    }
    cout << endl;

    // Resetiramo početak
    start = 0;

    // Ispis standardnih devijacija
    cout << "standard_deviation";
    for (size_t i = 0; i < names_metagenomic_sample.size(); ++i) {
        double mean = calculateMean(frag_sizes, start, read_counts[i]);
        double stddev = calculateStdDev(frag_sizes, start, read_counts[i], mean);
        cout << setw(15) << stddev;
        start += read_counts[i];
    }
    cout << endl;
    // Resetiramo početak
    start = 0;

    // Ispis broja ocitanja po bakteriji
    cout << "num_of_reads      ";
    for (size_t i = 0; i < names_metagenomic_sample.size(); ++i) {
        cout << setw(15) << read_counts[i];
        start += read_counts[i];
    }
    cout << endl;
    cout << endl;

    //----opis referentnog genoma sa više vektora---------------------------

vector<reference_wrapper<vector<double>>> partial_references_to_referent_distribution_vectors;
vector<shared_ptr<vector<double>>> real_vectors3;
vector<string> names_partial_refs;  //header csv datoteke

for (size_t ref_idx = 0; ref_idx < refs.size(); ++ref_idx) {
    auto data = refs[ref_idx][0]->data_;
    int frag_size = 500000; // veličina fragmenta
    int overlap = 250000; // preklapanje

    for (uint32_t start_idx = 0; start_idx < data.size(); start_idx += overlap) {
        string fragment;
        if (start_idx + frag_size > data.size()) {
            // Ako je kraj fragmenta izvan genoma, produži na početak
            fragment = data.substr(start_idx) + data.substr(0, (start_idx + frag_size) % data.size());
        } else {
            fragment = data.substr(start_idx, frag_size);
        }

        // Simulate creating a Sequence object from the fragment
        auto seq = make_shared<Sequence>("frag", 4, fragment.c_str(), fragment.size());

        // Convert shared_ptr to unique_ptr
        auto seq_unique = make_unique<Sequence>(*seq);

        // Stvaranje shared_ptr novog vektora za svaki fragment
        shared_ptr<vector<double>> fragment_distribution_vector = make_shared<vector<double>>(get_distribution_vector(seq_unique, k, true));

        //dodajem u vanjski vektor da bude vidljivo i izvan else bloka
        real_vectors3.push_back(fragment_distribution_vector);
        // Pohrana referenci na nove vektore
        partial_references_to_referent_distribution_vectors.push_back(ref(*fragment_distribution_vector));

        // Dodavanje imena bakterije u vektor imena
        names_partial_refs.push_back(names[ref_idx] + to_string(start_idx / overlap + 1));
    }
}

//ispis u csv datoteku
string filename4 = "PARTIAL_reference_vectors.csv";
string absolute_path4 = get_data_directory() + "/" + filename4;
write_csv(partial_references_to_referent_distribution_vectors, names_partial_refs, absolute_path4);

    
    return 0;
}