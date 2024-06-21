# Završni rad

**Pristupnik**: Domagoj Sviličić (0036540224)  
**Studij**: Elektrotehnika i informacijska tehnologija i Računarstvo  
**Modul**: Računarstvo  
**Mentor**: doc. dr. sc. Krešimir Križanović  
**Zadatak**: Analiza metagenomskog uzorka na temelju distribucije k-mera  

## Opis zadatka
U bioinformatici, metagenomski uzorak sadrži genetski materijal više različitih organizama. Osnovni korak u analizi metagenomskog uzorka je određivanje organizama koji se u njemu nalaze te njihove zastupljenosti u uzorku. U analizi nizova, k-mer je podniz duljine k. Za odabranu duljinu k, za svaku sekvencu je moguće odrediti distribuciju svih k-mera koji se u toj sekvenci pojavljuju. U sklopu ovog rada potrebno je ispitati mogućnost upotrebe distribucije k-mera za analizu metagenomskog uzorka. Koristiti nekoliko genoma iz javno dostupne baze podataka RefSeq. Simulirati očitanja te generirati metagenomski uzorak. Očitanja iz metagenomskog uzorka klasiﬁcirati na temelju sličnosti distribucije k-mera između očitanja i genoma. Programski kod je potrebno komentirati i pri pisanju pratiti neki od standardnih stilova. Napisati iscrpne upute za instalaciju i korištenje. Kompletno programsko rješenje postaviti na GitHub.

## Upute za pokretanje

```bash
git clone ssh-link
#(dodaj ssh link github repozitorija)
```
U direktoriju **data** nalaze se poddirektoriji **Reads**, **Referent_genomes** i **Referent_genomes_for_NO_ERROR_reads**. U direktorij **Reads** potrebno je dodati fasta ili fastq datoteke očitanja genomskih sekvenci od kojih je potrebno kreirati metagenomski uzorak. Direktorij **Referent_genomes** potrebno je napuniti fasta datotekama referentnih genoma koje želimo koristiti za postupak klasifikacije očitanja iz metagenomskog uzorka, a direktorij **Referent_genomes_for_NO_ERROR_reads** mora sadržavati fasta datoteke referentnih genoma od kojih će biti simulirana očitanja bez pogreške. Nakon toga je potrebno ispuniti CONFIGURATION_file.txt datoteku u kojoj se postavljaju potrebni parametri analize. (Nazivi fasta/fastq datoteka koji se tamo navode moraju odgovarati stvarnim nazivima datoteka dodanima u prije spomenute direktorije)

Pokretanje programske implementacije sastoji se od navigiranja do direktorija **build** unutar projekta i izvršavanja sljedeće naredbe naredbenog redka:
```bash
cmake .. && make && cd build && ./Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24
```

U direktoriju data će se nakon pokretanja stvoriti metagenomic_sample.fasta datoteka kreiranog metagenomskog uzorka, csv datoteke i pdf izveštaj provedene analize.

