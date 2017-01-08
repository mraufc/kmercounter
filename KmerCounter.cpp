/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   kmercounter.cpp
 * Author: Mehmet Rauf Celik
 *
 * Created on September 10, 2016, 5:31 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <ctime>
#include <map>
#include <stdint.h>
#include <vector>
#include <math.h>
#include <queue>
#include <algorithm> 

#define DEFAULT_KMER_LENGTH 30
#define DEFAULT_OUTPUT_SIZE 25
#define DEFAULT_BLOOMFILTER_SIZE 838860800
#define DEFAULT_NUMBER_OF_HASH_FUNCTIONS 11

using namespace std;

map<string,uint64_t> kmercounts; //kmrecounts maps a unique K mer to its count
map<string,uint64_t>::iterator kmerit;

/*
 * The purpose of this program is to count, sort and list non-unique K-mers of
 * a DNA sequence. The program requires at least one parameter of uncompressed
 * FASTQ file name. For more details of usage, type ./kmercounter help in 
 * command line.
 * 
 * In this program, there are more than one ways to count the most
 * frequent K mers given a FASTQ file and each has its own strengths and 
 * weaknesses, however they are also useful for confirming each other's findings
 * or testing against each other as well.
 */


/* BloomFilter class is implemented to be answer to large amounts of memory
 * requirements of more straightforward solutions.
 * two good tutorials to understand what a Bloom filter does and how to
 * implement one:
 * http://billmill.org/bloomfilter-tutorial/
 * http://blog.michaelschmatz.com/2016/04/11/how-to-write-a-bloom-filter-cpp/
 */
class BloomFilter {
    /*the default size of the bloom filter, which is defined by
     * DEFAULT_BLOOMFILTER_SIZE parameter is arbitrarily chosen.
     * the size of the filter is a tradeoff between memory and the overhead
     * brought by having too many false positive returns.
     * the default number of hash functions to be used is arbitrarily chosen as
     * well. Number of hash functions whose default is indicated by
     * DEFAULT_BLOOMFILTER_SIZE is a tradeoff between computational complexity
     * and bad filter distributiion.
     */
public:
    BloomFilter(uint64_t size, uint8_t numHashes) : bfbits(size), numBfHashes(numHashes) {
        cout<<"bfbits size is "<<size<<endl;
        cout<<"bfbits capacity is "<<bfbits.capacity()<<endl;
    };
    
    /*two very simple hash functions to generate necessary amount of hash functions
     * via double hashing.
     * hash functions can be improved in the future for better performance and/or
     * distribution
     */
    
    uint64_t hash1 (const char *data, size_t len) {
        uint64_t retVal=0;
        for(size_t i=0;i<len;i++) {
            switch(data[i]) {
                case 'A':
                    retVal=retVal+3*pow(4,i);
                    break;
                case 'G':
                    retVal=retVal+2*pow(4,i);
                    break;
                case 'C':
                    retVal=retVal+1*pow(4,i);
                    break;
                    /*Tymine is not checked, T is replaced by U in RNA sequences 
                     * and also there are small number of other values such as 
                     * N for Any, B for not-A anyway*/
                default:
                    break;
            }
        }
        return retVal;
    }
    uint64_t hash2 (const char *data, size_t len) {
        uint64_t retVal=0;
        for(size_t i=0;i<len;i++) {
            retVal=retVal*101+data[i];
        }      
        return retVal;
    }
    uint64_t hashf (char n, uint64_t hash1val, uint64_t hash2val, uint64_t filterSize) {
        return (hash1val + n * hash2val) % filterSize;
    }
    
    void add(const char *data, size_t len) {
        uint64_t hash1Val=hash1(data, len);
        uint64_t hash2Val=hash2(data, len);
        for (int i = 0; i < numBfHashes; i++) {
            bfbits[hashf(i, hash1Val, hash2Val, bfbits.size())] = true;
        }
    }
    /*a true value return from contains function does not necessarily
     * mean that the filter actually contains the value or the key has 
     * been seen before. Bloom filter guarantees the lack of false negatives
     * but does not guarantee the lack of false positives.
     */
    bool contains(const char *data, size_t len) {

        uint64_t hash1Val=hash1(data, len);
        uint64_t hash2Val=hash2(data, len);
        for (int i = 0; i < numBfHashes; i++) {
            if(bfbits[hashf(i, hash1Val, hash2Val, bfbits.size())] ==false)
                return false;
        }
        return true;
    }
private:
    vector<bool> bfbits;
    uint8_t numBfHashes;
    
};

/*the purpose of the listTopElements function is to sort and display top n
 * elements of the std::map kmercounts
 */

void listTopElements(uint32_t ols) {
    cout<<"listTopElements is called"<<endl;
    uint32_t listsize=ols;
    if(listsize>kmercounts.size())
        listsize=kmercounts.size();
    cout<<"listsize is "<<listsize<<endl;
    string kmerarray[listsize];
    uint64_t kmercountsarray[listsize];
    int ix=0;
    
    kmerit = kmercounts.begin();
    for (kmerit=kmercounts.begin(); kmerit!=kmercounts.end(); ++kmerit) {
        if(ix>=listsize)
            break;
        kmerarray[ix]=kmerit->first;
        kmercountsarray[ix]=kmerit->second;
        ix++;
    }
    int i=0;
    int j=0;
    int jmax;
    string tmpstr;
    uint64_t tmpval;
    for(int i=0;i<listsize;i++) {
        jmax=i;
        for(j=i+1;j<listsize;j++) {
            if(kmercountsarray[j]>kmercountsarray[jmax]) {
                jmax=j;
            }
        }   
        if(jmax!=i) {
            tmpval=kmercountsarray[i];
            kmercountsarray[i]=kmercountsarray[jmax];
            kmercountsarray[jmax]=tmpval; 
            tmpstr=kmerarray[i];
            kmerarray[i]=kmerarray[jmax];
            kmerarray[jmax]=tmpstr;   
        }
    }
    int steps;
    for (; kmerit!=kmercounts.end(); ++kmerit) {
        steps=0;
        tmpstr=kmerit->first;
        tmpval=kmerit->second;
        for (i=listsize-1;i>=0;i--) {
            if(tmpval<=kmercountsarray[i]) 
                break;
            steps++;
        }
        if (steps>1) {
            for(j=listsize-1;j>listsize-steps;j--) {
                kmercountsarray[j]=kmercountsarray[j-1];
                kmerarray[j]=kmerarray[j-1];
            }
            kmercountsarray[listsize-steps]=tmpval;
            kmerarray[listsize-steps]=tmpstr;
        }
    }
    for(int i=0;i<listsize;i++) {
        cout<<i+1<<". "<<kmerarray[i]<<": "<<kmercountsarray[i]<<endl;
    }
}

/*after the Bloom filter phase where we can't be certain if some of the non-
 unique values were labeled as unique, we need to open the file again and 
 * do a real count 
 */
void finalizeKmerCount(string fn, uint32_t kl, uint32_t ols) {
    cout<<"finalizeKmerCount begins"<<endl;  
    ifstream file (fn.c_str());
    string line;
    string identifier;
    
    uint32_t lineix,i;
    if(!file.is_open()) {
        cout<<"error: unable to open file or file does not exist "<<fn<<endl;
        return;
    }
    
    lineix=0;
    if(getline(file,identifier)) {
        lineix++;
        if(identifier.substr(0,1)=="@") {
            cout<<"FASTQ file open is successful, first sequence identifier is:"<<endl;
            cout<<identifier<<endl;
        } else {
            cout<<"error: missing sequence identifier in first line, is this not a FASTQ file?"<<endl;
            return;
        }  
    }
    
    clock_t begin = clock();

    /*for every K mer key, value is 0, increment this value
     * whenever the key is seen in file
     */
    while (getline(file,line)) {
        lineix++;
        if(lineix%4==2 && kl<=line.length()+1) {
            for(i=0;i<line.length()-kl+1;i++) {
                kmerit=kmercounts.find(line.substr(i,kl));
                if(kmerit != kmercounts.end())
                    kmercounts[line.substr(i,kl)]++;
            }
        }
    }
    /*after the file read ends, find and remove the ones with 
     * value equal to 1 as they are false positives of Bloom filter's
     contains function*/
    for (kmerit=kmercounts.begin(); kmerit!=kmercounts.end(); ++kmerit) {
        if(kmerit->second==1)
            kmercounts.erase(kmerit);           
    }
    
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout<<"finalizeKmerCount ends"<<endl;
    cout<<"elapsed time is "<<elapsed_secs<<endl;

    std::priority_queue<int, std::vector<int>, std::greater<int> > q2;
    file.close();
    
    listTopElements(ols);
    
}

/*this function counts using a Bloom filter for better memory efficiency
 choose this option for larger K values unless the file size is small.
 */
void kmerCountBloomFilter(string fn, uint32_t kl, uint32_t ols) {
    cout<<"kmerCountBloomFilter begins"<<endl;
    cout<<"target FASTQ file name is "<<fn<<endl;
    cout<<"k-mer length is "<<kl<<endl;
    cout<<"output size is "<<ols<<endl;
    
    ifstream file (fn.c_str());
    string line;
    string identifier;
    
    uint64_t bfsize=DEFAULT_BLOOMFILTER_SIZE;
    
    BloomFilter bf = BloomFilter(bfsize,DEFAULT_NUMBER_OF_HASH_FUNCTIONS);
    
    uint32_t lineix,i;
    if(!file.is_open()) {
        cout<<"error: unable to open file or file does not exist "<<fn<<endl;
        return;
    }
    
    lineix=0;
    if(getline(file,identifier)) {
        lineix++;
        if(identifier.substr(0,1)=="@") {
            cout<<"FASTQ file open is successful, first sequence identifier is:"<<endl;
            cout<<identifier<<endl;
        } else {
            cout<<"error: missing sequence identifier in first line, is this not a FASTQ file?"<<endl;
            return;
        }  
    }
    uint64_t hitCount=0;
    clock_t begin = clock();

    while (getline(file,line)) {
        lineix++;
        if(lineix%4==2 && kl<=line.length()+1) {
            for(i=0;i<line.length()-kl+1;i++) {
                /*if Bloom filter contains the K mer key (or if it thinks it does)
                 * then add it to possible unique values map. 
                 * if Bloom filter does not contain the K mer key, then insert the key
                 */
                if(bf.contains(line.substr(i,kl).c_str(),kl)) {
                    hitCount++;
                    /*there may be false positives so insert K mer key with 0 value.
                     *Bloom filter guarantees that there will not be any false negatives,
                     * so no need to check for insert fail
                     */
                    kmercounts.insert(pair<string,uint64_t>(line.substr(i,kl),0));
                } else {
                    bf.add(line.substr(i,kl).c_str(),kl);
                }
            }
        }
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout<<"kmerCountBloomFilter ends"<<endl;
    cout<<"elapsed time is "<<elapsed_secs<<endl;
    file.close();
    finalizeKmerCount(fn, kl, ols);
}
/*this is a very straightforward method that uses std::map
 * this method reads the file only once.
 * use this function when K is a small value no matter what input file size is
 * as the other methods use some sort of buffer (one of them uses a Bloom filter 
 * and the other uses a set) as these buffers become useless when input variance
 * is low.
 */
void kmerCountMapOnly(string fn, uint32_t kl, uint32_t ols) {
    cout<<"kmerCountMapOnly begins"<<endl;
    cout<<"target FASTQ file name is "<<fn<<endl;
    cout<<"k-mer length is "<<kl<<endl;
    cout<<"output size is "<<ols<<endl;
    
    ifstream file (fn.c_str());
    string line;
    string identifier;

    uint32_t lineix,i;
    if(!file.is_open()) {
        cout<<"error: unable to open file or file does not exist "<<fn<<endl;
        return;
    }
    
    lineix=0;
    if(getline(file,identifier)) {
        lineix++;
        if(identifier.substr(0,1)=="@") {
            cout<<"FASTQ file open is successful, first sequence identifier is:"<<endl;
            cout<<identifier<<endl;
        } else {
            cout<<"error: missing sequence identifier in first line, is this not a FASTQ file?"<<endl;
            return;
        }  
    }
    
    clock_t begin = clock();

    while (getline(file,line)) {
        lineix++;
        if(lineix%4==2 && kl<=line.length()+1) {
            for(i=0;i<line.length()-kl+1;i++) {
                if(!kmercounts.insert(pair<string,uint64_t>(line.substr(i,kl),0)).second) {
                    if(kmercounts[line.substr(i,kl)]==0)
                        kmercounts[line.substr(i,kl)]=2;
                    else
                        kmercounts[line.substr(i,kl)]++;
                }
            }
        }
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout<<"kmerCountMapOnly ends"<<endl;
    cout<<"elapsed time is "<<elapsed_secs<<endl;

    file.close();
    listTopElements(ols);
    
}


/**this is a middle of the road counting method that shines when input variance
 * is neither high nor low (i.e K=10) and file size is not large enough for 
 * enough K mer variance to fill up the memory
 **/
void kmerCountSetMap(string fn, uint32_t kl, uint32_t ols) {
    cout<<"kmerCountSetMap begins"<<endl;
    cout<<"target FASTQ file name is "<<fn<<endl;
    cout<<"k-mer length is "<<kl<<endl;
    cout<<"output size is "<<ols<<endl;
    
    ifstream file (fn.c_str());
    string line;
    string identifier;

    set<string> distinctkmer;


    uint32_t lineix,i;
    if(!file.is_open()) {
        cout<<"error: unable to open file or file does not exist "<<fn<<endl;
        return;
    }
    
    lineix=0;
    if(getline(file,identifier)) {
        lineix++;
        if(identifier.substr(0,1)=="@") {
            cout<<"FASTQ file open is successful, first sequence identifier is:"<<endl;
            cout<<identifier<<endl;
        } else {
            cout<<"error: missing sequence identifier in first line, is this not a FASTQ file?"<<endl;
            return;
        }  
    }
    
    clock_t begin = clock();
    
    while (getline(file,line)) {
        lineix++;
        if(lineix%4==2 && kl<=line.length()+1) {
            for(i=0;i<line.length()-kl+1;i++) {
                if(!distinctkmer.insert(line.substr(i,kl)).second) {
                    if (!kmercounts.insert(pair<string,long long>(line.substr(i,kl),2)).second) {
                        kmercounts[line.substr(i,kl)]++;
                    } 
                }
            }
        }
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout<<"kmerCountSetMap begins"<<endl;
    cout<<"elapsed time is "<<elapsed_secs<<endl;

    file.close();
    listTopElements(ols);
    
}
int main(int argc, char** argv) {    
    string filename="";
    uint32_t kmerlength=DEFAULT_KMER_LENGTH;
    uint32_t outpsize=DEFAULT_OUTPUT_SIZE;
    string implementation="bf";
    bool isParseOK=true;
    if(argc==1 || (argc==2 && string(argv[1])=="help") || (argc==2 && string(argv[1])=="?")) {
        cout<<"kmercounter help:"<<endl;
        cout<<"./kmercounter --f [filename] --k [k-mer length] --l [output list size] --i [implementation]"<<endl; 
        cout<<"a file name following option --f must be specified"<<endl;
        cout<<"default k-mer length is "<<DEFAULT_KMER_LENGTH<<endl;
        cout<<"default output list size is "<<DEFAULT_OUTPUT_SIZE<<endl;
        cout<<"default implementation is bf for Bloom Filter"<<endl;
        cout<<"other valid options are:"<<endl;
        cout<<"sm: a two step set map based implementation"<<endl;
        cout<<"mo: a 1 step map only implementation"<<endl;
        cout<<"general guidelines for choosing an impelentation:"<<endl;
        cout<<"for any K larger than 15, choose bloom filter unless the input file size is very small (<=50mb)"<<endl;
        cout<<"for any K smaller than 5, map only implementation is generally better even if the file size is very large"<<endl;
        cout<<"in between above K values, two step set map implementation pulls ahead"<<endl;
        cout<<"default Bloom filter size is "<<DEFAULT_BLOOMFILTER_SIZE<<" bits"<<endl;
        cout<<"default number of hash functions is "<<DEFAULT_NUMBER_OF_HASH_FUNCTIONS<<endl;
        
        
        isParseOK=false;
    } else if (argc==3 | argc==5 | argc==7 | argc==9) {
        for(uint32_t i=1;i<argc-1;i=i+2) {
            if(string(argv[i])=="--f") {
                filename=string(argv[i+1]);
                if(filename=="") {
                    cout<<"error: file name ./kmercounter help"<<endl;
                    isParseOK=false;
                }
            } else if(string(argv[i])=="--k") {
                char *strend;
                kmerlength = strtol(string(argv[i+1]).c_str(), &strend, 10);
                if(strend==string(argv[i+1])) {
                    cout<<"error: non-numeric k-mer length, try ./kmercounter help"<<endl;
                    isParseOK=false;
                }
                else if(kmerlength<=0) {
                    cout<<"error: non-positive k-mer length value, try ./kmercounter help"<<endl;
                    isParseOK=false;
                }
            } else if(string(argv[i])=="--l") {
                char *strend;
                outpsize = strtol(string(argv[i+1]).c_str(), &strend, 10);
                if(strend==string(argv[i+1])) {
                    cout<<"error: non-numeric output list size, try ./kmercounter help"<<endl;
                    isParseOK=false;
                }
                else if(outpsize<=0) {
                    cout<<"error: non-positive output list size value, try ./kmercounter help"<<endl;
                    isParseOK=false;
                }
            } else if(string(argv[i])=="--i") {
                implementation=string(argv[i+1]);
                if(implementation!="bf" && implementation != "sm" && implementation !="mo") {
                    cout<<"error: incorrect implementation value, try ./kmercounter help"<<endl;
                    isParseOK=false;
                }
            } else if(string(argv[i]).substr(0,2)=="--") {
                cout<<"error: unknown option "<<string(argv[i])<<", try ./kmercounter help"<<endl;
                isParseOK=false;
            }
        }
        if(filename=="") {
            cout<<"error: a file name following option --f must be specified"<<endl;
            isParseOK=false;
        }
    } else {
        cout<<"error: number of arguments is incorrect, try ./kmercounter help"<<endl;
        isParseOK=false;
    }
    if(isParseOK) {
        if(implementation=="bf")
            kmerCountBloomFilter(filename,kmerlength,outpsize);
        else if (implementation=="sm")
            kmerCountSetMap(filename,kmerlength,outpsize);
        else if (implementation=="mo")
            kmerCountMapOnly(filename,kmerlength,outpsize);
    }
        
    
    return 0;
}


